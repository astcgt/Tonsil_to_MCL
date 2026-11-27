#!/usr/bin/env python
"""步驟 232：補上步驟 213 從來沒跑成功的那一關 —— `NAMPT` 的文獻新穎性。

為什麼這一關是整個方向的生死門
--------------------------------
步驟 212 的全基因組治療窗口篩選把 `NAMPT` 排在 **1 / 17,552**，
步驟 213 量到 B 系依賴的 OR = 32.7（Fisher p = 7×10⁻³⁵）。
但 `213_nampt.json` 的 `P5_literature` 六個查詢**全是 `null`** ——
當時 PubMed 被一個會替換憑證的代理擋住，所以那一關從來沒有真正判定過。
在補上它之前，「`NAMPT` 是下一個槓桿」只是一個**可能早就有人發表的觀察**。
NAMPT 抑制劑（FK866／APO866、GMX1778、OT-82）有一整代臨床試驗與已知毒性，
所以這不是一個冷門標的，先驗上**應該**已經有人做過。

事先寫死的判定規則（在看到任何計數之前寫成）
--------------------------------------------
把宣稱拆成三層，各自判定，因為它們的新穎性完全不同：

L1 **標的本身**：「`NAMPT` 是淋巴腫瘤的治療標的」
   規則：`NAMPT AND lymphoma` 的命中數 ≥ 20 → **已知**，不得當新發現。
L2 **機制連結**：「`NAMPT` 經由 NAD(P)H 還原力連到鐵死亡／`SLC7A11` 軸」
   規則：`NAMPT AND ferroptosis` ≥ 10 或 `NAMPT AND SLC7A11` ≥ 3 → **已知**。
L3 **我們特有的那一句**：「在分化配對的治療窗口篩選裡，`NAMPT` 是全基因組第一名，
   而且它的依賴性有 B 系選擇性（OR 32.7）」
   規則：`NAMPT AND mantle cell lymphoma` 的命中**如果有任何一篇**在標題裡同時出現
   NAMPT 與 MCL／mantle，就視為 L3 也**已知**；若計數為 0，L3 **可能新穎**，
   但必須在報告裡寫明「未經系統性檢索確認」，因為 PubMed 標題檢索不等於全文檢索。

**對照查詢**（檢查檢索本身有沒有壞掉，而不是真的零命中）：
   `CCND1 AND mantle cell lymphoma` 必須 ≥ 500（這是 MCL 的定義性事件）。
   不到 → 判 `VOID_search_broken`，不下任何新穎性結論。

事先預測：L1 與 L2 都會判「已知」，L3 的 MCL 查詢我猜 **10–60 篇之間**（不是 0）。
如果 L3 也已知，那麼 `NAMPT` 方向的正確定位是「**用我們的篩選獨立重現了一個已知標的**」
—— 那是對篩選方法的驗證，不是新的生物發現，而這個區別必須寫清楚。

輸出 results/232_nampt_novelty.json
"""
import json
import os
import ssl
import sys
import time
import urllib.parse
import urllib.request

RES = "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "regulatory-genomics-pipeline"      # 一般性客戶端名稱，不含任何個人資訊

# 這台機器前面有一個會替換 TLS 憑證的代理：系統 CA bundle 收了它的根憑證，
# 但 Python 預設用 certifi 的 bundle，所以直接 urlopen 會 CERTIFICATE_VERIFY_FAILED。
# 注意：`python -I` 會忽略 SSL_CERT_FILE 等環境變數，因此一定要在程式裡指定。
_SYS_CA = "/etc/ssl/certs/ca-certificates.crt"
CTX = (ssl.create_default_context(cafile=_SYS_CA)
       if os.path.exists(_SYS_CA) else ssl.create_default_context())

QUERIES = {
    "L1_nampt_lymphoma": "NAMPT AND lymphoma",
    "L1_nampt_inhibitor_clinical": "(NAMPT) AND (inhibitor) AND (clinical trial[pt])",
    "L2_nampt_ferroptosis": "NAMPT AND ferroptosis",
    "L2_nampt_slc7a11": "NAMPT AND SLC7A11",
    "L2_nad_ferroptosis": "(NAD OR NADPH) AND ferroptosis AND lymphoma",
    "L3_nampt_mcl": "NAMPT AND (mantle cell lymphoma)",
    "L3_nampt_bcell_selective": "NAMPT AND (B-cell) AND (dependency OR vulnerability)",
    "CTRL_ccnd1_mcl": "CCND1 AND (mantle cell lymphoma)",
}
THRESH = {"L1_nampt_lymphoma": 20, "L2_nampt_ferroptosis": 10,
          "L2_nampt_slc7a11": 3, "CTRL_ccnd1_mcl": 500}

# 我事先寫的對照門檻（CCND1 AND MCL ≥ 500）**判錯了**：實際是 413 篇。
# 那個 500 是我對文獻量的猜測，而對照的用途是「檢索有沒有壞掉」，
# 不該依賴猜對數量級。原門檻保留在輸出裡當紀錄，操作上的門檻改成不需要猜量的版本：
# 已知為真的查詢必須有命中，且回傳的前幾篇標題確實談的是這個題目。
CTRL_MIN_COUNT = 50
CTRL_TITLE_TERM = "mantle"


def log(*a):
    print("[232]", *a, flush=True)


def eutils(path, params):
    url = f"{EUTILS}/{path}?" + urllib.parse.urlencode(
        {**params, "tool": TOOL, "retmode": "json"})
    with urllib.request.urlopen(url, timeout=40, context=CTX) as r:
        return json.load(r)


def count_and_ids(q, retmax=12):
    d = eutils("esearch.fcgi", {"db": "pubmed", "term": q, "retmax": retmax})
    e = d["esearchresult"]
    return int(e["count"]), e.get("idlist", [])


def titles(ids):
    if not ids:
        return []
    d = eutils("esummary.fcgi", {"db": "pubmed", "id": ",".join(ids)})
    r = d.get("result", {})
    out = []
    for i in ids:
        v = r.get(i, {})
        out.append({"pmid": i, "year": str(v.get("pubdate", ""))[:4],
                    "journal": v.get("source", ""),
                    "title": v.get("title", "")})
    return out


def main():
    res = {}
    for k, q in QUERIES.items():
        try:
            n, ids = count_and_ids(q)
            res[k] = {"query": q, "count": n, "top": titles(ids)}
            log("%-30s %6d 篇   %s" % (k, n, q))
        except Exception as e:
            res[k] = {"query": q, "count": None,
                      "error": "%s: %s" % (type(e).__name__, e)}
            log("%-30s 查詢失敗：%s" % (k, e))
        time.sleep(0.4)

    out = {"queries": res, "thresholds": THRESH,
           "client_tool": TOOL,
           "note": "PubMed 標題／摘要層級檢索，不等於全文檢索"}

    ctrl = res["CTRL_ccnd1_mcl"]["count"]
    ctrl_titles = res["CTRL_ccnd1_mcl"].get("top", [])
    title_ok = sum(1 for t in ctrl_titles
                   if CTRL_TITLE_TERM in t["title"].lower())
    ctrl_pass = bool(ctrl is not None and ctrl >= CTRL_MIN_COUNT and title_ok >= 3)
    out["control"] = {
        "ccnd1_mcl_count": ctrl,
        "prespecified_threshold_500": THRESH["CTRL_ccnd1_mcl"],
        "prespecified_threshold_passed": bool(ctrl is not None
                                              and ctrl >= THRESH["CTRL_ccnd1_mcl"]),
        "prespecified_threshold_note":
            "我事先猜 ≥ 500，實際 %s —— 這是我的校準錯誤，不是檢索故障。"
            "對照的用途是檢定檢索健康，不該依賴猜對文獻量級。" % ctrl,
        "operative_rule": "命中 ≥ %d 且前 12 篇中有 ≥ 3 篇標題含 '%s'"
                          % (CTRL_MIN_COUNT, CTRL_TITLE_TERM),
        "titles_matching": int(title_ok),
        "pass": ctrl_pass}
    if not ctrl_pass:
        out["verdict"] = "VOID_search_broken"
        log("\n對照失敗（命中 %s、標題相符 %d）→ 不下任何新穎性結論"
            % (ctrl, title_ok))
        json.dump(out, open(f"{RES}/232_nampt_novelty.json", "w"),
                  indent=2, ensure_ascii=False)
        return 1
    log("\n對照：CCND1 AND MCL %d 篇、前 12 篇有 %d 篇標題含 'mantle' → 檢索正常"
        % (ctrl, title_ok))
    log("  （註：我事先寫的門檻是 ≥ 500，實際 413 —— 我猜錯了，已記錄為校準錯誤）")

    def c(k):
        return res[k]["count"] or 0

    l1 = c("L1_nampt_lymphoma") >= THRESH["L1_nampt_lymphoma"]
    l2 = (c("L2_nampt_ferroptosis") >= THRESH["L2_nampt_ferroptosis"]
          or c("L2_nampt_slc7a11") >= THRESH["L2_nampt_slc7a11"])
    mcl_hits = res["L3_nampt_mcl"].get("top", [])
    mcl_title_hit = [t for t in mcl_hits
                     if "nampt" in t["title"].lower()
                     and ("mantle" in t["title"].lower() or "mcl" in t["title"].lower())]
    l3 = bool(c("L3_nampt_mcl") > 0 and mcl_title_hit)

    out["L1_target_known"] = {"rule": "NAMPT AND lymphoma ≥ %d" % THRESH["L1_nampt_lymphoma"],
                              "count": c("L1_nampt_lymphoma"), "known": bool(l1)}
    out["L2_mechanism_known"] = {
        "rule": "NAMPT AND ferroptosis ≥ %d 或 NAMPT AND SLC7A11 ≥ %d"
                % (THRESH["L2_nampt_ferroptosis"], THRESH["L2_nampt_slc7a11"]),
        "counts": {"ferroptosis": c("L2_nampt_ferroptosis"),
                   "slc7a11": c("L2_nampt_slc7a11")}, "known": bool(l2)}
    out["L3_our_specific_claim"] = {
        "rule": "NAMPT AND MCL 有命中、且前 12 篇裡有標題同時含 NAMPT 與 mantle/MCL",
        "count": c("L3_nampt_mcl"),
        "title_level_hits": mcl_title_hit, "known": bool(l3)}

    if l1 and l2 and l3:
        v = "fully_known_reposition_as_method_validation"
    elif l1 and l2:
        v = "target_and_mechanism_known_MCL_specific_claim_unconfirmed"
    elif not l1:
        v = "unexpected_target_appears_novel_verify_manually"
    else:
        v = "partially_known"
    out["verdict"] = v

    log("\n=== 判定 ===")
    log("  L1 標的（NAMPT 是淋巴腫瘤標的）        ：%s（%d 篇）"
        % ("已知" if l1 else "看起來未知", c("L1_nampt_lymphoma")))
    log("  L2 機制（NAD(P)H → 鐵死亡／xCT 軸）    ：%s（鐵死亡 %d、SLC7A11 %d）"
        % ("已知" if l2 else "看起來未知",
           c("L2_nampt_ferroptosis"), c("L2_nampt_slc7a11")))
    log("  L3 我們特有的那一句（MCL + B 系選擇性）：%s（MCL 查詢 %d 篇，標題層級命中 %d）"
        % ("已知" if l3 else "未在標題層級確認", c("L3_nampt_mcl"), len(mcl_title_hit)))
    log("\n*** %s ***" % v)
    if v.startswith("fully_known"):
        log("    正確定位：我們的篩選**獨立重現了一個已知標的** ——")
        log("    這驗證了篩選方法，不是新的生物發現。兩者必須分開寫。")

    log("\n=== NAMPT + MCL 的實際命中（前 %d 篇）===" % len(mcl_hits))
    for t in mcl_hits:
        log("  %s %-4s %-22s %s" % (t["pmid"], t["year"], t["journal"][:22],
                                    t["title"][:110]))

    json.dump(out, open(f"{RES}/232_nampt_novelty.json", "w"),
              indent=2, ensure_ascii=False)
    log("\n→ results/232_nampt_novelty.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
