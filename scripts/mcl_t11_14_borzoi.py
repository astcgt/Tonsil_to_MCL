#!/usr/bin/env python
"""
Plan B for AlphaGenome: Borzoi (via gReLU) for t(11;14) enhancer hijacking in MCL.

NOT EXECUTED -- API matched line-by-line against gReLU 1.1.0 sdist + official docs:
  grelu.resources.load_model            https://genentech.github.io/gReLU/autoapi/grelu/resources/index.html
  LightningModel.predict_on_seqs        https://genentech.github.io/gReLU/autoapi/grelu/lightning/index.html
  LightningModel.input_coord_to_output_bin   (same page)
  grelu.variant.predict_variant_effects https://genentech.github.io/gReLU/autoapi/grelu/variant/index.html
  ckpt + task metadata                  https://huggingface.co/Genentech/borzoi-model
  geometry (crop_len 5120 -> 6144 bins) https://github.com/calico/borzoi/blob/main/examples/params.json

Geometry: input 524,288 bp -> output 6,144 bins x 32 bp = 196,608 bp, centred.
"""

import numpy as np
import pandas as pd
import torch

import grelu.resources
from grelu.sequence.format import intervals_to_strings

# ---------------------------------------------------------------- config ----
# "hg38" lets genomepy download it; a local path needs a .fai (samtools faidx).
GENOME  = "hg38"
SEQ_LEN = 524_288
HALF    = SEQ_LEN // 2
DEVICE  = 0

# CCND1, GRCh38, + strand. Verified: Ensembl REST lookup/symbol/homo_sapiens/CCND1
CCND1_CHROM, CCND1_TSS, CCND1_END = "chr11", 69_641_070, 69_654_474

# chr11 breakpoint, 5' (centromeric) of CCND1. ~120 kb is a weakly sourced
# literature ballpark for the major translocation cluster. USE YOUR OWN CALL.
BP_UPSTREAM = 120_000

# IGH, chr14, minus strand. Verified anchors (Ensembl REST):
#   IGHM 105,851,705-105,856,218   IGHJ6 105,863,198-105,863,258
# Emu sits between them, so Emu ~ chr14:105,856,218-105,863,198.
IGH_CHROM, IGH_BP = "chr14", 105_863_300
EMU_LO, EMU_HI    = 105_856_218, 105_863_198

# der(14) = chr14 pter..JH joined to chr11 breakpoint..qter, both forward on the
# chr14/chr11 plus strand, which puts Emu immediately 5' of the junction and
# CCND1 downstream in its normal + orientation, giving Emu -> TSS ~120-127 kb.
# "revcomp" is a NEGATIVE CONTROL: it flips the block so Emu lands at the far
# end (offset 102-7,082), i.e. ~255-262 kb from the TSS, with the same base
# composition. Effect should largely vanish.
IGH_ORIENT = "forward"      # "forward" (der(14)) | "revcomp" (control)


def revcomp(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def fetch(chrom, start, end):
    iv = pd.DataFrame({"chrom": [chrom], "start": [start], "end": [end]})
    return intervals_to_strings(iv, genome=GENOME)[0].upper()


# ------------------------------------------- 1. build REF / ALT windows ----
# CCND1 TSS sits at the window centre in BOTH, so output bins align.
win_start = CCND1_TSS - HALF
win_end   = win_start + SEQ_LEN

ref_seq = fetch(CCND1_CHROM, win_start, win_end)

keep_chr11 = BP_UPSTREAM + HALF        # chr11 bases kept: junction -> window end
n_igh      = SEQ_LEN - keep_chr11      # IGH bases spliced in (142,144 here)

chr11_part = fetch(CCND1_CHROM, CCND1_TSS - BP_UPSTREAM, win_end)
igh_core   = fetch(IGH_CHROM, IGH_BP - n_igh, IGH_BP)   # Emu at its 3' end
igh_part   = igh_core if IGH_ORIENT == "forward" else revcomp(igh_core)

alt_seq = igh_part + chr11_part

assert len(ref_seq) == len(alt_seq) == SEQ_LEN
assert ref_seq[-keep_chr11:] == alt_seq[-keep_chr11:], "CCND1 side must match"

# Where did Emu land? Junction is at window offset n_igh.
igh_seg_start = IGH_BP - n_igh
if IGH_ORIENT == "forward":
    emu_off = (EMU_LO - igh_seg_start, EMU_HI - igh_seg_start)
else:
    emu_off = (n_igh - (EMU_HI - igh_seg_start), n_igh - (EMU_LO - igh_seg_start))
print(f"junction at offset {n_igh:,}; Emu at {emu_off[0]:,}-{emu_off[1]:,}")
print(f"Emu -> CCND1 TSS: {HALF - emu_off[1]:,} to {HALF - emu_off[0]:,} bp")

# ------------------------------------------------------- 2. load Borzoi ----
model = grelu.resources.load_model(
    repo_id="Genentech/borzoi-model", filename="human_rep0.ckpt",
)
model.eval()

tasks = pd.DataFrame(model.data_params["tasks"])     # name/assay/sample/scale/...
print(tasks.columns.tolist(), len(tasks))

# ---------------------------------------- 3. select lymphoid / B tracks ----
B_PAT = (r"GM12878|B.?lymphoblastoid|CD19\+|Memory B|Plasma B|B Lymphocyte"
         r"|spleen|tonsil|lymph node")
is_b = tasks["sample"].str.contains(B_PAT, case=False, regex=True, na=False)
# CCND1 is on +, so drop antisense tracks (identifier ending in '-').
is_sense = ~tasks["name"].astype(str).str.endswith("-")

groups = {a: tasks.index[is_b & is_sense & (tasks["assay"] == a)].to_numpy()
          for a in ("RNA", "CAGE", "ATAC", "DNASE")}
for a, idx in groups.items():
    print(f"  {a:6s}: {len(idx)} tracks")

# ---------------------------------------------------------- 4. predict ----
# One sequence at a time: a 524 kb forward pass is memory-hungry and batching 2
# roughly doubles peak activations. Calico's own example notebook states the
# gradient step is "runnable only on 40GB GPU cards, e.g. A100".
with torch.no_grad():
    ref_pred = model.predict_on_seqs(ref_seq, device=DEVICE)[0]   # (7611, 6144)
    alt_pred = model.predict_on_seqs(alt_seq, device=DEVICE)[0]
print("pred shape:", ref_pred.shape,
      f"| peak VRAM {torch.cuda.max_memory_allocated(DEVICE)/2**30:.1f} GiB")


def unsquash(x, clip_soft=384.0, track_scale=0.01, track_transform=0.75):
    """Reverse Borzoi's squashed (sqrt-ish) output scale.

    Borzoi targets use sum_stat=sum_sqrt with a soft clip, so raw outputs are
    NOT linear coverage and a raw ratio is not a true fold change. The transform
    is monotonic, so the DIRECTION of a WT-vs-ALT effect is preserved either way.
    Values follow Calico's example notebook (untransform_old=True).
    """
    x = np.asarray(x, dtype=np.float64) / track_scale
    m = x > clip_soft
    x[m] = (x[m] - clip_soft) ** 2 + clip_soft
    return x ** (1.0 / track_transform)


# Bin ranges. crop_len=5120, bin_size=32 -> TSS lands on centre bin 3072.
tss_bin  = model.input_coord_to_output_bin(CCND1_TSS, start_pos=win_start)
end_bin  = model.input_coord_to_output_bin(CCND1_END, start_pos=win_start)
prom_lo, prom_hi = tss_bin - 16, tss_bin + 16      # ~1 kb around the TSS
print(f"CCND1 bins: TSS {tss_bin}, end {end_bin} (of 6144)")

for assay, idx in groups.items():
    if len(idx) == 0:
        continue
    # RNA over the gene body; CAGE/ATAC/DNASE over the promoter.
    lo, hi = (tss_bin, end_bin + 1) if assay == "RNA" else (prom_lo, prom_hi)
    r = unsquash(ref_pred[idx, lo:hi]).sum(axis=1)
    a = unsquash(alt_pred[idx, lo:hi]).sum(axis=1)
    lfc = np.log2((a + 1.0) / (r + 1.0))
    out = pd.DataFrame({
        "track":  tasks.loc[idx, "name"].to_numpy(),
        "sample": tasks.loc[idx, "sample"].to_numpy(),
        "ref": r, "alt": a, "log2FC": lfc,
    }).sort_values("log2FC", ascending=False)
    print(f"\n===== {assay} at CCND1: WT vs t(11;14) =====")
    print(out.to_string(index=False))
    print(f"median log2FC = {np.median(lfc):+.3f}")
