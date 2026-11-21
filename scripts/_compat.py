"""
相容層：scanpy 1.11.5 會對 scipy 的 *sparse array*（csr_array 等）呼叫 getnnz()，
但 scipy 已把 get* 系列方法從 sparse array API 移除（只保留在 sparse matrix）。
這裡補回最小可用的 getnnz / getformat，避免 sc.pp.neighbors 直接崩掉。

用法：在任何 import scanpy 的腳本最上方加
    import _compat  # noqa
"""
import numpy as np
import scipy.sparse as _sp


def _getnnz(self, axis=None):
    if axis is None:
        return int(self.nnz)
    if self.format == "csr":
        if axis in (1, -1):
            return np.diff(self.indptr)
        return np.bincount(self.indices, minlength=self.shape[1])
    if self.format == "csc":
        if axis in (0, -2):
            return np.diff(self.indptr)
        return np.bincount(self.indices, minlength=self.shape[0])
    return self.tocsr().getnnz(axis=axis)


def _getformat(self):
    return self.format


_patched = []
for _cls_name in ("csr_array", "csc_array", "coo_array", "lil_array",
                  "dia_array", "bsr_array", "dok_array"):
    _cls = getattr(_sp, _cls_name, None)
    if _cls is None:
        continue
    if not hasattr(_cls, "getnnz"):
        _cls.getnnz = _getnnz
        _patched.append(_cls_name)
    if not hasattr(_cls, "getformat"):
        _cls.getformat = _getformat

if _patched:
    print(f"[_compat] 已為 {', '.join(_patched)} 補上 getnnz()", flush=True)


# ---------------------------------------------------------------------------
# scanpy 的 _get_indices_distances_from_sparse_matrix 在切片後會回傳
# 非 C-contiguous 的 view，umap 的 numba kernel 只吃 C-contiguous → TypingError。
# 這裡包一層強制轉連續。
# ---------------------------------------------------------------------------
try:
    from scanpy.neighbors import _common as _sc_common

    _orig = _sc_common._get_indices_distances_from_sparse_matrix

    def _wrapped(d, n_neighbors):
        ind, dist = _orig(d, n_neighbors)
        return np.ascontiguousarray(ind), np.ascontiguousarray(dist)

    if getattr(_orig, "__name__", "") != "_wrapped":
        _sc_common._get_indices_distances_from_sparse_matrix = _wrapped
        # scanpy.neighbors 以 from-import 取用，需一併覆寫
        import scanpy.neighbors as _sc_nb
        if hasattr(_sc_nb, "_get_indices_distances_from_sparse_matrix"):
            _sc_nb._get_indices_distances_from_sparse_matrix = _wrapped
        print("[_compat] 已修補 scanpy 的 kNN 連續性問題", flush=True)
except Exception as _e:  # pragma: no cover
    print("[_compat] scanpy 修補略過:", _e, flush=True)
