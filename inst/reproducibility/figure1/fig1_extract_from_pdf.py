"""Recover the plotted data of Figure 1 (chessboard toy example) from the
original R-generated vector PDFs.  No data are simulated: every value is read
from the PDF drawing operators (filled rectangles for image(), polylines for
the u/v strips).

Outputs (CSV):
  fig1_heatmaps_from_pdf.csv : image, row (1 = top), col (1 = left), gray (0=black..1=white)
  fig1_vectors_from_pdf.csv  : vector (u/v), method, role (truth/estimate), index,
                               value_rel  (position inside the original strip's plot region, 0=bottom,1=top)
                               value_std  (per-strip affine map so the true vector's two levels are 0 and 1)
                               page_x, page_y (pt, top-down), stroke_rgb, lwd_pt
Usage: python fig1_extract_from_pdf.py <dir_with_pdfs_or_paths.json> <outdir>
"""
import sys, json, os, numpy as np, pandas as pd, pymupdf

def load_paths(arg):
    if arg.endswith(".json"):
        return json.load(open(arg))
    names = dict(True_="True.pdf", Noisy="NoisyChessboard.pdf", SVD="svd_Reconstructed.pdf",
                 Smooth="smooth_reconstructed.pdf", Lasso="Lasso_reconstructed.pdf",
                 Both="pic.compress.pdf", u1="u1plots.pdf", v1="v1plots.pdf")
    return {k.rstrip("_"): os.path.join(arg, v) for k, v in names.items()}

IMAGES = [("True", "True"), ("Noisy", "Noisy"), ("SVD", "SVD"), ("Smooth", "Smooth"),
          ("Lasso", "Group lasso"), ("Both", "Smooth + group lasso")]
STRIP_METHODS = ["SVD", "Smooth", "Group lasso", "Smooth + group lasso"]  # top -> bottom (from strip labels)

def extract_image(path):
    pg = pymupdf.open(path)[0]
    fills = [d for d in pg.get_drawings() if d["type"] == "f"]
    R = np.array([[d["rect"].x0, d["rect"].y0, d["rect"].x1, d["rect"].y1] for d in fills])
    C = np.array([d["fill"] for d in fills])
    assert np.allclose(C[:, 0], C[:, 1]) and np.allclose(C[:, 1], C[:, 2]), "non-gray fill"
    assert all(d.get("fill_opacity", 1) == 1 for d in fills)
    xs = np.unique(np.round(R[:, 0], 2)); ys = np.unique(np.round(R[:, 1], 2))
    ix = np.searchsorted(xs, np.round(R[:, 0], 2)); iy = np.searchsorted(ys, np.round(R[:, 1], 2))
    M = np.full((len(ys), len(xs)), np.nan)
    M[iy, ix] = C[:, 0]
    assert not np.isnan(M).any() and len(fills) == M.size, "grid not complete / duplicated cells"
    # R's pdf() writes 3-decimal sRGB; snap to the 255-level gray actually written
    M = np.round(M * 255) / 255
    frame = (R[:, 0].min(), R[:, 1].min(), R[:, 2].max(), R[:, 3].max())
    return M, frame

def extract_strips(path, vec):
    pg = pymupdf.open(path)[0]
    H = pg.rect.height
    raw = b"".join(pg.parent.xref_stream(x) for x in pg.get_contents()).decode("latin1")
    import re
    # plot-region clip rectangles of the 4 strips (PDF bottom-up coords: x y w h)
    clips = [tuple(map(float, c)) for c in re.findall(r"([\d.]+) ([\d.]+) ([\d.]+) ([\d.]+) re W n", raw)]
    regions = sorted({c for c in clips if c[0] > 0}, key=lambda c: -c[1])  # top strip first
    assert len(regions) == 4, regions
    strokes = [d for d in pg.get_drawings() if d["type"] == "s"]
    assert len(strokes) == 8
    rows = []
    for j, d in enumerate(strokes):
        pts = [d["items"][0][1]] + [it[2] for it in d["items"]]
        assert all(it[0] == "l" for it in d["items"])
        xy = np.array([[p.x, p.y] for p in pts]); n = len(xy)
        # assign to strip by vertical position
        ymid = xy[:, 1].mean()
        s = [k for k, (x0, y0, w, h) in enumerate(regions) if H - (y0 + h) - 0.05 <= ymid <= H - y0 + 0.05][0]
        x0, y0, w, h = regions[s]; ytop = H - (y0 + h)
        gray = np.allclose(d["color"][0], d["color"][1]) and np.allclose(d["color"][1], d["color"][2])
        role = "truth" if gray else "estimate"
        idx = 1 + (xy[:, 0] - x0) / w * (n - 1)
        assert np.abs(idx - np.round(idx)).max() < 0.01, "x is not an integer index grid"
        rel = (ytop + h - xy[:, 1]) / h
        rgb = "#%02X%02X%02X" % tuple(int(round(c * 255)) for c in d["color"])
        for i in range(n):
            rows.append(dict(vector=vec, method=STRIP_METHODS[s], role=role, index=int(round(idx[i])),
                             value_rel=rel[i], page_x=xy[i, 0], page_y=xy[i, 1], stroke_rgb=rgb, lwd_pt=d["width"]))
    return pd.DataFrame(rows)

def main(arg, outdir):
    P = load_paths(arg); os.makedirs(outdir, exist_ok=True)
    recs = []; summary = []
    for key, label in IMAGES:
        M, frame = extract_image(P[key])
        lv = np.unique(M)
        summary.append(dict(image=label, nrow=M.shape[0], ncol=M.shape[1], n_gray_levels=len(lv),
                            gray_levels=" ".join(f"{v:.4f}" for v in lv)))
        r, c = np.indices(M.shape)
        recs.append(pd.DataFrame(dict(image=label, row=r.ravel() + 1, col=c.ravel() + 1, gray=M.ravel())))
    H = pd.concat(recs); H.to_csv(os.path.join(outdir, "fig1_heatmaps_from_pdf.csv"), index=False)
    L = pd.concat([extract_strips(P["u1"], "u"), extract_strips(P["v1"], "v")], ignore_index=True)
    L["value_std"] = np.nan
    for (vec, m), g in L.groupby(["vector", "method"]):
        t = g.loc[g.role == "truth", "value_rel"]
        L.loc[g.index, "value_std"] = (g["value_rel"] - t.min()) / (t.max() - t.min())
    L.to_csv(os.path.join(outdir, "fig1_vectors_from_pdf.csv"), index=False)
    S = pd.DataFrame(summary); print(S.to_string(index=False))
    print(L.groupby(["vector", "method", "role"]).size().unstack().to_string())
    nr, nc = S.loc[0, "nrow"], S.loc[0, "ncol"]
    nu = L.query("vector=='u'").groupby(["method", "role"]).size().unique()
    nv = L.query("vector=='v'").groupby(["method", "role"]).size().unique()
    print(f"u lengths {nu} vs nrow {nr}; v lengths {nv} vs ncol {nc}; match:",
          list(nu) == [nr] and list(nv) == [nc])

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
