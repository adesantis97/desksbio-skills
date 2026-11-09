#!/usr/bin/env python3
"""Invasion assay (Transwell, crystal violet) — invaded cells per field, with human review.

The primary measure is the **number of cells per field** (cells), counted by shape: the stained
mask is split into its lobes (peaks of the distance transform), one point per cell, so a large
spread cell counts once. Every photo is brought to the same width (P["work_w"]) before counting,
so fields shot at different resolutions are comparable.
On a confluent field (area_pct > P["confluent"]) cells merge into a mesh: the count is a minimum.
The fraction of the field covered by the stain (area_pct) is kept as a second measure.

Commands:
  detect    <CONDITION_DIR> [--thr X]          measure every 20X photo of the condition, write
                                               _invasion/fields.csv + overlays + review.html
  review    <CONDITION_DIR>                    slider page: one threshold for the whole condition,
                                               exclude the bad fields; save and recompute
  apply     <CONDITION_DIR>                    re-apply _invasion/corrections.json -> fields.csv
  run       <ROOT> [--thr X]                   detect on every condition of the tree
  aggregate <ROOT> [--ctrl CTRL] [--plot]      mean per condition, fold change, mean over replicates

Expected tree (three levels under the root, names are free):
  ROOT/<replicate>/<line>/<condition>/*.bmp

The threshold is the only subjective parameter: keep it **the same within one replicate+line**,
otherwise the percentages across conditions are not comparable.
"""
import argparse, warnings, base64, io, json, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage as ndi
from skimage import feature, measure, morphology, transform

# scikit-image 0.26 renamed a parameter of remove_small_objects; the result is the same
warnings.filterwarnings("ignore", category=FutureWarning)

HERE = Path(__file__).resolve().parent
OUTDIR = "_invasion"
EXTS = {".bmp", ".tif", ".tiff", ".png", ".jpg", ".jpeg"}

P = dict(
    thr=0.10,          # threshold on the stain index; 0.10 calibrated on 20X crystal-violet membranes
    min_obj=20,        # smaller objects = debris / membrane pores (px at full resolution)
    smooth=1.5,        # filter sigma before thresholding (px)
    confluent=60.0,    # above this area_pct the field is confluent: the cell count is a minimum
    work_w=1536,       # every photo is resized to this width before counting cells
    cell_r=22,         # radius of one cell at work_w (px); 22 calibrated on 20X A375 membranes
    hole=400,          # holes smaller than this (pores) are filled inside a cell (px at work_w)
    min_cell=500,      # smaller objects are debris (px at work_w)
    preview_w=760,     # thumbnail width on the review page
)
# scale used to save the stain index as an 8-bit PNG: v = (s - S0) / S1 * 255
S0, S1 = -0.20, 2.20
LADDER = [round(x, 3) for x in np.arange(0.00, 0.605, 0.01)]  # precomputed thresholds, exact values


# ---------------------------------------------------------------- image
def flatfield(ch):
    """Vignetting = quadratic surface fitted on the bright pixels only (membrane, not cells)."""
    h, w = ch.shape
    yy, xx = np.mgrid[0:h, 0:w]
    y = (yy / h - .5).ravel(); x = (xx / w - .5).ravel(); v = ch.ravel()
    sel = v >= np.percentile(v, 70)
    A = np.stack([np.ones_like(x), x, y, x * x, y * y, x * y], 1)
    c, *_ = np.linalg.lstsq(A[sel], v[sel], rcond=None)
    return (A @ c).reshape(h, w)


def stain_index(path):
    """Crystal violet absorbs green and transmits blue -> OD(G) - OD(B).
    The optical-density difference is insensitive to the photo's white balance."""
    rgb = np.asarray(Image.open(path).convert("RGB")).astype(np.float32)
    flat = np.empty_like(rgb)
    for k in range(3):
        flat[..., k] = np.clip(rgb[..., k] / np.maximum(flatfield(rgb[..., k]), 1e-3) * 255., 0, 255)
    od = -np.log10((flat + 1.) / 256.)
    return ndi.gaussian_filter(od[..., 1] - od[..., 2], P["smooth"]), rgb


def mask_at(s, thr):
    m = s > thr
    return morphology.remove_small_objects(m, P["min_obj"])


def cells_at(s, thr):
    """Cell centres at threshold thr, in work_w coordinates: one point per lobe of the stained mask."""
    k = P["work_w"] / s.shape[1]
    if abs(k - 1) > 1e-3:
        s = transform.resize(s, (round(s.shape[0] * k), P["work_w"]), anti_aliasing=True, preserve_range=True)
    m = morphology.binary_closing(s > thr, morphology.disk(2))
    m = morphology.remove_small_holes(m, P["hole"])
    m = morphology.binary_opening(m, morphology.disk(3))
    m = morphology.remove_small_objects(m, P["min_cell"])
    lab = measure.label(m)
    r = P["cell_r"]
    d = ndi.gaussian_filter(ndi.distance_transform_edt(m), 2)
    pk = feature.peak_local_max(d, min_distance=r, threshold_abs=r * 0.35, labels=lab, exclude_border=False)
    have = set(lab[tuple(pk.T)].tolist()) if len(pk) else set()
    pts = [[int(x), int(y)] for y, x in pk]
    pts += [[int(p.centroid[1]), int(p.centroid[0])] for p in measure.regionprops(lab) if p.label not in have]
    return sorted(pts, key=lambda q: (q[1], q[0])), list(m.shape)


# ---------------------------------------------------------------- folder I/O
def photos(cond_dir, mag="20"):
    """Photos of the condition. Skips the overview magnification (10X) and the results folder."""
    out = []
    for p in sorted(Path(cond_dir).iterdir()):
        if p.is_dir() or p.suffix.lower() not in EXTS or p.name.startswith("."):
            continue
        tags = re.findall(r"(\d+)\s*[xX]", p.name)
        if tags and mag and tags[-1] != mag:
            continue
        out.append(p)
    return out


def conditions(root):
    """ROOT/<replicate>/<line>/<condition> — the leaf folders with at least one photo."""
    out = []
    for rep in sorted(p for p in Path(root).iterdir() if p.is_dir() and not p.name.startswith(("_", "."))):
        for line in sorted(p for p in rep.iterdir() if p.is_dir() and not p.name.startswith(("_", "."))):
            for cond in sorted(p for p in line.iterdir() if p.is_dir() and not p.name.startswith(("_", "."))):
                if photos(cond):
                    out.append((rep.name, line.name, cond.name, cond))
    return out


def corrections(out):
    p = out / "corrections.json"
    return json.loads(p.read_text()) if p.exists() else {}


# ---------------------------------------------------------------- detect
def cmd_detect(a):
    cond = Path(a.cond).resolve()
    out = cond / OUTDIR
    out.mkdir(exist_ok=True)
    fs = photos(cond, a.mag)
    if not fs:
        sys.exit(f"no {a.mag}X photos in {cond}")
    thr = a.thr if a.thr is not None else P["thr"]
    store, rows = {}, []
    for f in fs:
        s, rgb = stain_index(f)
        curve = {f"{t:.2f}": round(100. * mask_at(s, t).mean(), 2) for t in LADDER}
        pts, wshape = cells_at(s, thr)
        store[f.name] = dict(curve=curve, shape=list(s.shape), stem=f.stem,
                             iod=round(float(s.mean()), 4), cells=pts, cells_thr=thr, work_shape=wshape)
        rows.append(len(pts))
        _preview(rgb, out / f"prev_{f.stem}.jpg")
        _stainmap(s, out / f"stain_{f.stem}.png")
        _overlay(rgb, mask_at(s, thr), pts, out / f"{f.stem}_overlay.jpg")
    (out / "detect.json").write_text(json.dumps(dict(threshold=thr, mag=a.mag, params=P,
                                                     ladder=LADDER, fields=store)))
    _write_fields(cond, out, thr, store)
    build_review_html(out)
    print(f"{cond.parent.parent.name}/{cond.parent.name}/{cond.name}: {len(fs)} fields, "
          f"threshold {thr:.2f}, mean {np.mean(rows):.0f} cells/field")


def _preview(rgb, path):
    im = Image.fromarray(rgb.astype(np.uint8))
    im.thumbnail((P["preview_w"], P["preview_w"]))
    im.save(path, quality=72)


def _stainmap(s, path):
    """8-bit stain index at the thumbnail scale: the page applies the threshold to it.
    The numbers on the page are still the exact ones, computed at full resolution in `curve`."""
    v = np.clip((s - S0) / (S1 - S0) * 255., 0, 255).astype(np.uint8)
    im = Image.fromarray(v)
    im.thumbnail((P["preview_w"], P["preview_w"]), Image.BILINEAR)
    im.save(path, optimize=True)


def _overlay(rgb, m, pts, path):
    from PIL import ImageDraw
    ov = rgb.copy()
    ov[m] = 0.55 * ov[m] + 0.45 * np.array([0., 230., 90.])
    im = Image.fromarray(ov.astype(np.uint8))
    im.thumbnail((1200, 1200))
    k = im.width / P["work_w"]; dr = ImageDraw.Draw(im)
    for x, y in pts:
        x, y = x * k, y * k
        dr.ellipse([x - 6, y - 6, x + 6, y + 6], fill=(255, 225, 0), outline=(0, 0, 0))
    im.save(path, quality=80)


def final_cells(cond, out, store, thr, corr, write=True):
    """Points of every field at thr, after the clicks of the review page.
    Clicks are kept only if they were made on points computed at the same threshold."""
    changed = False
    res = {}
    edits = corr.get("cells", {}) if corr.get("cells_thr") == thr else {}
    for name, f in store.items():
        if f.get("cells_thr") != thr or "cells" not in f:
            s, _ = stain_index(cond / name)
            f["cells"], f["work_shape"] = cells_at(s, thr); f["cells_thr"] = thr; changed = True
        e = edits.get(name, {})
        rem = set(e.get("remove", []))
        pts = [p for i, p in enumerate(f["cells"]) if i not in rem] + [list(map(int, q)) for q in e.get("add", [])]
        res[name] = dict(pts=pts, auto=len(f["cells"]), removed=len(rem & set(range(len(f["cells"])))),
                         added=len(e.get("add", [])))
    if changed and write:
        d = json.loads((out / "detect.json").read_text()); d["fields"] = store
        (out / "detect.json").write_text(json.dumps(d))
    return res


def _write_fields(cond, out, thr, store):
    corr = corrections(out)
    thr = corr.get("threshold", thr)
    drop = set(corr.get("exclude", []))
    recs = []
    cells = final_cells(cond, out, store, thr, corr)
    for name, f in store.items():
        area = f["curve"].get(f"{thr:.2f}")
        if area is None:                     # threshold outside the precomputed ladder: re-measure
            s, _ = stain_index(cond / name)
            area = round(100. * mask_at(s, thr).mean(), 2)
        c = cells[name]
        note = corr.get("notes", {}).get(name, "")
        if area > P["confluent"]:
            note = (note + "; " if note else "") + "confluent: cell count is a minimum"
        recs.append(dict(photo=name, cells=len(c["pts"]), cells_auto=c["auto"], removed=c["removed"],
                         added=c["added"], area_pct=area, iod=f["iod"], used=name not in drop, note=note))
    df = pd.DataFrame(recs).sort_values("photo")
    df.insert(0, "condition", cond.name)
    df.insert(0, "line", cond.parent.name)
    df.insert(0, "replicate", cond.parent.parent.name)
    df["threshold"] = thr
    df.to_csv(out / "fields.csv", index=False)
    k = df[df.used]
    pd.DataFrame([dict(replicate=cond.parent.parent.name, line=cond.parent.name, condition=cond.name,
                       n_fields=len(k), excluded=len(df) - len(k), threshold=thr,
                       cells_mean=round(k.cells.mean(), 1),
                       cells_sd=round(k.cells.std(ddof=1), 1) if len(k) > 1 else np.nan,
                       cells_cv_pct=round(100 * k.cells.std(ddof=1) / k.cells.mean(), 1) if len(k) > 1 and k.cells.mean() else np.nan,
                       confluent_fields=int((k.area_pct > P["confluent"]).sum()),
                       area_pct_mean=round(k.area_pct.mean(), 2),
                       area_pct_sd=round(k.area_pct.std(ddof=1), 2) if len(k) > 1 else np.nan,
                       cv_pct=round(100 * k.area_pct.std(ddof=1) / k.area_pct.mean(), 1) if len(k) > 1 else np.nan,
                       iod_mean=round(k.iod.mean(), 4))]).to_csv(out / "condition.csv", index=False)
    return df


def cmd_apply(a):
    cond = Path(a.cond).resolve(); out = cond / OUTDIR
    store = json.loads((out / "detect.json").read_text())
    df = _write_fields(cond, out, store["threshold"], store["fields"])
    corr = corrections(out); thr = corr.get("threshold", store["threshold"])
    cells = final_cells(cond, out, store["fields"], thr, corr, write=False)
    for f in photos(cond, store.get("mag", "20")):
        if f.name in df.photo.values:
            s, rgb = stain_index(f)
            _overlay(rgb, mask_at(s, thr), cells[f.name]["pts"], out / f"{f.stem}_overlay.jpg")
    build_review_html(out)
    k = df[df.used]
    print(f"{cond.name}: {len(k)}/{len(df)} fields, threshold {thr:.2f}, "
          f"{k.cells.mean():.0f} ± {k.cells.std(ddof=1):.0f} cells/field, area {k.area_pct.mean():.1f} %")
    return df


# ---------------------------------------------------------------- run / aggregate
def cmd_run(a):
    for rep, line, cname, cond in conditions(a.root):
        cmd_detect(argparse.Namespace(cond=str(cond), thr=a.thr, mag=a.mag))


def cmd_aggregate(a):
    root = Path(a.root).resolve()
    fields, conds = [], []
    for rep, line, cname, cond in conditions(root):
        fp = cond / OUTDIR / "fields.csv"
        if not fp.exists():
            print(f"⚠️  missing {fp} — run detect first", file=sys.stderr); continue
        fields.append(pd.read_csv(fp))
        conds.append(pd.read_csv(cond / OUTDIR / "condition.csv"))
    if not conds:
        sys.exit("no measured condition")
    F = pd.concat(fields, ignore_index=True)
    C = pd.concat(conds, ignore_index=True)

    ctrl = (C[C.condition.str.upper() == a.ctrl.upper()]
            .set_index(["replicate", "line"])[["cells_mean", "area_pct_mean"]]
            .rename(columns={"cells_mean": "ctrl_cells", "area_pct_mean": "ctrl_area_pct"}))
    C = C.join(ctrl, on=["replicate", "line"])
    C["fold_change"] = (C.cells_mean / C.ctrl_cells).round(3)
    C["area_fold_change"] = (C.area_pct_mean / C.ctrl_area_pct).round(3)
    C = C.sort_values(["line", "condition", "replicate"])
    C.to_csv(root / "invasion_conditions.csv", index=False)
    F.to_csv(root / "invasion_fields.csv", index=False)

    g = C.groupby(["line", "condition"], sort=False)
    S = g.agg(n_rep=("fold_change", "count"),
              cells_mean=("cells_mean", "mean"),
              confluent_fields=("confluent_fields", "sum"),
              fc_mean=("fold_change", "mean"),
              fc_sd=("fold_change", lambda v: v.std(ddof=1)),
              fc_sem=("fold_change", lambda v: v.std(ddof=1) / np.sqrt(v.count()))).reset_index()
    S["p_vs_ctrl"] = [_ttest1(C[(C.line == l) & (C.condition == c)].fold_change.dropna().values, 1.0)
                      for l, c in zip(S.line, S.condition)]
    S["area_pct_mean"] = g["area_pct_mean"].mean().values
    S["area_fc_mean"] = g["area_fold_change"].mean().values
    for col in ["cells_mean", "area_pct_mean", "fc_mean", "fc_sd", "fc_sem", "area_fc_mean"]:
        S[col] = S[col].round(3)
    S.to_csv(root / "invasion_summary.csv", index=False)
    print(S.to_string(index=False))
    xl = ""
    try:  # one Excel file with the three tables, if openpyxl is there
        with pd.ExcelWriter(root / "invasion.xlsx") as w:
            for name, df in (("summary", S), ("conditions", C), ("fields", F)):
                df.to_excel(w, sheet_name=name, index=False)
                ws = w.sheets[name]; ws.freeze_panes = "A2"
                for col in ws.columns:
                    ws.column_dimensions[col[0].column_letter].width = max(10, min(40, max(len(str(c.value or "")) for c in col) + 2))
        xl = " · invasion.xlsx"
    except ImportError:
        print("(pip install openpyxl for invasion.xlsx)", file=sys.stderr)
    print(f"\nwrote  invasion_fields.csv · invasion_conditions.csv · invasion_summary.csv{xl}  to {root}")
    if a.plot:
        import subprocess
        subprocess.run(["Rscript", str(HERE / "plot_invasion.R"), str(root / "invasion_conditions.csv"),
                        str(root), a.colour], check=True)


def _ttest1(v, mu):
    """One-sample t-test against mu. Without scipy.stats: p from Student's t via betainc."""
    n = len(v)
    if n < 3 or np.std(v, ddof=1) == 0:
        return np.nan
    from scipy import special
    t = (np.mean(v) - mu) / (np.std(v, ddof=1) / np.sqrt(n))
    df = n - 1
    return round(float(special.betainc(df / 2, 0.5, df / (df + t * t))), 4)


# ---------------------------------------------------------------- review
def _b64(path, mime):
    return f"data:{mime};base64," + base64.b64encode(Path(path).read_bytes()).decode()


def build_review_html(out):
    store = json.loads((out / "detect.json").read_text())
    cond = out.parent
    imgs = {name: dict(prev=_b64(out / f"prev_{f['stem']}.jpg", "image/jpeg"),
                       stain=_b64(out / f"stain_{f['stem']}.png", "image/png"))
            for name, f in store["fields"].items()}
    tpl = (HERE / "review_template.html").read_text()
    data = dict(store=store, images=imgs, corrections=corrections(out),
                work_w=P["work_w"], cell_r=P["cell_r"], confluent=P["confluent"],
                title=f"{cond.parent.parent.name} · {cond.parent.name} · {cond.name}")
    (out / "review.html").write_text(tpl.replace("/*__DATA__*/null", json.dumps(data)))


def cmd_review(a):
    import http.server, threading, webbrowser
    cond = Path(a.cond).resolve(); out = cond / OUTDIR
    build_review_html(out)

    class H(http.server.BaseHTTPRequestHandler):
        def log_message(self, *args): pass

        def _send(self, code, body, ctype="text/plain; charset=utf-8"):
            b = body.encode() if isinstance(body, str) else body
            self.send_response(code); self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(b))); self.end_headers(); self.wfile.write(b)

        def do_GET(self):
            self._send(200, (out / "review.html").read_bytes(), "text/html; charset=utf-8")

        def do_POST(self):
            body = self.rfile.read(int(self.headers.get("Content-Length") or 0)).decode()
            if self.path == "/save":
                (out / "corrections.json").write_text(body)
                buf, old = io.StringIO(), sys.stdout
                sys.stdout = buf
                try:
                    cmd_apply(argparse.Namespace(cond=str(cond)))
                finally:
                    sys.stdout = old
                print(buf.getvalue().strip(), flush=True)
                self._send(200, buf.getvalue())
            elif self.path == "/done":
                self._send(200, "ok"); threading.Thread(target=srv.shutdown).start()

    srv = http.server.HTTPServer(("127.0.0.1", a.port), H)
    url = f"http://127.0.0.1:{srv.server_port}/"
    print(f"review open at {url} — 'Save and recompute', then 'Done'", flush=True)
    webbrowser.open(url)
    srv.serve_forever()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    d = sp.add_parser("detect"); d.add_argument("cond"); d.add_argument("--thr", type=float)
    d.add_argument("--mag", default="20"); d.set_defaults(fn=cmd_detect)
    p = sp.add_parser("apply"); p.add_argument("cond"); p.set_defaults(fn=cmd_apply)
    r = sp.add_parser("review"); r.add_argument("cond"); r.add_argument("--port", type=int, default=0)
    r.set_defaults(fn=cmd_review)
    u = sp.add_parser("run"); u.add_argument("root"); u.add_argument("--thr", type=float)
    u.add_argument("--mag", default="20"); u.set_defaults(fn=cmd_run)
    g = sp.add_parser("aggregate"); g.add_argument("root"); g.add_argument("--ctrl", default="CTRL")
    g.add_argument("--plot", action="store_true"); g.add_argument("--colour", default="darkred")
    g.set_defaults(fn=cmd_aggregate)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
