#!/usr/bin/env python3
"""Colony formation assay — colony counting with human review.

Commands:
  detect  <photo> --design design.csv --out DIR   find dishes + colonies, write
                                                  detections.json, numbered overlays, review.html
  review  DIR                                     open the review page in the browser; "Save" writes
                                                  DIR/corrections.json and runs apply by itself
  apply   DIR [--corrections corrections.json]    apply the corrections (removed/added)
                                                  and write counts.csv: colonies + mean size + fold change;
                                                  colonies.csv: one row per colony (size in px)
                                                  --plot [--colour X]: ggplot figures (optional)
  (several photos in the same DIR: run detect once per photo, apply once)

design.csv (one row per dish, order = rows from the top, then left to right):
  image,dish,line,condition,dose
  plate1.jpg,1,HeLa,CTRL,10nM

Corrections (from review.html or written by hand in chat):
  {"<image>": {"<dish>": {"remove": [3, 17], "add": [[x, y, r], ...]}}}
"""
import argparse, warnings, base64, io, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage import feature, measure, segmentation, transform, morphology, filters

# scikit-image 0.26 renamed a parameter of remove_small_objects; the result is the same
warnings.filterwarnings("ignore", category=FutureWarning)

HERE = Path(__file__).resolve().parent

# ---------- default parameters (overridable from the CLI) ----------
P = dict(
    rim_frac=0.90,        # keep the inner 90% of the radius: the dish rim is not a colony
    stain_thr=None,       # threshold on (B - G); None = automatic (Otsu inside the dishes, min 25)
    min_area_frac=4e-5,   # minimum colony area, as a fraction of the dish area
    min_peak_dist=4,      # minimum distance between two centres when splitting merged colonies (px)
)

FONTS = ["/System/Library/Fonts/Supplemental/Arial.ttf", "DejaVuSans.ttf", "arial.ttf"]
FONTS_BOLD = ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "DejaVuSans-Bold.ttf", "arialbd.ttf"]


def load_font(names, size):
    for n in names:
        try:
            return ImageFont.truetype(n, size)
        except OSError:
            pass
    return ImageFont.load_default()


def load_rgb(path):
    return np.asarray(Image.open(path).convert("RGB")).astype(np.int16)


def stain_map(rgb):
    """Crystal violet: high blue, low green. The grey background has B ~ G."""
    return np.clip(rgb[..., 2] - rgb[..., 1], 0, 255).astype(np.float32)


def find_dishes(rgb, n_expected):
    """Dishes = large circles. Hough transform on the edges of the greyscale image."""
    gray = rgb.mean(axis=2) / 255.0
    h, w = gray.shape
    edges = feature.canny(gray, sigma=2.5)
    # dishes from ~1/8 of the short side (many dishes per photo) up to the full short side (one dish)
    rmin, rmax = int(min(h, w) * 0.12), int(min(h, w) * 0.50)
    radii = np.arange(rmin, rmax, max(2, (rmax - rmin) // 120))
    acc = transform.hough_circle(edges, radii)
    _, cx, cy, rr = transform.hough_circle_peaks(
        acc, radii, total_num_peaks=n_expected, min_xdistance=rmin, min_ydistance=rmin)
    dishes = [dict(cx=int(x), cy=int(y), r=int(r)) for x, y, r in zip(cx, cy, rr)]
    # reading order: by rows (tolerance of one minimum radius), then left to right
    dishes.sort(key=lambda d: (round(d["cy"] / (rmin)), d["cx"]))
    return dishes


def detect_in_dish(stain, dish, thr):
    h, w = stain.shape
    yy, xx = np.ogrid[:h, :w]
    inner = (xx - dish["cx"]) ** 2 + (yy - dish["cy"]) ** 2 <= (dish["r"] * P["rim_frac"]) ** 2
    s = ndi.gaussian_filter(stain, 1.0)
    mask = (s > thr) & inner
    mask = morphology.remove_small_objects(mask, 3)
    min_area = max(4, P["min_area_frac"] * np.pi * dish["r"] ** 2)
    # split touching colonies: watershed on the distance transform
    dist = ndi.distance_transform_edt(mask)
    peaks = feature.peak_local_max(dist, min_distance=P["min_peak_dist"], labels=measure.label(mask),
                                   exclude_border=False)
    markers = np.zeros_like(mask, dtype=np.int32)
    for i, (py, px) in enumerate(peaks, 1):
        markers[py, px] = i
    lab = segmentation.watershed(-dist, markers, mask=mask)
    cols = []
    for rp in measure.regionprops(lab, intensity_image=s):
        if rp.area < min_area:
            continue
        cy, cx = rp.centroid
        cols.append(dict(x=round(float(cx), 1), y=round(float(cy), 1),
                         r=round(float(np.sqrt(rp.area / np.pi)), 1),
                         area=int(rp.area), stain=round(float(rp.intensity_mean), 1)))
    cols.sort(key=lambda c: (c["y"], c["x"]))
    for i, c in enumerate(cols, 1):
        c["id"] = i
    return cols


def overlay_png(rgb, dishes, dets, out_png, corrections=None):
    im = Image.fromarray(rgb.astype(np.uint8)).convert("RGB")
    scale = 2  # enlarged overlay: the numbers must be readable
    im = im.resize((im.width * scale, im.height * scale), Image.LANCZOS)
    dr = ImageDraw.Draw(im)
    font, big = load_font(FONTS, 11), load_font(FONTS_BOLD, 28)
    for d in dishes:
        k = str(d["dish"])
        cr = (corrections or {}).get(k) or (corrections or {}).get("P" + k, {})  # "1" or "P1"
        removed = set(cr.get("remove", []))
        x, y, r = d["cx"] * scale, d["cy"] * scale, d["r"] * scale
        dr.ellipse([x - r, y - r, x + r, y + r], outline=(0, 160, 255), width=2)
        dr.text((x - r, y - r - 30), f"P{k}", fill=(0, 120, 220), font=big)
        for c in dets[k]:
            cx, cy, cr_ = c["x"] * scale, c["y"] * scale, max(c["r"] * scale, 5) + 2
            col = (220, 30, 30) if c["id"] in removed else (0, 200, 60)
            dr.ellipse([cx - cr_, cy - cr_, cx + cr_, cy + cr_], outline=col, width=2)
            dr.text((cx + cr_, cy - cr_), str(c["id"]), fill=col, font=font)
        for a in cr.get("add", []):
            ax, ay, ar = a[0] * scale, a[1] * scale, max(a[2] * scale if len(a) > 2 else 8, 6)
            dr.ellipse([ax - ar, ay - ar, ax + ar, ay + ar], outline=(255, 150, 0), width=2)
    im.save(out_png)


def cmd_detect(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    design = pd.read_csv(a.design)
    img = Path(a.image)
    rows = design[design.image == img.name]
    if rows.empty:
        sys.exit(f"design.csv has no rows for {img.name}")
    if a.thr is not None:
        P["stain_thr"] = a.thr
    rgb = load_rgb(img)
    stain = stain_map(rgb)
    dishes = find_dishes(rgb, len(rows))
    h, w = stain.shape
    yy, xx = np.ogrid[:h, :w]
    inside = np.zeros((h, w), bool)
    for d in dishes:
        inside |= (xx - d["cx"]) ** 2 + (yy - d["cy"]) ** 2 <= (d["r"] * P["rim_frac"]) ** 2
    thr = P["stain_thr"] if P["stain_thr"] is not None else max(25.0, float(filters.threshold_otsu(stain[inside])))
    dets = {}
    for i, d in enumerate(dishes, 1):
        d["dish"] = i
        dets[str(i)] = detect_in_dish(stain, d, thr)
    # copy the image next to the results: DIR stays self-contained
    if img.resolve() != (out / img.name).resolve():
        (out / img.name).write_bytes(img.read_bytes())
    if not (out / "design.csv").exists():
        design.to_csv(out / "design.csv", index=False)
    store_p = out / "detections.json"
    store = json.loads(store_p.read_text()) if store_p.exists() else {}
    store[img.name] = dict(threshold=thr, params=P, dishes=dishes, colonies=dets)
    store_p.write_text(json.dumps(store, indent=1))
    overlay_png(rgb, dishes, dets, out / f"{img.stem}_overlay.png")
    build_review_html(out)
    for d in dishes:
        print(f"{img.name} P{d['dish']}: {len(dets[str(d['dish'])])} colonies (threshold {thr:.0f})")


def build_review_html(out):
    store = json.loads((out / "detections.json").read_text())
    imgs = {name: "data:image/jpeg;base64," + base64.b64encode((out / name).read_bytes()).decode()
            for name in store}
    tpl = (HERE / "review_template.html").read_text()
    cp = out / "corrections.json"
    saved = json.loads(cp.read_text()) if cp.exists() else {}
    html = tpl.replace("/*__DATA__*/null", json.dumps(dict(store=store, images=imgs, corrections=saved)))
    (out / "review.html").write_text(html)


def load_corrections(out, path):
    p = Path(path) if path else out / "corrections.json"
    return json.loads(p.read_text()) if p.exists() else {}


def added_sizes(rgb, thr, adds):
    """Area (px) of hand-added colonies: the violet spot under the click,
    even if it was below the counting minimum area. No spot within 6 px -> the click disc."""
    mask = ndi.gaussian_filter(stain_map(rgb), 1.0) > thr
    lab = measure.label(mask)
    sizes = np.bincount(lab.ravel())
    h, w = mask.shape
    out = []
    for a in adds:
        x, y = int(round(a[0])), int(round(a[1]))
        y0, y1, x0, x1 = max(y - 6, 0), min(y + 7, h), max(x - 6, 0), min(x + 7, w)
        win = lab[y0:y1, x0:x1]
        ids = win[win > 0]
        if ids.size:
            # the spot closest to the click
            yy, xx = np.nonzero(win)
            i = np.argmin((yy + y0 - y) ** 2 + (xx + x0 - x) ** 2)
            out.append(int(sizes[win[yy[i], xx[i]]]))
        else:
            out.append(int(np.pi * (a[2] if len(a) > 2 else 6) ** 2))
    return out


def cmd_apply(a):
    out = Path(a.out)
    store = json.loads((out / "detections.json").read_text())
    corr = load_corrections(out, a.corrections)
    design = pd.read_csv(a.design) if a.design else pd.read_csv(out / "design.csv")
    recs, per_col = [], []
    for name, s in store.items():
        cimg = corr.get(name, {})
        rgb = load_rgb(out / name)
        for d in s["dishes"]:
            k = str(d["dish"])
            cr = cimg.get(k) or cimg.get("P" + k, {})  # "1" or "P1"
            rem = set(cr.get("remove", []))
            adds = cr.get("add", [])
            sizes = [c["area"] for c in s["colonies"][k] if c["id"] not in rem]
            sizes += added_sizes(rgb, s["threshold"], adds)
            per_col += [dict(image=name, dish=d["dish"], size_px=v) for v in sizes]
            recs.append(dict(image=name, dish=d["dish"], auto=len(s["colonies"][k]), removed=len(rem),
                             added=len(adds), colonies=len(sizes),
                             mean_size_px=round(float(np.mean(sizes)), 1) if sizes else np.nan))
        overlay_png(rgb, s["dishes"], s["colonies"], out / f"{Path(name).stem}_overlay.png", cimg)
    df = design.merge(pd.DataFrame(recs), on=["image", "dish"], how="left")
    # fold change: vs the CTRL of the same line (and same dose, if the column is filled in)
    keys = ["line"] + (["dose"] if "dose" in df and df["dose"].notna().any() else [])
    ctrl = (df[df.condition.str.upper() == "CTRL"].groupby(keys, dropna=False)[["colonies", "mean_size_px"]]
            .mean().rename(columns={"colonies": "ctrl_colonies", "mean_size_px": "ctrl_mean_size_px"})
            .reset_index())
    df = df.merge(ctrl, on=keys, how="left")
    df["fold_change"] = (df["colonies"] / df["ctrl_colonies"]).round(3)
    df["size_fold_change"] = (df["mean_size_px"] / df["ctrl_mean_size_px"]).round(3)
    df.to_csv(out / "counts.csv", index=False)
    design.merge(pd.DataFrame(per_col), on=["image", "dish"]).to_csv(out / "colonies.csv", index=False)
    show = [c for c in ["line", "condition", "dose", "auto", "removed", "added", "colonies", "fold_change",
                        "mean_size_px", "size_fold_change"] if c in df]
    print(df[show].to_string(index=False))
    if a.plot:
        import subprocess
        subprocess.run(["Rscript", str(HERE / "plot_colonies.R"), str(out / "counts.csv"),
                        str(out), a.colour], check=True)


def cmd_review(a):
    """Serve review.html locally: "Save" writes corrections.json into DIR and recomputes."""
    import http.server, threading, webbrowser
    out = Path(a.out).resolve()
    build_review_html(out)

    class H(http.server.BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

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
                buf = io.StringIO(); old = sys.stdout; sys.stdout = buf
                try:
                    cmd_apply(argparse.Namespace(out=str(out), corrections=None, design=None, plot=False, colour=None))
                finally:
                    sys.stdout = old
                print(f"saved {out / 'corrections.json'} and recomputed", flush=True)
                self._send(200, buf.getvalue())
            elif self.path == "/done":
                self._send(200, "ok")
                threading.Thread(target=srv.shutdown).start()

    srv = http.server.HTTPServer(("127.0.0.1", a.port), H)
    url = f"http://127.0.0.1:{srv.server_port}/"
    print(f"review open at {url} — press 'Save' on the page, then 'Done'", flush=True)
    webbrowser.open(url)
    srv.serve_forever()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    d = sp.add_parser("detect"); d.add_argument("image"); d.add_argument("--design", required=True)
    d.add_argument("--out", required=True); d.add_argument("--thr", type=float)
    d.set_defaults(fn=cmd_detect)
    p = sp.add_parser("apply"); p.add_argument("out"); p.add_argument("--corrections")
    p.add_argument("--design")
    p.add_argument("--plot", action="store_true", help="ggplot theme_bw figures (off by default)")
    p.add_argument("--colour", default="darkred", help="colour of the treated condition in the figures")
    p.set_defaults(fn=cmd_apply)
    r = sp.add_parser("review"); r.add_argument("out"); r.add_argument("--port", type=int, default=0)
    r.set_defaults(fn=cmd_review)
    a = ap.parse_args(); a.fn(a)


if __name__ == "__main__":
    main()
