#!/usr/bin/env python3
"""Make the generated figures journal-ready and viewable.

Two jobs, both idempotent, run after make_figs.py / flow_figure.py:

1. Journal files in figures_submission/ (PNG and TIFF) are rewritten as plain
   RGB with an explicit DPI of at least 300.  Matplotlib writes RGBA by default;
   the alpha channel is fully opaque for these figures, so dropping it is
   lossless and removes a common cause of figures failing to open in viewers,
   and the 300 dpi floor matches the journal's minimum.  Pixel data is never
   resampled - only the declared resolution changes.

2. figures_preview/ is refreshed with a full-resolution JPEG (quality 88) and a
   lossless RGB PNG for every figure, plus a multi-page AllFigures.pdf, a
   ContactSheet.jpg overview and a README.  These formats open in any viewer,
   which EPS and TIFF do not.

Usage:  NHANES_ROOT=/path/to/workspace python3 normalise_figures.py
"""
import os
import glob
from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None

ROOT = os.environ.get("NHANES_ROOT") or os.getcwd()
SRC = os.path.join(ROOT, "figures_submission")
PREV = os.path.join(ROOT, "figures_preview")
os.makedirs(PREV, exist_ok=True)

MIN_DPI = 300
PDF_LONG_EDGE = 2400
SHEET_W = 1750

FIGS = [
    ("fig1_forest", "Fig1_forest", "Figure 1. Forest plot, single-metal associations"),
    ("fig2a_qgcomp_weights_3m", "Fig2a_qgcomp_weights", "Figure 2a. Quantile g-computation weights (3 metals)"),
    ("fig3_quartiles", "Fig3_quartiles", "Figure 3. Quartile dose-response"),
    ("fig4_rcs", "Fig4_rcs", "Figure 4. Restricted cubic splines (3 metals)"),
    ("fig5_rcs_se_mn", "Fig5_rcs_se_mn", "Figure 5. Restricted cubic splines (Se / Mn)"),
    ("figS1_flow", "FigS1_flow", "Figure S1. Participant flow chart"),
    ("graphical_abstract", "GraphicalAbstract", "Graphical abstract"),
]


def journal_ready():
    print("== 1. journal files: plain RGB with an explicit >= 300 dpi ==")
    for stem, _, _ in FIGS:
        png, tif = os.path.join(SRC, stem + ".png"), os.path.join(SRC, stem + ".tiff")
        if not os.path.exists(png):
            print(f"   skip {stem}: not found")
            continue
        im = Image.open(png)
        cur = im.info.get("dpi", (72, 72))
        dpi = (max(MIN_DPI, round(cur[0])), max(MIN_DPI, round(cur[1])))
        rgb = im.convert("RGB")
        rgb.save(png, "PNG", dpi=dpi, optimize=True)
        if os.path.exists(tif):
            rgb.save(tif, "TIFF", compression="tiff_lzw", dpi=dpi)
        print(f"   {stem:26s} {im.mode} -> RGB  dpi {tuple(round(x) for x in cur)} -> {dpi}  px {im.size}")


def previews():
    print("== 2. figures_preview/: viewable copies, PDF and contact sheet ==")
    pages = []
    for stem, out, _ in FIGS:
        png = os.path.join(SRC, stem + ".png")
        if not os.path.exists(png):
            continue
        im = Image.open(png).convert("RGB")
        dpi = im.info.get("dpi", (300, 300))
        im.save(os.path.join(PREV, out + ".jpg"), "JPEG", quality=88, dpi=dpi, optimize=True)
        im.save(os.path.join(PREV, out + ".png"), "PNG", dpi=dpi, optimize=True)
        w, h = im.size
        s = PDF_LONG_EDGE / max(w, h)
        pages.append(im.resize((max(1, round(w * s)), max(1, round(h * s))), Image.LANCZOS) if s < 1 else im.copy())
        print(f"   {out}.jpg / .png  ({w}x{h})")
    if pages:
        pages[0].save(os.path.join(PREV, "AllFigures.pdf"), "PDF", resolution=300.0,
                      save_all=True, append_images=pages[1:], quality=88)
        print(f"   AllFigures.pdf  ({len(pages)} pages)")
    try:
        font = ImageFont.truetype("C:/Windows/Fonts/msyh.ttc", 30)
    except Exception:
        try:
            font = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 30)
        except Exception:
            font = ImageFont.load_default()
    lab_h, pad, gap = 44, 28, 34
    rows = []
    for im, (_, _, label) in zip(pages, FIGS):
        w, h = im.size
        rows.append((im.resize((SHEET_W, max(1, round(h * SHEET_W / w))), Image.LANCZOS), label))
    sheet = Image.new("RGB", (SHEET_W + 2 * pad, pad + sum(r.height + lab_h + gap for r, _ in rows) + pad), "white")
    d = ImageDraw.Draw(sheet)
    y = pad
    for r, label in rows:
        d.text((pad, y), label, fill="black", font=font)
        y += lab_h
        sheet.paste(r, (pad, y))
        y += r.height
        d.line([(pad, y + gap // 2), (pad + SHEET_W, y + gap // 2)], fill=(200, 200, 200), width=2)
        y += gap
    sheet.save(os.path.join(PREV, "ContactSheet.jpg"), "JPEG", quality=86, optimize=True)
    print(f"   ContactSheet.jpg  ({sheet.width}x{sheet.height})")


if __name__ == "__main__":
    journal_ready()
    previews()
    print("done")
