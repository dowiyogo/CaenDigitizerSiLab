#!/usr/bin/python3 -I
# set_png_dpi.py — SOCHIFI 2026-11-13: escribe 300 ppp en los metadatos (pHYs) de los PNG del póster, sin cambiar los píxeles (2400 × 1800 =
# 20.3 × 15.2 cm a 300 ppp). Uso: set_png_dpi.py <png>...
import sys

from PIL import Image

for p in sys.argv[1:]:
    im = Image.open(p)
    im.load()
    im.save(p, dpi=(300, 300))
    print(p, Image.open(p).info.get("dpi"), Image.open(p).size)
