"""Warren+2012 の Figure 1-3 から inter-moss 箱（緑枠）の座標を実測する。

論文は箱の座標を書いていない。しかし
  - 図のパネルは "a 400″ × 400″ region centered on the NOAA active region
    coordinates"（本文 p.3）
  - その中心座標 Xcen, Ycen は Table 1 に載っている
ので、図中の緑枠のピクセル位置から solar 座標を逆算できる。

    python scripts/extract_paper_boxes.py papers/Warren_2012_*.pdf

注意:
  - 枠線の太さ（数 arcsec 相当）ぶん、外接矩形は真の箱より大きめに出る。
  - パネルが本当に 400″ 四方で、中心が Table 1 の値ちょうどだという前提。
    region 7 について EIS データで検証した結果は docs/00_log.md を参照。
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

# Table 1: region -> (NOAA, date, Xcen, Ycen)
TABLE1 = {
    1: ("1082", "2010-06-19 01:57:44", -306.4, 439.3),
    2: ("1158", "2011-02-12 15:32:13", -248.4, -211.8),
    3: ("1082", "2010-06-21 01:46:37", 162.9, 405.2),
    4: ("1259", "2011-07-25 09:36:09", 224.7, 323.4),
    5: ("1150", "2011-01-31 11:25:19", -470.9, -250.6),
    6: ("1147", "2011-01-21 14:10:50", 26.6, 476.5),
    7: ("1243", "2011-07-02 03:38:08", -299.0, 216.6),
    8: ("1089", "2010-07-23 15:03:07", -363.4, -453.6),
    9: ("1109", "2010-09-29 23:51:36", 361.5, 261.5),
    10: ("1193", "2011-04-19 13:32:20", 36.3, 363.5),
    11: ("1190", "2011-04-11 12:00:42", -492.6, 281.0),
    12: ("1271", "2011-08-21 12:25:42", -50.8, 150.8),
    13: ("1190", "2011-04-15 01:17:19", 218.1, 304.4),
    14: ("1339", "2011-11-08 19:14:27", 88.1, 258.4),
    15: ("1339", "2011-11-10 11:33:19", 406.0, 266.8),
}
# Figure 1 = page 3 (regions 1-5), Figure 2 = page 4 (6-10), Figure 3 = page 5 (11-15)
PAGES = {3: [1, 2, 3, 4, 5], 4: [6, 7, 8, 9, 10], 5: [11, 12, 13, 14, 15]}
FOV = 400.0


def panels(mask):
    """白でない大きな正方形領域＝画像パネルを見つける。"""
    lab, _ = ndimage.label(ndimage.binary_closing(mask, np.ones((9, 9))))
    out = []
    for sl in ndimage.find_objects(lab):
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if h > 300 and w > 300 and 0.9 < h / w < 1.1:
            out.append((sl[0].start, sl[0].stop, sl[1].start, sl[1].stop))
    out.sort()
    return out


def main():
    pdf = sys.argv[1]
    tmp = Path(tempfile.mkdtemp())
    print(f"{'reg':>3} {'panel':>8}  {'solar X':>18} {'solar Y':>18}  {'size (arcsec)':>16}")
    for page, regs in PAGES.items():
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "400",
                        "-png", pdf, str(tmp / f"p{page}")], check=True)
        png = sorted(tmp.glob(f"p{page}-*.png"))[0]
        a = np.asarray(Image.open(png).convert("RGB")).astype(int)
        r, g, b = a[:, :, 0], a[:, :, 1], a[:, :, 2]
        nonwhite = (r < 245) | (g < 245) | (b < 245)
        green = (g > 90) & (g > r + 40) & (g > b + 40)
        ps = panels(nonwhite)
        if len(ps) != 25:
            print(f"  ! page {page}: パネルが {len(ps)} 個 (25 を期待)")
        rows = sorted({p[0] for p in ps})
        cols = sorted({p[2] for p in ps})
        for ir, y0 in enumerate(rows):
            if ir >= len(regs):
                continue
            reg = regs[ir]
            xc, yc = TABLE1[reg][2], TABLE1[reg][3]
            for ic, x0 in enumerate(cols):
                pan = [p for p in ps if p[0] == y0 and p[2] == x0]
                if not pan:
                    continue
                yy0, yy1, xx0, xx1 = pan[0]
                sub = green[yy0:yy1, xx0:xx1]
                if sub.sum() < 30:
                    continue
                # 複数の箱がありうる（region 13 は 2 つ）ので連結成分に分ける
                lab, n = ndimage.label(ndimage.binary_dilation(sub, np.ones((5, 5))))
                H, W = yy1 - yy0, xx1 - xx0
                for k in range(1, n + 1):
                    ys, xs = np.nonzero(lab == k)
                    if len(ys) < 30:
                        continue
                    sx0 = xc - FOV / 2 + xs.min() / W * FOV
                    sx1 = xc - FOV / 2 + xs.max() / W * FOV
                    sy1 = yc + FOV / 2 - ys.min() / H * FOV
                    sy0 = yc + FOV / 2 - ys.max() / H * FOV
                    print(f"{reg:3d} {ic:8d}  [{sx0:8.1f},{sx1:8.1f}] "
                          f"[{sy0:8.1f},{sy1:8.1f}]  {sx1-sx0:7.1f} x {sy1-sy0:6.1f}")
                break   # 1 パネルぶんで十分


if __name__ == "__main__":
    main()
