"""論文 Table 2 の「AIA 94 Å = 7.20 DN/s」を使って箱を独立に検証する。

Table 2 の最終行は EIS 輝線ではなく **AIA 94 Å（Fe XVIII 分離後）の強度**で、
同じ inter-moss 箱で測った値。EIS の輝線比とは完全に独立な情報なので、
「箱の位置が正しいか」を別ルートで確かめられる。

データは JSOC の synoptic アーカイブ（1024x1024, 2.4″/px, **登録不要**）。
    http://jsoc.stanford.edu/data/aia/synoptic/YYYY/MM/DD/HHHH/AIAyyyymmdd_hhmm_wwww.fits
フルディスクの level-1 (4096², 65 MB) は sdo7.nascom.nasa.gov が遅くて
タイムアウトするので、こちらを使う。2.4″/px でも 15×23″ の箱なら
6×10 画素あり、平均値の確認には十分。

    python scripts/check_fe18_box.py
"""
import glob
import sys

import numpy as np
import astropy.units as u
import sunpy.map

sys.path.insert(0, "scripts")
from aia_fe18 import fe18 as fe18_from_aia  # noqa: E402

# 論文 Figure 2 の緑枠から実測した region 7 の箱（scripts/extract_paper_boxes.py）
PAPER_BOX = dict(x0=-321.8, x1=-306.4, y0=202.6, y1=226.2)
PAPER_VALUE = 7.20     # Table 2 の AIA 94 Å I_obs [DN/s]
PAPER_SIGMA = 1.40


def load(w):
    f = glob.glob(f"data/sdo/synoptic/AIA*_{w:04d}.fits")
    if not f:
        raise SystemExit(f"AIA {w} が data/sdo/synoptic/ に無い")
    m = sunpy.map.Map(f[0])
    exp = m.meta.get("exptime", 1.0)
    print(f"  AIA {w:4d}: {m.data.shape}  exptime={exp:.3f}s  "
          f"lvl={m.meta.get('lvl_num')}  scale={m.scale[0]:.3f}")
    return sunpy.map.Map(m.data / exp, m.meta)      # DN/s にする


def box_mean(m, b):
    from astropy.coordinates import SkyCoord
    bl = SkyCoord(b["x0"] * u.arcsec, b["y0"] * u.arcsec, frame=m.coordinate_frame)
    tr = SkyCoord(b["x1"] * u.arcsec, b["y1"] * u.arcsec, frame=m.coordinate_frame)
    sub = m.submap(bl, top_right=tr)
    return float(np.nanmean(sub.data)), sub.data.shape


def main():
    print("AIA synoptic (1024x1024) を読む:")
    a94, a171, a193 = load(94), load(171), load(193)

    fe18 = fe18_from_aia(a94.data, a171.data, a193.data)
    m18 = sunpy.map.Map(fe18, a94.meta)

    print(f"\n論文の箱 X=[{PAPER_BOX['x0']}, {PAPER_BOX['x1']}] "
          f"Y=[{PAPER_BOX['y0']}, {PAPER_BOX['y1']}]")
    v, shp = box_mean(m18, PAPER_BOX)
    print(f"  Fe XVIII 平均 = {v:8.3f} DN/s   ({shp[0]}x{shp[1]} px)")
    print(f"  論文 Table 2  = {PAPER_VALUE:8.3f} +- {PAPER_SIGMA:.2f} DN/s")
    print(f"  ratio         = {v/PAPER_VALUE:8.3f}")

    # 周辺を走査して 7.20 DN/s になる場所を探す
    print("\n近傍で Fe XVIII の箱平均が 7.20 DN/s になる場所:")
    print(f"{'Xc':>8} {'Yc':>8} {'FeXVIII':>9} {'ratio':>7}")
    dx = PAPER_BOX["x1"] - PAPER_BOX["x0"]
    dy = PAPER_BOX["y1"] - PAPER_BOX["y0"]
    hits = []
    for xc in np.arange(-360, -270, 5.0):
        for yc in np.arange(170, 270, 5.0):
            b = dict(x0=xc - dx / 2, x1=xc + dx / 2, y0=yc - dy / 2, y1=yc + dy / 2)
            try:
                vv, _ = box_mean(m18, b)
            except Exception:
                continue
            if np.isfinite(vv):
                hits.append((abs(np.log10(max(vv, 1e-3) / PAPER_VALUE)), xc, yc, vv))
    hits.sort()
    for d, xc, yc, vv in hits[:10]:
        print(f"{xc:8.1f} {yc:8.1f} {vv:9.3f} {vv/PAPER_VALUE:7.3f}")


if __name__ == "__main__":
    main()
