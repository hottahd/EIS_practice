"""箱を総当たりで動かして、論文 Table 2 に一番よく合う場所を探す。

背景: IDL(SSW) でも Python(eispac) でも、箱 y=[244:274], x=[32:40] だと
median ratio が 0.88-0.89 にしかならない。しかも波長（=較正）依存が無く、
IDL と Python が小数第 2 位まで一致するので、原因は **箱の位置** に絞られる。
論文は inter-moss 領域を目視で選んでおり、座標が書かれていない。

指標:
  med   : ratio の中央値。1.0 に近いほど良い
  scat  : log10(ratio) の標準偏差。小さいほど良い（全線が揃って合う）
  slope : log10(ratio) vs logT の傾き。0 に近いほど温度組成が論文と同じ
  score : |log10(med)| + scat  ← これを最小化する

使い方:
    python scripts/scan_boxes.py data/eis/eis_20110702_030712.data.h5 \
           --dy 30 --dx 8 --ystep 15 --xstep 6 --out work/box_scan.csv
"""
import argparse
import sys
import warnings

import numpy as np

sys.path.insert(0, "scripts")
from compare_table2 import fit_box  # noqa: E402

warnings.filterwarnings("ignore")

# 温度（log T [K]、CHIANTI の形成温度のおおよそ）。傾きの計算に使う。
LOGT = {
    "Si VII": 5.80, "Fe IX": 5.90, "Fe X": 6.05, "Fe XI": 6.15,
    "S X": 6.15, "Si X": 6.15, "Fe XII": 6.20, "Fe XIII": 6.25,
    "Fe XIV": 6.30, "Fe XV": 6.35, "S XIII": 6.40, "Fe XVI": 6.45,
    "Ar XIV": 6.50, "Ca XIV": 6.55, "Ca XV": 6.65, "Ca XVI": 6.70,
    "Ca XVII": 6.75,
}


def score_box(rows):
    """Ca XVII は eispac のテンプレートが壊れているので必ず除く。"""
    r = [x for x in rows if x["ion"] != "Ca XVII" and np.isfinite(x["ratio"]) and x["ratio"] > 0]
    if len(r) < 15:
        return None
    ratios = np.array([x["ratio"] for x in r])
    logt = np.array([LOGT.get(x["ion"], 6.2) for x in r])
    med = float(np.median(ratios))
    scat = float(np.std(np.log10(ratios)))
    slope = float(np.polyfit(logt, np.log10(ratios), 1)[0])
    return dict(med=med, scat=scat, slope=slope,
                score=abs(np.log10(med)) + scat, n=len(r))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datafile")
    ap.add_argument("--dy", type=int, default=30)
    ap.add_argument("--dx", type=int, default=8)
    ap.add_argument("--ystep", type=int, default=15)
    ap.add_argument("--xstep", type=int, default=6)
    ap.add_argument("--ymin", type=int, default=120)
    ap.add_argument("--ymax", type=int, default=430)
    ap.add_argument("--xmin", type=int, default=6)
    ap.add_argument("--xmax", type=int, default=52)
    ap.add_argument("--out", default="work/box_scan.csv")
    a = ap.parse_args()

    results = []
    with open(a.out, "w") as f:
        f.write("y0,y1,x0,x1,n,median,scatter,slope,score\n")
        for y0 in range(a.ymin, a.ymax - a.dy + 1, a.ystep):
            for x0 in range(a.xmin, a.xmax - a.dx + 1, a.xstep):
                rows = fit_box(a.datafile, y0, y0 + a.dy, x0, x0 + a.dx)
                s = score_box(rows)
                if s is None:
                    continue
                f.write(f"{y0},{y0+a.dy},{x0},{x0+a.dx},{s['n']},"
                        f"{s['med']:.4f},{s['scat']:.4f},{s['slope']:+.4f},{s['score']:.4f}\n")
                f.flush()
                results.append((s["score"], y0, x0, s))
                print(f"y=[{y0}:{y0+a.dy}] x=[{x0}:{x0+a.dx}]  "
                      f"med={s['med']:.3f} scat={s['scat']:.3f} "
                      f"slope={s['slope']:+.3f} score={s['score']:.3f}", flush=True)

    results.sort()
    print("\n=== 上位 10 箱 ===")
    for sc, y0, x0, s in results[:10]:
        print(f"  y=[{y0}:{y0+a.dy}] x=[{x0}:{x0+a.dx}]  "
              f"med={s['med']:.3f} scat={s['scat']:.3f} slope={s['slope']:+.3f} score={sc:.3f}")


if __name__ == "__main__":
    main()
