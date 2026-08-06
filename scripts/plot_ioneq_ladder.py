"""コロナの電離平衡を描く —— 「なぜ輝線が温度計になるか」の図。

【この図が言いたいこと】

コロナの電離平衡は **密度に依らず温度だけで決まる**。光球の Saha とは違う:

  光球(LTE): 電離 ∝ n_e n_ion  vs  三体再結合 ∝ n_e^2 n_ion+1
             → n_e が残る = 密度が効く
  コロナ    : 電離 ∝ n_e n_ion  vs  輻射/二電子性再結合 ∝ n_e n_ion+1
             → **n_e が約分される** = 温度だけで決まる

結果として各イオンは log T で 0.3-0.5 dex の狭い範囲にしか存在しない。
だから「Fe XII が見える = 1-2 MK のプラズマがある」と言える。

上段: 鉄の電離段が温度で入れ替わる様子
下段: 論文が使う 22 輝線（+ AIA Fe XVIII）のイオンが log T 5.8-6.9 を
      敷き詰めていること（Warren 本人が設計した観測なので当然だが）

    XUVTOP=/opt/ssw/packages/chianti/dbase python scripts/plot_ioneq_ladder.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROM = {7: "VII", 9: "IX", 10: "X", 11: "XI", 12: "XII", 13: "XIII",
       14: "XIV", 15: "XV", 16: "XVI", 17: "XVII", 18: "XVIII"}

# 論文の 22 輝線が使うイオン（+ AIA の Fe XVIII）
USE = [(14, 7, "Si VII 275.4"), (26, 9, "Fe IX 188.5"), (26, 10, "Fe X 184.5"),
       (26, 11, "Fe XI 180.4"), (16, 10, "S X 264.2"), (14, 10, "Si X 258.4"),
       (26, 12, "Fe XII 195.1"), (26, 13, "Fe XIII 202.0"), (26, 14, "Fe XIV 264.8"),
       (26, 15, "Fe XV 284.2"), (16, 13, "S XIII 256.7"), (26, 16, "Fe XVI 263.0"),
       (18, 14, "Ar XIV 194.4"), (20, 14, "Ca XIV 193.9"), (20, 15, "Ca XV 201.0"),
       (20, 16, "Ca XVI 208.6"), (20, 17, "Ca XVII 192.9"),
       (26, 18, "Fe XVIII 93.9 (AIA)")]


def read_ioneq(path):
    lines = open(path).read().split("\n")
    n1, _ = (int(x) for x in lines[0].split())
    logt = np.array([float(x) for x in
                     " ".join(lines[1:1 + n1 // 10 + 1]).split()][:n1])
    data = {}
    for l in lines[1:]:
        p = l.split()
        if len(p) == n1 + 2:
            data[(int(p[0]), int(p[1]))] = np.array([float(x) for x in p[2:]])
    return logt, data


def main():
    xuv = os.environ.get("XUVTOP", "/opt/ssw/packages/chianti/dbase")
    path = os.path.join(xuv, "ioneq", "chianti.ioneq")
    if not os.path.exists(path):
        sys.exit(f"{path} が無い。XUVTOP を設定するか CHIANTI を用意すること")
    logt, data = read_ioneq(path)
    out = sys.argv[1] if len(sys.argv) > 1 else "figures/ioneq_ladder.png"

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 7.2), sharex=True)
    cm = plt.get_cmap("turbo")

    ax = axes[0]
    ions = [9, 10, 11, 12, 13, 14, 15, 16, 18]
    for k, i in enumerate(ions):
        y = data.get((26, i))
        if y is None:
            continue
        c = cm(k / (len(ions) - 1))
        ax.plot(logt, y, lw=2, color=c)
        j = int(np.argmax(y))
        ax.annotate(f"Fe {ROM[i]}", (logt[j], y[j]), fontsize=8,
                    ha="center", va="bottom", color=c)
    ax.set_ylabel(r"ionisation fraction  $N_{\rm ion}/N_{\rm Fe}$")
    ax.set_title("Coronal ionisation equilibrium (CHIANTI 9.0.1)\n"
                 "Each ion exists over only ~0.3-0.5 dex in log T "
                 "-> each line is a thermometer", fontsize=10)
    ax.set_ylim(0, 1.0)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for k, (z, i, lab) in enumerate(USE):
        y = data.get((z, i))
        if y is None:
            continue
        c = cm(k / (len(USE) - 1))
        ax.plot(logt, y, lw=1.8, color=c)
        j = int(np.argmax(y))
        ax.annotate(lab, (logt[j], y[j]), fontsize=6.5, rotation=90,
                    ha="center", va="bottom", color=c)
    ax.set_xlim(5.4, 7.3)
    ax.set_ylim(0, 0.85)
    ax.set_xlabel(r"$\log T$ [K]")
    ax.set_ylabel("ionisation fraction")
    ax.set_title("The 22 lines used by Warren+2012 (+ AIA Fe XVIII): "
                 "they tile log T = 5.8 - 6.9", fontsize=10)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    for z, i, lab in [(26, 12, "Fe XII"), (26, 18, "Fe XVIII"), (20, 17, "Ca XVII")]:
        y = data[(z, i)]
        j = int(np.argmax(y))
        w = logt[y > 0.1 * y.max()]
        print(f"  {lab:9s} peak logT={logt[j]:.2f}  max fraction={y[j]:.3f}  "
              f"width(>10%)={w.min():.2f}-{w.max():.2f} ({w.max()-w.min():.2f} dex)")


if __name__ == "__main__":
    main()
