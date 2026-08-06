"""MCMC_DEM の結果を論文 Figure 6-8 と同じ形式でプロットする。

論文は DEM そのものではなく **EM 分布 xi(Te) dTe** を描いている（式 (3)）。
ここでも DEM[cm^-5/logK] に T を掛けたもの（= xi dT に比例）を描く。

    python scripts/plot_dem.py work/mcmc_dem_result.txt figures/dem_region7.png
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 論文 Table 1 region 7
ALPHA_PAPER, BETA_PAPER = 2.9, 9.0
# 論文 Table 2 の輝線が効く温度域（G(T) のピーク、scripts/idl/09_gofnt.pro の出力）
GPEAK = {
    "Si VII": 5.80, "Fe IX": 5.90, "Fe X": 6.05, "Fe XI": 6.10, "S X": 6.20,
    "Si X": 6.15, "Fe XII": 6.20, "Fe XIII": 6.25, "Fe XIV": 6.30, "Fe XV": 6.35,
    "S XIII": 6.40, "Fe XVI": 6.45, "Ar XIV": 6.55, "Ca XIV": 6.55,
    "Ca XV": 6.65, "Ca XVI": 6.70, "Ca XVII": 6.75,
}


def slope(logt, y, t0, t1):
    m = (logt >= t0) & (logt <= t1) & (y > 0)
    if m.sum() < 3:
        return np.nan
    return np.polyfit(logt[m], np.log10(y[m]), 1)[0]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "work/mcmc_dem_result.txt"
    out = sys.argv[2] if len(sys.argv) > 2 else "figures/dem_region7.png"
    d = np.loadtxt(src)
    logt = d[:, 0]
    best, med, q25, q75 = d[:, 1], d[:, 2], d[:, 3], d[:, 4]

    # 論文の EM 分布 = DEM * T
    T = 10 ** logt
    em, emlo, emhi, embest = med * T, q25 * T, q75 * T, best * T

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.fill_between(logt, emlo, emhi, color="0.75", label="MCMC 50% range")
    ax.plot(logt, em, "k-", lw=2, label="MCMC median")
    ax.plot(logt, embest, "r--", lw=1, label=r"best fit (min $\chi^2$)")

    for ion, t in GPEAK.items():
        ax.axvline(t, color="0.9", lw=0.6, zorder=0)

    a = slope(logt, em, 6.0, 6.6)
    b = -slope(logt, em, 6.6, 7.0)
    i = int(np.nanargmax(em))
    ax.set_yscale("log")
    ax.set_xlabel(r"$\log T$ [K]")
    ax.set_ylabel(r"$\xi(T_e)\,dT_e \propto \mathrm{DEM}\times T$  [cm$^{-5}$]")
    ax.set_title(
        f"Warren+2012 region 7 (2011-07-02, NOAA 1243), inter-moss box\n"
        f"peak log T = {logt[i]:.2f} ({T[i]/1e6:.1f} MK),  "
        f"$\\alpha$ = {a:.2f} (paper {ALPHA_PAPER}),  "
        f"$\\beta$ = {b:.2f} (paper {BETA_PAPER})",
        fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    print(f"peak logT={logt[i]:.2f} ({T[i]/1e6:.2f} MK)  alpha={a:.2f} beta={b:.2f}")


if __name__ == "__main__":
    main()
