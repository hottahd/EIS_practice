"""MCMC_DEM（PINTofALE）と demregpy（正則化）の EM 分布を重ねて描く。

講習会モジュール 7 の図。「手法によって DEM がどう変わるか」を一目で見せる。

★ 単位（一度間違えたので明記）
  PINTofALE の DEM は [cm^-5 / logK]  → EM(ビン) = DEM × ΔlogT
  demregpy の DEM は  [cm^-5 / K]     → EM(ビン) = DEM × ΔT = DEM × T ln10 ΔlogT
  ここを揃えないと 1e7 ずれる。

    python scripts/plot_dem_compare.py [出力png]
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MCMC = "work/mcmc_dem_result.txt"
DEMREG = "work/demregpy_result.txt"
ALPHA_PAPER, BETA_PAPER = 2.9, 9.0


def slope(lt, y, t0, t1):
    m = (lt >= t0) & (lt <= t1) & (y > 0)
    return np.polyfit(lt[m], np.log10(y[m]), 1)[0] if m.sum() > 2 else np.nan


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "figures/dem_compare.png"

    m = np.loadtxt(MCMC)
    lt_m = m[:, 0]
    dlt_m = float(np.median(np.diff(lt_m)))
    em_m = m[:, 2] * dlt_m                      # DEM[/logK] × ΔlogT
    lo_m = m[:, 3] * dlt_m
    hi_m = m[:, 4] * dlt_m

    d = np.loadtxt(DEMREG)
    lt_d = d[:, 0]
    em_d = d[:, 4]                              # すでに EM(ビン) で保存してある
    e_d = d[:, 2] * 10 ** lt_d * np.log(10) * float(np.median(np.diff(lt_d)))

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.fill_between(lt_m, lo_m, hi_m, color="0.8", label="MCMC 50% range")
    ax.plot(lt_m, em_m, "k-", lw=2.2, marker="o", ms=3.5,
            label="PINTofALE MCMC_DEM")
    ax.errorbar(lt_d, em_d, yerr=e_d, fmt="s-", color="tab:red", lw=1.4, ms=3.5,
                capsize=2, label="demregpy (regularised)")

    a_m, b_m = slope(lt_m, em_m, 6.0, 6.6), -slope(lt_m, em_m, 6.6, 7.0)
    a_d, b_d = slope(lt_d, em_d, 6.0, 6.6), -slope(lt_d, em_d, 6.6, 7.0)
    ip_m = int(np.nanargmax(em_m))
    ip_d = int(np.nanargmax(em_d))

    ax.axvline(6.60, color="tab:blue", ls=":", lw=1)
    ax.text(6.62, ax.get_ylim()[0], " 4 MK", color="tab:blue", fontsize=8, va="bottom")

    ax.set_yscale("log")
    ax.set_xlabel(r"$\log T$ [K]")
    ax.set_ylabel(r"$\xi(T_e)\,dT_e$  [cm$^{-5}$]")
    ax.set_title(
        "Warren+2012 region 7, inter-moss box:  MCMC vs regularised inversion\n"
        f"peak  MCMC {lt_m[ip_m]:.2f} ({10**lt_m[ip_m]/1e6:.1f} MK)   "
        f"demreg {lt_d[ip_d]:.2f} ({10**lt_d[ip_d]/1e6:.1f} MK)   |   "
        rf"$\alpha$  MCMC {a_m:.2f}  demreg {a_d:.2f}  (paper {ALPHA_PAPER})",
        fontsize=9)
    ax.legend(fontsize=8, loc="lower center")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    print(f"MCMC    : peak logT={lt_m[ip_m]:.2f}  alpha={a_m:.2f}  beta={b_m:.2f}")
    print(f"demregpy: peak logT={lt_d[ip_d]:.2f}  alpha={a_d:.2f}  beta={b_d:.2f}")
    print(f"paper   : alpha={ALPHA_PAPER}  beta={BETA_PAPER}")


if __name__ == "__main__":
    main()
