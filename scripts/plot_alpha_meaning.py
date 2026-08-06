"""傾き alpha の意味と、その信頼度を 1 枚で示す図。

【この図が言いたいこと】

1. **alpha はピークより低温側（log T 6.0-6.6）で測る。**
   つまり「4 MK より冷たいプラズマがどれだけ少ないか」を測っている。
   EM 分布は「プラズマが各温度で過ごす時間」を n^2 で重み付けしたものなので、
   alpha が急 = 冷たいプラズマが少ない = **ループが冷える前に再加熱されている**。

2. **論文の alpha=2.9 はナノフレア（低頻度加熱）の予測 2.0-2.3 より急**。
   → 加熱は高頻度、プラズマは熱平衡に近い、という主張になる。

3. **★ ところが我々の測定値 alpha=2.30 はちょうど境界上。**
   同じ観測データから、論文はナノフレアを否定でき、我々はできない。
   差は解析の選択（箱の位置・DEM の手法・平滑化）だけ。

4. 一方 **ピークの 4 MK はどの手法でも動かない**。
   同じ図から読む 2 つの量で、信頼度がまったく違う。

    python scripts/plot_alpha_meaning.py [出力png]
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MCMC = "work/mcmc_dem_result.txt"
DEMREG = "work/demregpy_result.txt"
ALPHA_PAPER = 2.9
ALPHA_NANOFLARE_MAX = 2.3      # Mulu-Moore et al. 2011（論文 p.11 が引く値）


def slope(x, y, a, b):
    k = (x >= a) & (x <= b) & (y > 0)
    return np.polyfit(x[k], np.log10(y[k]), 1)[0]


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "figures/alpha_meaning.png"

    m = np.loadtxt(MCMC)
    lt = m[:, 0]
    dlt = float(np.median(np.diff(lt)))
    em, lo, hi = m[:, 2] * dlt, m[:, 3] * dlt, m[:, 4] * dlt

    al = slope(lt, em, 6.0, 6.6)
    be = -slope(lt, em, 6.6, 7.0)
    ip = int(np.argmax(em))
    pk = em[ip]

    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    ax.axvspan(6.0, 6.6, color="tab:blue", alpha=0.07)
    ax.axvspan(6.6, 7.0, color="tab:red", alpha=0.07)
    ax.fill_between(lt, lo, hi, color="0.82", label="MCMC 50% range")
    ax.plot(lt, em, "k-", lw=2.4, marker="o", ms=4, label="MCMC (PINTofALE)")

    try:
        d = np.loadtxt(DEMREG)
        ald = slope(d[:, 0], d[:, 4], 6.0, 6.6)
        ax.plot(d[:, 0], d[:, 4], "s--", color="tab:red", lw=1.4, ms=3.5,
                label=f"demregpy (regularised), alpha={ald:.2f}")
    except Exception:
        ald = np.nan

    x = np.linspace(6.0, 6.6, 10)
    for a, c, lab in [(ALPHA_PAPER, "tab:green", f"paper  alpha={ALPHA_PAPER}"),
                      (ALPHA_NANOFLARE_MAX, "tab:orange",
                       f"nanoflare sim. max  alpha={ALPHA_NANOFLARE_MAX}")]:
        ax.plot(x, pk * 10 ** (a * (x - lt[ip])), ":", color=c, lw=2, label=lab)
    ax.plot(x, pk * 10 ** (al * (x - lt[ip])), "-", color="k", lw=1.2, alpha=0.6,
            label=f"this work  alpha={al:.2f}")

    ax.set_yscale("log")
    ax.set_xlim(5.6, 7.15)
    ax.set_ylim(1e24, 6e27)
    ax.set_xlabel(r"$\log T$ [K]")
    ax.set_ylabel(r"$\xi(T)\,\Delta\log T$   [cm$^{-5}$]")
    ax.set_title(
        "Warren+2012 region 7: what the slope alpha means, and how much to trust it\n"
        "blue band: alpha is measured BELOW the peak (the cooling tail)   "
        "red band: beta", fontsize=9)
    ax.text(6.62, 3e27, f"peak {10**lt[ip]/1e6:.1f} MK\n(robust)",
            fontsize=8, color="tab:blue")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    print(f"this work : alpha={al:.2f}  beta={be:.2f}  peak logT={lt[ip]:.2f} "
          f"({10**lt[ip]/1e6:.2f} MK)")
    print(f"demregpy  : alpha={ald:.2f}")
    print(f"paper     : alpha={ALPHA_PAPER}   nanoflare sim max: {ALPHA_NANOFLARE_MAX}")


if __name__ == "__main__":
    main()
