"""EM loci 図を描く —— DEM を理解する最短路。

【この図が何を表すか】

輝線 λ の観測強度 I_λ は
    I_λ = (1/4π) ∫ G_λ(T) ξ(T) dT
で、ξ(T) が求めたい温度分布（DEM）。

ここで「**もし視線上のプラズマが全部ちょうど温度 T にあったら**、
必要な EM はいくらか」を考える:

    EM_loci,λ(T) = 4π I_λ / G_λ(T)

- これは各温度における **EM の上限**。実際のプラズマは複数の温度に散らばって
  いるので、真の EM は必ずこの曲線より **下** にある。
- 全輝線の曲線を重ねると、**下側の包絡線**が DEM の目安になる。
- **等温プラズマなら全曲線が 1 点で交わる。** 交わらなければ多温度。

→ 逆問題を解く前にこれを描けば「答えの見当」がつく。
   解が変になったときの検算にもなる。

    python scripts/plot_emloci.py [出力png]
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GOFNT = "work/gofnt_chianti901_005.txt"
INTEN = "work/idl_intensities_tied.csv"
CAINT = "work/idl_ca_intensities.csv"
MCMC = "work/mcmc_dem_result.txt"


def read_gofnt(path):
    lines = open(path).readlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# nT nline"))
    nT, nline = (int(v) for v in lines[i + 1].split())
    k = i + 2

    def skip(tag):
        nonlocal k
        while not lines[k].startswith(tag):
            k += 1
        k += 1

    def take(n):
        nonlocal k
        out = []
        while len(out) < n:
            out += [float(x) for x in lines[k].split()]
            k += 1
        return np.array(out[:n])

    skip("# logT")
    logT = take(nT)
    skip("# ion")
    names, wtar = [], []
    for _ in range(nline):
        p = lines[k].split()
        names.append(" ".join(p[:-2]))
        wtar.append(float(p[-2]))
        k += 1
    skip("# G(T)")
    G = np.array([take(nT) for _ in range(nline)])
    return logT, names, np.array(wtar), G


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "figures/emloci_region7.png"
    logT, names, wtar, G = read_gofnt(GOFNT)

    iobs = np.zeros(len(names))
    for r in csv.DictReader(open(INTEN)):
        w = float(r["wvl"])
        k = int(np.argmin(np.abs(wtar - w)))
        if abs(wtar[k] - w) < 0.01:
            iobs[k] = float(r["I_idl"])
    for line in open(CAINT):
        p = [x.strip() for x in line.split(",")]
        try:
            w, ical = float(p[0]), float(p[1])
        except (ValueError, IndexError):
            continue
        if abs(w - 208.604) < 0.01 or abs(w - 192.858) < 0.01:
            iobs[int(np.argmin(np.abs(wtar - w)))] = ical

    ok = np.where(iobs > 0)[0]
    fig, ax = plt.subplots(figsize=(8.4, 5.8))

    cmap = plt.get_cmap("turbo")
    tpk = np.array([logT[int(np.argmax(G[k]))] for k in ok])
    norm = plt.Normalize(tpk.min(), tpk.max())

    env = np.full(len(logT), np.inf)
    for k, tp in zip(ok, tpk):
        g = G[k]
        m = g > g.max() * 1e-3          # G が十分ある温度だけ描く
        loci = 4 * np.pi * iobs[k] / np.where(g > 0, g, np.nan)
        ax.plot(logT[m], loci[m], color=cmap(norm(tp)), lw=1.3, alpha=0.9)
        # ラベルは曲線の最小値（＝その線が最も強く効く温度）に置く
        j = int(np.nanargmin(np.where(m, loci, np.inf)))
        ax.annotate(f"{names[k]} {wtar[k]:.1f}", (logT[j], loci[j]),
                    fontsize=6.5, color=cmap(norm(tp)),
                    xytext=(2, 2), textcoords="offset points")
        env = np.minimum(env, np.where(m, loci, np.inf))

    ax.plot(logT, env, "k--", lw=1.6, label="lower envelope (upper limit on EM)")

    # MCMC で解いた EM を重ねる
    try:
        m = np.loadtxt(MCMC)
        dlt_m = float(np.median(np.diff(m[:, 0])))
        ax.plot(m[:, 0], m[:, 2] * dlt_m, "k-", lw=2.6, marker="o", ms=4,
                label="EM distribution from MCMC")
    except Exception:
        pass

    ax.set_yscale("log")
    ax.set_xlim(5.4, 7.2)
    ax.set_ylim(1e25, 1e31)
    ax.set_xlabel(r"$\log T$ [K]")
    ax.set_ylabel(r"EM$_{\rm loci} = 4\pi I_\lambda / G_\lambda(T)$   [cm$^{-5}$]")
    ax.set_title(
        "EM loci: EM required if ALL plasma were at temperature T\n"
        "Each curve is an upper limit; the true EM lies below. Isothermal -> all cross at one point.",
        fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=ax, label=r"peak of $G(T)$ for that line, $\log T$")
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")

    j = int(np.argmin(env))
    # 包絡線の「最小値」は低温側の弱い線で決まってしまうので、
    # 高温側（logT >= 6.0）で最大になる場所＝プラズマが最も多い温度を見る
    fin = np.isfinite(env)
    jj = int(np.argmin(np.where(fin, env, np.inf)))
    print(f"包絡線の最小: logT = {logT[jj]:.2f} ({10**logT[jj]/1e6:.2f} MK), "
          f"EM = {env[jj]:.3e} cm^-5")
    print("  ※ 包絡線は EM の上限なので、最小になる温度が"
          "「その温度だけでは説明しきれない = 最も強い拘束」を与える")


if __name__ == "__main__":
    main()
