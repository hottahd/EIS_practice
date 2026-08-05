"""inter-moss 領域（＝解析する箱）を選ぶ。

Warren et al. (2012) の条件:
  「AIA Fe XVIII で明るく、AIA 171 の moss（足元放射）を含まない場所」
  狙いはループの **頂上付近** の性質を測ること。足元（moss）や冷却中のループが
  混ざると、活動領域全体のモデルが無いと解釈できなくなる。

論文は目で見て手で選んでいる。ここでは同じ判断を数値化して候補を出し、
最後は図を見て人間が決める（自動化しきらないのが正しい。物理的判断だから）。

指標:
    score = median(Fe XVIII) / median(AIA 171)
  Fe XVIII が明るく 171 が暗いほど大きい。

使い方:
    python select_intermoss.py [aia_on_eis_grid.npz]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 箱の大きさの候補 (ny, nx) [EIS ピクセル]
# EIS の y は 1"/pix、この観測の x は 2"/step なので (40,10) は 40" x 20"
BOX_SIZES = [(30, 8), (40, 10), (50, 12), (60, 15)]


def scan(fe, a171, ny, nx, step=4, min_fe=3.0):
    """箱をずらしながら score を計算し、候補を返す。"""
    H, W = fe.shape
    out = []
    for y0 in range(0, H - ny, step):
        for x0 in range(0, W - nx, step):
            f = fe[y0:y0 + ny, x0:x0 + nx]
            m = a171[y0:y0 + ny, x0:x0 + nx]
            if not np.isfinite(f).all() or not np.isfinite(m).all():
                continue
            fmed, mmed = np.median(f), np.median(m)
            if fmed < min_fe:          # Fe XVIII が暗すぎる場所は除外
                continue
            out.append((fmed / mmed, fmed, mmed, y0, y0 + ny, x0, x0 + nx))
    out.sort(reverse=True)
    return out


def main(npz="data/cache/aia_on_eis_grid.npz", outpng="figures/intermoss_region7.png"):
    d = np.load(npz)
    fe, a171, eis = d["fe18"], d["aia171"], d["eis_fe12"]

    print(f"Fe XVIII: median={np.nanmedian(fe):.2f}  max={np.nanmax(fe):.1f} DN/s")
    print(f"AIA 171 : median={np.nanmedian(a171):.0f} DN/s\n")

    best = []
    for ny, nx in BOX_SIZES:
        cands = scan(fe, a171, ny, nx)
        print(f"--- 箱 {ny} x {nx} pix ({ny}\" x {nx*2}\") : 候補 {len(cands)} ---")
        for c in cands[:3]:
            score, fmed, mmed, y0, y1, x0, x1 = c
            print(f"  score={score:6.4f}  FeXVIII={fmed:6.2f}  AIA171={mmed:7.0f}"
                  f"   y=[{y0}:{y1}] x=[{x0}:{x1}]")
        if cands:
            best.append(cands[0])
    print()

    if not best:
        print("候補なし。min_fe を下げること。")
        return

    top = max(best)
    score, fmed, mmed, y0, y1, x0, x1 = top
    print(f"★ 採用候補: y=[{y0}:{y1}] x=[{x0}:{x1}]  "
          f"(FeXVIII={fmed:.2f} DN/s, AIA171={mmed:.0f} DN/s, score={score:.4f})")

    fig, axes = plt.subplots(1, 3, figsize=(13, 9))
    for ax, (img, t) in zip(axes, [(fe, "AIA Fe XVIII"),
                                   (a171, "AIA 171 (moss)"),
                                   (eis, "EIS Fe XII 195.1")]):
        v = np.sqrt(np.clip(img, 0, None))
        lo, hi = np.nanpercentile(v, [1, 99.5])
        ax.imshow(v, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi)
        for c in best:
            _, _, _, yy0, yy1, xx0, xx1 = c
            ax.plot([xx0, xx1, xx1, xx0, xx0], [yy0, yy0, yy1, yy1, yy0],
                    color="cyan", lw=0.8, alpha=0.6)
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color="white", lw=2)
        ax.set_title(t, fontsize=10)
        ax.set_xlabel("EIS x [pix]")
    axes[0].set_ylabel("EIS y [pix]")
    fig.suptitle("inter-moss candidates (white = adopted)")
    fig.tight_layout()
    fig.savefig(outpng, dpi=110)
    print("wrote", outpng)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/cache/aia_on_eis_grid.npz")
