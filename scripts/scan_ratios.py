"""「箱を変えれば Si VII / S XIII / Fe XVI が論文に合う場所があるか」を直接探す。

背景:
  IDL(SSW) と Python(eispac) は 22 輝線すべてで一致した（median 0.89）。
  つまり実装の問題ではない。しかし Si VII 275.368 = 0.39、
  S XIII 256.686 = 0.72、Fe XVI 262.984 = 0.62 が論文より低いまま動かない。

  Fe XII 195.119 に対する比で見ると:
      Si VII/FeXII   論文 0.058  我々 0.024
      Fe XVI/FeXII   論文 0.550  我々 0.362
      S XIII/FeXII   論文 0.403  我々 0.308
  論文の場所は中温に対して冷たい方も熱い方も多い = DEM が広い。
  moss（足元）と高温コアの両方を少し含む箱の特徴。

  前回のスキャンは 30x8 画素の 1 サイズしか試しておらず、
  しかも「全線の合いの良さ」を最適化していたのでこの可能性を潰せていない。
  ここでは **箱の大きさも変え**、上の 3 つの比を直接見る。

速度のため、フィットではなく **窓積分**（波長方向に足すだけ）で比を出す。
比を見るだけなので背景の絶対値以外は十分。候補が出たら本フィットで確かめる。

    python scripts/scan_ratios.py data/eis/eis_20110702_030712.data.h5
"""
import sys

import numpy as np
import eispac

# (ラベル, 波長, 積分に使う波長幅[A], 論文 Table 2 の I)
LINES = [
    ("Si VII 275.368", 275.368, 0.15, 66.85),
    ("Fe XII 195.119", 195.119, 0.15, 1147.35),
    ("S XIII 256.686", 256.686, 0.15, 462.30),
    ("Fe XVI 262.984", 262.984, 0.15, 630.81),
]


def window_sum(datafile, wvl, halfwidth):
    """輝線の窓を積分した強度マップ（背景は窓の両端から線形に引く）。"""
    c = eispac.read_cube(datafile, wvl)
    w = np.asarray(c.wavelength)
    if w.ndim == 3:
        w = np.nanmean(w, axis=(0, 1))
    d = c.data
    core = (w > wvl - halfwidth) & (w < wvl + halfwidth)
    bg = ((w > wvl - 3 * halfwidth) & (w < wvl - 1.6 * halfwidth)) | \
         ((w > wvl + 1.6 * halfwidth) & (w < wvl + 3 * halfwidth))
    if bg.sum() < 2:
        bg = ~core
    dw = float(np.median(np.diff(w)))
    base = np.nanmean(d[:, :, bg], axis=2)
    return (np.nansum(d[:, :, core], axis=2) - base * core.sum()) * dw


def main():
    datafile = sys.argv[1]
    maps = {}
    for lab, wvl, hw, ipap in LINES:
        maps[lab] = window_sum(datafile, wvl, hw)
        print(f"  {lab}: map {maps[lab].shape}", flush=True)

    ny, nx = maps["Fe XII 195.119"].shape
    tgt = {lab: ipap / 1147.35 for lab, _, _, ipap in LINES if lab != "Fe XII 195.119"}
    print("\n論文の比 (対 Fe XII 195.119):",
          "  ".join(f"{k.split()[0]+k.split()[1]}={v:.4f}" for k, v in tgt.items()))

    rows = []
    for dy, dx in [(10, 4), (20, 6), (30, 8), (45, 12), (60, 16), (90, 24)]:
        for y0 in range(100, ny - dy, max(dy // 3, 5)):
            for x0 in range(2, nx - dx, max(dx // 3, 2)):
                sl = (slice(y0, y0 + dy), slice(x0, x0 + dx))
                f12 = np.nanmean(maps["Fe XII 195.119"][sl])
                if not np.isfinite(f12) or f12 <= 0:
                    continue
                r = {lab: np.nanmean(maps[lab][sl]) / f12 for lab in tgt}
                # 3 比を同時にどれだけ再現できるか（log 距離）
                dist = np.sqrt(sum((np.log10(max(r[k], 1e-9) / v)) ** 2 for k, v in tgt.items()) / 3)
                rows.append((dist, dy, dx, y0, x0, r))

    rows.sort(key=lambda t: t[0])
    print("\n=== 3 比を同時に最もよく再現する箱 上位 15 ===")
    print(f"{'dist':>6} {'size':>8} {'y0':>5} {'x0':>4}  " +
          "  ".join(f"{k.split()[0]+k.split()[1]:>12}" for k in tgt))
    for dist, dy, dx, y0, x0, r in rows[:15]:
        print(f"{dist:6.3f} {dy:4d}x{dx:<3d} {y0:5d} {x0:4d}  " +
              "  ".join(f"{r[k]:12.4f}" for k in tgt))

    # Si VII 比だけを最大化する箱（moss にどれだけ寄れば届くか）
    rows.sort(key=lambda t: -t[5]["Si VII 275.368"])
    print("\n=== Si VII/Fe XII が最大の箱 上位 8（論文は 0.0583）===")
    for dist, dy, dx, y0, x0, r in rows[:8]:
        print(f"  {dy:3d}x{dx:<3d} y0={y0:4d} x0={x0:3d}  "
              f"SiVII/FeXII={r['Si VII 275.368']:.4f}  "
              f"FeXVI/FeXII={r['Fe XVI 262.984']:.4f}  "
              f"SXIII/FeXII={r['S XIII 256.686']:.4f}")


if __name__ == "__main__":
    main()
