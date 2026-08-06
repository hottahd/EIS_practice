"""論文 §2 の順番どおり「1 画素ずつフィット → 箱内で強度を平均」する。

Warren+2012 §2:
    "Intensities were then determined for each emission line of interest at
     every spatial pixel using Gaussian fits to the line profiles."

これまで我々は「箱内でスペクトルを平均 → 1 回フィット」していた（論文 §2 の
別の箇所にある "averaged them together to create high signal-to-noise line
profiles" に引きずられた）。**順番が逆**だと弱い輝線で系統差が出る:

  - 1 画素ずつフィットすると振幅に >=0 の拘束がかかるので、ノイズに
    埋もれた画素でも正の値が出る（ノイズの整流）。平均は上に偏る。
  - 先にスペクトルを平均すると S/N が上がってその偏りが消える。

→ 箱の中で最も弱い線が最も影響を受ける。この箱では Si VII 275.368 が該当。

使い方:
    python scripts/fit_perpixel_box.py data/eis/eis_..._data.h5 244 274 32 40
"""
import sys

import numpy as np

sys.path.insert(0, "scripts")
from lines_warren2012 import LINES, pick_component  # noqa: E402
from fit_box_spectra import average_spectrum  # noqa: E402
import eispac  # noqa: E402


def main():
    datafile = sys.argv[1]
    y0, y1, x0, x1 = (int(v) for v in sys.argv[2:6])
    ncpu = int(sys.argv[6]) if len(sys.argv) > 6 else 8

    print(f"box y=[{y0}:{y1}] x=[{x0}:{x1}]   ({y1-y0} x {x1-x0} = {(y1-y0)*(x1-x0)} px)")
    print(f"{'line':<16s} {'A:平均→fit':>12s} {'B:fit→平均':>12s} {'B/A':>6s} "
          f"{'I_paper':>9s} {'A/pap':>6s} {'B/pap':>6s}")

    rows = []
    for ion, wvl, tname, i_paper, sig_paper in LINES:
        t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
        comp, _ = pick_component(t, wvl)

        # --- A: 箱内平均スペクトルを 1 回フィット（これまでのやり方）
        wave, inten, sig, _ = average_spectrum(datafile, wvl, y0, y1, x0, x1)
        fa = eispac.fit_spectra(inten, t, wave=wave, errs=sig, ncpu=1,
                                ignore_warnings=True)
        i_a = float(np.atleast_1d(fa.fit["int"][..., comp]).ravel()[0]) if fa else np.nan

        # --- B: 1 画素ずつフィットして強度を平均（論文の順番）
        cube = eispac.read_cube(datafile, wvl)
        sub = cube[y0:y1, x0:x1]
        fb = eispac.fit_spectra(sub, t, ncpu=ncpu, ignore_warnings=True)
        if fb is None:
            i_b = np.nan
        else:
            arr = np.asarray(fb.fit["int"][..., comp], dtype=float)
            chi2 = np.asarray(fb.fit.get("chi2", np.ones_like(arr)), dtype=float)
            good = np.isfinite(arr) & (arr > -1e29) & np.isfinite(chi2)
            i_b = float(np.nanmean(arr[good])) if good.any() else np.nan

        rows.append((ion, wvl, i_a, i_b, i_paper))
        ratio = i_b / i_a if (np.isfinite(i_a) and i_a != 0) else np.nan
        print(f"{ion+' '+format(wvl,'.3f'):<16s} {i_a:12.2f} {i_b:12.2f} {ratio:6.2f} "
              f"{i_paper:9.2f} {i_a/i_paper:6.2f} {i_b/i_paper:6.2f}")

    ra = np.array([r[2] / r[4] for r in rows if r[0] != "Ca XVII"])
    rb = np.array([r[3] / r[4] for r in rows if r[0] != "Ca XVII"])
    ra, rb = ra[np.isfinite(ra)], rb[np.isfinite(rb)]
    print()
    print(f"A (平均→fit): median={np.median(ra):.3f}  15%以内 "
          f"{int(((ra>=0.85)&(ra<=1.15)).sum())}/{len(ra)}  最大ずれ "
          f"{max(max(x,1/x) for x in ra if x>0):.2f} 倍")
    print(f"B (fit→平均): median={np.median(rb):.3f}  15%以内 "
          f"{int(((rb>=0.85)&(rb<=1.15)).sum())}/{len(rb)}  最大ずれ "
          f"{max(max(x,1/x) for x in rb if x>0):.2f} 倍")


if __name__ == "__main__":
    main()
