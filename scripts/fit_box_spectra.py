"""指定した箱（inter-moss 領域）の中で EIS スペクトルを平均してから輝線フィットする。

これは Warren et al. (2012) の Step 5 そのもの:
  「選んだ視野内の EIS データを各スペクトルウィンドウから抜き出して平均し
    （欠損データは平均に入れない）、高 S/N の線プロファイルを作ってから
    単一ガウシアンでフィットする」

全ラスターを 22 輝線分フィットすると 30 分以上かかるが、先に平均して
1 本のプロファイルにしてしまえば 22 本でも数秒で終わる。
講習会（特に Google Colab）ではこの順番が本質的に重要。

使い方:
    python fit_box_spectra.py <eis_..._data.h5> <y0> <y1> <x0> <x1> [出力csv]
"""
import sys
import numpy as np
import eispac

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from lines_warren2012 import LINES, pick_component


def average_spectrum(datafile, wvl, y0, y1, x0, x1):
    """箱の中でスペクトルを平均する。NaN（欠損・不良ピクセル）は平均に含めない。

    返り値: (波長配列, 平均強度, 平均強度の誤差)
    """
    cube = eispac.read_cube(datafile, wvl)

    data = cube.data[y0:y1, x0:x1, :]              # (ny, nx, nwvl)
    errs = cube.uncertainty.array[y0:y1, x0:x1, :]
    wave = cube.wavelength[y0:y1, x0:x1, :]

    # EIS は波長が空間位置ごとに僅かにずれる（スリット傾き・軌道変動）。
    # 平均プロファイルを作るときは共通の波長グリッドに載せ替える必要がある。
    # ここでは箱内の波長ずれが 1 ピクセル未満であることを確認したうえで
    # 単純平均する。ずれが大きい場合は補間が必要。
    wmean = np.nanmean(wave, axis=(0, 1))
    wspread = np.nanmax(np.nanstd(wave, axis=(0, 1)))
    dw = np.median(np.diff(wmean))
    if wspread > 0.3 * abs(dw):
        print(f"  ! 警告: 箱内の波長ずれが大きい (std={wspread:.4f} A, "
              f"1pix={abs(dw):.4f} A)。補間を検討すること")

    good = np.isfinite(data)
    n = good.sum(axis=(0, 1))                       # 各波長で有効なピクセル数
    inten = np.nansum(np.where(good, data, 0), axis=(0, 1)) / np.maximum(n, 1)
    # 平均の誤差 = sqrt(Σσ²)/N
    sig = np.sqrt(np.nansum(np.where(good, errs**2, 0), axis=(0, 1))) / np.maximum(n, 1)

    inten[n == 0] = np.nan
    sig[n == 0] = np.nan
    return wmean, inten, sig, int(np.median(n))


def fit_all_lines(datafile, y0, y1, x0, x1, outcsv=None):
    rows = []
    print(f"箱: y=[{y0}:{y1}] x=[{x0}:{x1}]  ({y1-y0} x {x1-x0} = {(y1-y0)*(x1-x0)} pixels)")
    print(f"{'line':<16s} {'I_fit':>10s} {'sig_fit':>9s} {'I_Warren':>10s} {'ratio':>7s}  "
          f"{'chi2':>8s}  component")
    print("-" * 88)

    for ion, wvl, tmplt_name, i_paper, sig_paper in LINES:
        wave, inten, sig, npix = average_spectrum(datafile, wvl, y0, y1, x0, x1)
        tmplt = eispac.read_template(eispac.data.get_fit_template_filepath(tmplt_name))
        comp, ids = pick_component(tmplt, wvl)

        fit = eispac.fit_spectra(inten, tmplt, wave=wave, errs=sig,
                                 ncpu=1, ignore_warnings=True)
        if fit is None:
            print(f"{ion+' '+f'{wvl:.3f}':<16s} {'FIT FAILED':>10s}")
            continue

        i_fit = float(np.atleast_1d(fit.fit["int"][..., comp]).ravel()[0])
        s_fit = float(np.atleast_1d(fit.fit["err_int"][..., comp]).ravel()[0])
        chi2 = float(np.atleast_1d(fit.fit["chi2"]).ravel()[0])
        ratio = i_fit / i_paper

        rows.append((ion, wvl, i_fit, s_fit, i_paper, sig_paper, ratio, chi2, npix))
        print(f"{ion+' '+f'{wvl:.3f}':<16s} {i_fit:10.2f} {s_fit:9.2f} "
              f"{i_paper:10.2f} {ratio:7.2f}  {chi2:8.1f}  [{comp}] {ids[comp]}")

    if outcsv:
        import csv
        with open(outcsv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ion", "wvl", "I_fit", "sig_fit", "I_Warren2012",
                        "sig_Warren2012", "ratio", "chi2", "npix_used"])
            w.writerows(rows)
        print("\nwrote", outcsv)
    return rows


if __name__ == "__main__":
    f = sys.argv[1]
    y0, y1, x0, x1 = (int(v) for v in sys.argv[2:6])
    out = sys.argv[6] if len(sys.argv) > 6 else None
    fit_all_lines(f, y0, y1, x0, x1, out)
