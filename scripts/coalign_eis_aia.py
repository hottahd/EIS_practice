"""EIS ラスターと AIA の座標を合わせ、AIA 画像を EIS のピクセル格子に載せ替える。

なぜ必要か:
  EIS のポインティング情報には数 arcsec〜十数 arcsec の系統誤差があることが
  知られている。ヘッダの xcen/ycen をそのまま信じて AIA と重ねると
  inter-moss 領域の選択がずれる。Warren et al. (2012) は箱を「AIA Fe XVIII で
  明るく AIA 171 の moss が無い場所」として選んでいるので、
  この座標合わせの精度がそのまま強度の精度になる。

方法:
  EIS Fe XII 195.119 の強度ラスター（1.6 MK）と AIA 193（1.6 MK）は
  形態がよく似ているので、この 2 枚の相互相関でずれを測る。

出力:
  - EIS 格子に載せ替えた AIA 171 / Fe XVIII （npz）
  - 確認用の図

使い方:
    python coalign_eis_aia.py <eis_..._data.h5> [sdoディレクトリ]
"""
import os
import sys
import glob
import numpy as np
import astropy.units as u
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sunpy.map
from astropy.coordinates import SkyCoord

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from aia_fe18 import fe18
from make_fe18_map import load_aia

CACHE = "data/cache"


def eis_intensity_map(eisfile, wvl=195.119, tmplt_name="fe_12_195_119.2c.template.h5"):
    """EIS の輝線強度マップ（sunpy Map）を作る。結果は FITS にキャッシュする。

    全ラスターのフィットは 1 分半ほどかかるので、2 回目以降はキャッシュを読む。
    """
    import eispac
    os.makedirs(CACHE, exist_ok=True)
    tag = os.path.basename(eisfile).replace(".data.h5", "")
    cache = f"{CACHE}/{tag}_{tmplt_name.split('.')[0]}_int.fits"
    if os.path.exists(cache):
        print(f"キャッシュを使用: {cache}")
        return sunpy.map.Map(cache)

    print(f"EIS 全ラスターをフィット中（初回のみ、1-2 分）: {tmplt_name}")
    tmplt = eispac.read_template(eispac.data.get_fit_template_filepath(tmplt_name))
    cube = eispac.read_cube(eisfile, tmplt.central_wave)
    fit = eispac.fit_spectra(cube, tmplt, ncpu=1, ignore_warnings=True)
    m = fit.get_map(component=0, measurement="intensity")
    m.save(cache, overwrite=True)
    print(f"キャッシュ作成: {cache}")
    return m


def cross_correlate_shift(ref, img, max_shift=30):
    """img を ref に合わせるためのピクセルシフト (dy, dx) を相互相関で求める。

    どちらも同じ形の 2 次元配列。NaN は中央値で埋める。
    整数ピクセル精度（EIS のピクセルは 1" x 1" 相当なのでこれで十分）。
    """
    def prep(a):
        a = np.array(a, float)
        a[~np.isfinite(a)] = np.nanmedian(a)
        # 明るいコアだけで相関が決まらないよう、強度の平方根をとって
        # ダイナミックレンジを圧縮する
        return np.sqrt(np.clip(a, 0, None))

    r, i = prep(ref), prep(img)
    ny, nx = r.shape
    best, bdy, bdx = -np.inf, 0, 0
    for dy in range(-max_shift, max_shift + 1):
        for dx in range(-max_shift, max_shift + 1):
            rs = r[max(0, dy):ny + min(0, dy), max(0, dx):nx + min(0, dx)]
            is_ = i[max(0, -dy):ny + min(0, -dy), max(0, -dx):nx + min(0, -dx)]
            if rs.size < 0.5 * r.size:
                continue
            # 重なり領域ごとに正規化した Pearson 相関を使う。
            # 配列全体で一度だけ正規化すると、重なりが小さいシフトほど
            # 見かけの相関が高くなってしまい、ずれを過大評価する。
            rc, ic = rs - rs.mean(), is_ - is_.mean()
            denom = np.sqrt((rc**2).sum() * (ic**2).sum())
            c = float((rc * ic).sum() / denom) if denom > 0 else -np.inf
            if c > best:
                best, bdy, bdx = c, dy, dx
    return bdy, bdx, best


def main(eisfile, sdodir="data/sdo", outnpz="data/cache/aia_on_eis_grid.npz",
         outpng="figures/coalign_region7.png"):
    os.makedirs(CACHE, exist_ok=True)

    m_eis = eis_intensity_map(eisfile)
    print(f"EIS 強度マップ: {m_eis.data.shape}")

    print("AIA を読み込み中 ...")
    m94, m171, m193 = (load_aia(sdodir, w) for w in (94, 171, 193))

    # AIA を EIS の WCS に載せ替える（EIS の 1x1 arcsec 格子になる）
    print("AIA を EIS 格子へ再投影中 ...")
    a193 = m193.reproject_to(m_eis.wcs)
    a171 = m171.reproject_to(m_eis.wcs)
    a94 = m94.reproject_to(m_eis.wcs)
    fe = fe18(a94.data, a171.data, a193.data)

    # EIS Fe XII 195 と AIA 193 でずれを測る
    dy, dx, cc = cross_correlate_shift(m_eis.data, a193.data, max_shift=25)
    print(f"相互相関で求めたずれ: dy={dy} pix, dx={dx} pix  (相関 {cc:.3f})")
    print(f"  → EIS ヘッダの座標に対して AIA を dy={dy}, dx={dx} ずらすと合う")

    def shift(a, dy, dx):
        out = np.full_like(a, np.nan, dtype=float)
        ny, nx = a.shape
        out[max(0, dy):ny + min(0, dy), max(0, dx):nx + min(0, dx)] = \
            a[max(0, -dy):ny + min(0, -dy), max(0, -dx):nx + min(0, -dx)]
        return out

    a171s, fes, a193s = (shift(v, dy, dx) for v in (a171.data, fe, a193.data))

    np.savez(outnpz, aia171=a171s, aia193=a193s, fe18=fes,
             eis_fe12=m_eis.data, dy=dy, dx=dx)
    print("wrote", outnpz)

    # 確認図
    fig, axes = plt.subplots(1, 4, figsize=(16, 9))
    for ax, (d, t) in zip(axes, [(m_eis.data, "EIS Fe XII 195.1"),
                                 (a193s, "AIA 193 (shifted)"),
                                 (a171s, "AIA 171 (moss tracer)"),
                                 (fes, "AIA Fe XVIII")]):
        v = np.sqrt(np.clip(d, 0, None))
        lo, hi = np.nanpercentile(v, [1, 99.5])
        ax.imshow(v, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi)
        ax.set_title(t, fontsize=10)
        ax.set_xlabel("EIS x [pix]")
    axes[0].set_ylabel("EIS y [pix]")
    fig.suptitle(f"EIS/AIA coalignment  (shift dy={dy}, dx={dx} pix)")
    fig.tight_layout()
    fig.savefig(outpng, dpi=110)
    print("wrote", outpng)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "data/sdo")
