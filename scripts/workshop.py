"""講習会ノートブックの共通ヘルパ。

**なぜ要るか**: Colab はノートブック 1 冊ごとに新しい VM が立ち上がる。
モジュール 2 が書いた `work/box_intensities.csv` も、モジュール 4 が書いた
`data/cache/aia_on_eis_grid.npz` も、**次のノートを開いた時点では存在しない**。

受講者が途中のモジュールから始めても、あるいは Colab のセッションが切れても
続きができるように、「無ければその場で作る」関数をここにまとめる。
どれも数秒〜十数秒で終わる。

各関数の中身は、対応するモジュールのノートで**表示して説明している処理そのもの**
（ノート側は教材なので式や罠を見せながら書き、こちらは再利用のために畳んである）。
"""
import os
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EIS_FILE = "data/eis/eis_20110702_030712.data.h5"
EIS_BASE = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
AIA_DIR = "data/sdo/synoptic"
AIA_BASE = "http://jsoc.stanford.edu/data/aia/synoptic/2011/07/02/H0300"

# 論文 Table 2 の region 7 で採用した inter-moss 箱（モジュール 4 で選ぶ）
BOX = dict(y0=244, y1=274, x0=32, x1=40)

# AIA 94 A の Fe XVIII 分離（論文 Appendix。指数の誤植は修正済み。モジュール 3）
_A = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]


def download(url, path):
    """path が無ければ url から落とす。あれば何もしない。"""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"  downloading {os.path.basename(path)} ...")
    urllib.request.urlretrieve(url, path)
    return path


def ensure_eis():
    """EIS level-1 HDF5（94 MB + 421 KB）を用意する。"""
    for ext in ("data", "head"):
        download(f"{EIS_BASE}/eis_20110702_030712.{ext}.h5",
                 f"data/eis/eis_20110702_030712.{ext}.h5")
    return EIS_FILE


def ensure_aia(wave):
    """AIA synoptic（1024x1024, 2.4"/px, 登録不要）を用意して DN/s の Map を返す。"""
    import sunpy.map
    f = download(f"{AIA_BASE}/AIA20110702_0338_{wave:04d}.fits",
                 f"{AIA_DIR}/AIA20110702_0338_{wave:04d}.fits")
    m = sunpy.map.Map(f)
    return sunpy.map.Map(m.data / m.meta["exptime"], m.meta)


def fe18(i94, i171, i193, f=0.31):
    """AIA 94 から Fe XVIII 成分 [DN/s] を取り出す（モジュール 3）。"""
    x = np.clip((f * np.asarray(i171, float) + (1 - f) * np.asarray(i193, float))
                / 116.54, 0.0, 30.0)
    return np.asarray(i94, float) - 0.39 * sum(a * x**i for i, a in enumerate(_A))


def box_intensities(path="work/box_intensities.csv", box=None, force=False):
    """箱平均 → 22 輝線フィット（モジュール 2）。CSV があればそれを読む。

    返り値: [(ion, wvl, I_fit, I_paper, ratio), ...]
    """
    import csv
    if os.path.exists(path) and not force:
        with open(path) as f:
            return [(r["ion"], float(r["wvl"]), float(r["I_fit"]),
                     float(r["I_paper"]), float(r["ratio"]))
                    for r in csv.DictReader(f)]

    print("work/box_intensities.csv が無いので 22 輝線をフィットする（10 秒ほど）")
    import eispac
    import sys
    sys.path.insert(0, HERE)
    from lines_warren2012 import LINES, pick_component
    from fit_box_spectra import average_spectrum

    ensure_eis()
    box = box or BOX
    rows = []
    for ion, wvl, tname, i_paper, _ in LINES:
        w, I, s, _n = average_spectrum(EIS_FILE, wvl, **box)
        t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
        comp, _ids = pick_component(t, wvl)
        fit = eispac.fit_spectra(I, t, wave=w, errs=s, ncpu=1, ignore_warnings=True)
        i_fit = float(np.atleast_1d(fit.fit["int"][..., comp]).ravel()[0])
        rows.append((ion, wvl, i_fit, i_paper, i_fit / i_paper))

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w_ = csv.writer(f)
        w_.writerow(["ion", "wvl", "I_fit", "I_paper", "ratio"])
        w_.writerows(rows)
    print(f"wrote {path}")
    return rows


def eis_raster_map(wvl=195.119, datafile=None):
    """EIS のウィンドウを波長方向に積んで sunpy Map にする（モジュール 4）。

    フィットは不要。欠損は cube.mask で落とす（大きな負のフラグ値なので
    np.isfinite では防げない。モジュール 1 参照）。
    """
    import astropy.units as u
    import eispac
    import sunpy.map
    from astropy.coordinates import SkyCoord

    datafile = datafile or ensure_eis()
    c = eispac.read_cube(datafile, wvl)
    d = np.where(np.asarray(c.mask, dtype=bool), np.nan, c.data)
    img = np.nanmean(d, axis=2) * d.shape[2]

    h, p = c.meta["index"], c.meta["pointing"]
    ref = SkyCoord(p["xcen"] * u.arcsec, p["ycen"] * u.arcsec,
                   obstime=h["date_obs"], observer="earth", frame="helioprojective")
    hdr = sunpy.map.make_fitswcs_header(
        img, ref, scale=[p["x_scale"], p["y_scale"]] * u.arcsec / u.pix,
        instrument="EIS", wavelength=wvl * u.angstrom)
    hdr["measrmnt"] = "intensity"        # eispac の EISMap が要求するキー
    return sunpy.map.Map(img, hdr)


def cross_correlate_shift(ref, img, max_shift=25):
    """img を ref に合わせる (dy, dx) を整数画素で求める（モジュール 4）。

    重なり領域ごとに正規化した Pearson 相関を使うこと。全体で一度だけ
    正規化すると、重なりの小さいシフトほど見かけの相関が上がる。
    """
    def prep(a):
        a = np.array(a, float)
        a[~np.isfinite(a)] = np.nanmedian(a)
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
            rc, ic = rs - rs.mean(), is_ - is_.mean()
            den = np.sqrt((rc**2).sum() * (ic**2).sum())
            c = float((rc * ic).sum() / den) if den > 0 else -np.inf
            if c > best:
                best, bdy, bdx = c, dy, dx
    return bdy, bdx, best


def aia_on_eis_grid(path="data/cache/aia_on_eis_grid.npz", force=False):
    """AIA を EIS 格子に載せ替えて座標合わせしたもの（モジュール 4）。

    npz があればそれを読む。無ければ作る（AIA 3 MB のダウンロード込みで十数秒）。
    """
    if os.path.exists(path) and not force:
        return dict(np.load(path))

    print("data/cache/aia_on_eis_grid.npz が無いので作る（十数秒）")
    m_eis = eis_raster_map()
    a94, a171, a193 = (ensure_aia(w) for w in (94, 171, 193))
    r94, r171, r193 = (m.reproject_to(m_eis.wcs) for m in (a94, a171, a193))
    fe = fe18(r94.data, r171.data, r193.data)

    dy, dx, cc = cross_correlate_shift(m_eis.data, r193.data)
    print(f"  相互相関のずれ: dy={dy}, dx={dx} pix (相関 {cc:.3f})")

    def shift(a):
        out = np.full_like(np.asarray(a, float), np.nan)
        ny, nx = a.shape
        out[max(0, dy):ny + min(0, dy), max(0, dx):nx + min(0, dx)] = \
            a[max(0, -dy):ny + min(0, -dy), max(0, -dx):nx + min(0, -dx)]
        return out

    out = dict(aia171=shift(r171.data), aia193=shift(r193.data), fe18=shift(fe),
               eis_fe12=m_eis.data, dy=dy, dx=dx)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **out)
    print(f"wrote {path}")
    return out


# --- ラスターの一部をフィットする（第 3・4 章で共有する）------------------
_FIT_CACHE = {}


def fit_region(wvl=195.119, tmplt_name="fe_12_195_119.2c.template.h5",
               y0=180, y1=340, ncpu=2):
    """ラスターの y=[y0:y1] を丸ごとフィットして結果を返す（結果は使い回す）。

    速度（第 3 章）と線幅（第 4 章）は**同じフィット**から出るので、
    2 度フィットしないようにキャッシュする。
    """
    key = (wvl, y0, y1)
    if key in _FIT_CACHE:
        return _FIT_CACHE[key]

    import eispac
    ensure_eis()
    tmplt = eispac.read_template(eispac.data.get_fit_template_filepath(tmplt_name))
    cube = eispac.read_cube(EIS_FILE, tmplt.central_wave)
    print(f"y=[{y0}:{y1}] をフィット中（{(y1-y0)*cube.data.shape[1]} スペクトル）...")
    fit = eispac.fit_spectra(cube[y0:y1, :, :], tmplt, ncpu=ncpu,
                             ignore_warnings=True)
    _FIT_CACHE[key] = (fit, cube)
    return fit, cube
