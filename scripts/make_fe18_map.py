"""AIA 94/171/193 から Fe XVIII マップを作り、EIS ラスターの視野と重ねて表示する。

Warren et al. (2012) の Figure 1-3 に相当する図を作る。この図を見て
"inter-moss" 領域（Fe XVIII で明るく、171 A の moss が無い場所）を選ぶ。

処理の流れ:
  1. AIA level-1 を level-1.5 に (aiapy: update_pointing + register)
  2. 露光時間で割って DN/s にする
  3. Fe XVIII = I_94 - I_94warm(I_171, I_193)   ← scripts/aia_fe18.py
  4. HMI を AIA と同じ向き・グリッドに合わせる
  5. EIS ラスターの視野を四角で重ねる

【AIA の較正について】
  aiapy.calibrate.degradation() による感度劣化補正は **かけない**。
  論文の Fe XVIII 経験式が劣化補正なしのデータで導かれているため。

使い方:
    python make_fe18_map.py <eis_..._data.h5> [sdoディレクトリ] [出力png]
"""
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


def load_aia(sdodir, wave):
    """AIA level-1 を読み、level-1.5 化して DN/s にする。"""
    from aiapy.calibrate import register, update_pointing
    files = sorted(glob.glob(f"{sdodir}/aia.lev1.{wave}A_*.fits"))
    if not files:
        raise FileNotFoundError(f"AIA {wave} A のファイルが {sdodir} にありません")
    m = sunpy.map.Map(files[0])
    m = update_pointing(m)         # JSOC の master pointing で姿勢を更新
    m = register(m)                # 回転・スケール合わせ (= IDL の aia_prep)
    # DN/s に規格化
    m = sunpy.map.Map(m.data / m.exposure_time.to_value(u.s), m.meta)
    return m


def eis_fov(eisfile):
    """EIS ラスターの視野 (xcen, ycen, fovx, fovy) を arcsec で返す。"""
    import eispac
    cube = eispac.read_cube(eisfile, 195.119)
    idx = cube.meta["index"]
    return idx["xcen"], idx["ycen"], idx["fovx"], idx["fovy"], idx["date_obs"]


def main(eisfile, sdodir="data/sdo", outpng="figures/fe18_region7.png"):
    print("AIA を読み込み・level-1.5 化中 ...")
    m94 = load_aia(sdodir, 94)
    m171 = load_aia(sdodir, 171)
    m193 = load_aia(sdodir, 193)

    xc, yc, fovx, fovy, t_eis = eis_fov(eisfile)
    print(f"EIS 視野: xcen={xc:.1f}\" ycen={yc:.1f}\" {fovx:.0f}\" x {fovy:.0f}\"  @{t_eis}")

    # EIS 視野の周りに切り出す（余白 150"）。
    # 先に切り出してから再投影するとフルディスクを再投影せずに済み、
    # Colab でも数秒で終わる。
    pad = 150 * u.arcsec
    bl = SkyCoord((xc - fovx / 2) * u.arcsec - pad, (yc - fovy / 2) * u.arcsec - pad,
                  frame=m94.coordinate_frame)
    tr = SkyCoord((xc + fovx / 2) * u.arcsec + pad, (yc + fovy / 2) * u.arcsec + pad,
                  frame=m94.coordinate_frame)

    # level-1.5 化しても波長ごとに配列サイズ・姿勢が微妙に違う。
    # 94 A のグリッドを基準にして 171/193 を再投影して揃える。
    # （揃えないと Fe XVIII の引き算がピクセル単位でずれる）
    s94 = m94.submap(bl, top_right=tr)
    print(f"切り出しサイズ: {s94.data.shape}")
    s171 = m171.submap(bl, top_right=tr).reproject_to(s94.wcs)
    s193 = m193.submap(bl, top_right=tr).reproject_to(s94.wcs)

    print("Fe XVIII を分離中 ...")
    s_fe18 = sunpy.map.Map(fe18(s94.data, s171.data, s193.data), s94.meta)

    # HMI: AIA と向きが 180 度違うので回転してから AIA グリッドに合わせる
    hmifiles = sorted(glob.glob(f"{sdodir}/hmi.m_45s*.fits"))
    s_hmi = None
    if hmifiles:
        print("HMI を AIA グリッドに合わせ中 ...")
        h = sunpy.map.Map(hmifiles[0]).rotate(order=3)
        s_hmi = h.submap(bl, top_right=tr).reproject_to(s94.wcs)

    panels = [("HMI B_los [G]", s_hmi, "gray", (-500, 500), False),
              ("AIA 171", s171, "sdoaia171", None, True),
              ("AIA 193", s193, "sdoaia193", None, True),
              ("AIA 94", s94, "sdoaia94", None, True),
              ("AIA Fe XVIII", s_fe18, "inferno", None, True)]
    panels = [p for p in panels if p[1] is not None]

    fig = plt.figure(figsize=(4.2 * len(panels), 9))
    for k, (title, sub, cmap, vlim, sqrt) in enumerate(panels):
        ax = fig.add_subplot(1, len(panels), k + 1, projection=sub)
        d = sub.data
        if sqrt:
            d = np.sqrt(np.clip(d, 0, None))
            vmin, vmax = np.nanpercentile(d, [1, 99.7])
        else:
            vmin, vmax = vlim
        ax.imshow(d, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=11)
        ax.grid(False)

        # EIS ラスターの視野を白枠で表示
        corners = SkyCoord([xc - fovx / 2, xc + fovx / 2, xc + fovx / 2,
                            xc - fovx / 2, xc - fovx / 2] * u.arcsec,
                           [yc - fovy / 2, yc - fovy / 2, yc + fovy / 2,
                            yc + fovy / 2, yc - fovy / 2] * u.arcsec,
                           frame=sub.coordinate_frame)
        ax.plot_coord(corners, color="white", lw=1.2)
        if k > 0:
            ax.coords[1].set_ticklabel_visible(False)
            ax.set_ylabel("")

    fig.suptitle(f"NOAA 1243  {t_eis}   (white box = EIS raster FOV)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outpng, dpi=110)
    print("wrote", outpng)

    # Fe XVIII の総強度（論文 Table 1 の I_hot と比較する量）
    tot = np.nansum(s_fe18.data[s_fe18.data > 2.0])
    print(f"Fe XVIII 総強度 (>2 DN/s のピクセルのみ) = {tot:.3e} DN/s")
    print(f"  論文 Table 1 region 7 の I_hot = 6.18e4 DN/s")
    print("  ※ 論文は活動領域全体で積算しており、積算範囲が違えば値は変わる")


if __name__ == "__main__":
    eisf = sys.argv[1]
    sdod = sys.argv[2] if len(sys.argv) > 2 else "data/sdo"
    outp = sys.argv[3] if len(sys.argv) > 3 else "figures/fe18_region7.png"
    main(eisf, sdod, outp)
