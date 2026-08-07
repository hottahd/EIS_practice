# %% [markdown]
# # 第 1 章: EIS のデータを見る
#
# **20 分**
#
# EIS が何を観測しているのかを、実際のデータを触って確認します。

# %% [markdown]
# ## 1-1. EIS はスペクトルを撮る。画像は自分で作る
#
# EIS はスリット（1″ × 512″）を太陽に当て、その 1 次元の像を
# 波長分解して CCD に落とします。**1 回の露出で得られるのは
# (空間 512) × (波長) の 2 次元**です。
#
# 2 次元の画像がほしければ**スリットでスキャンします。**。これが**ラスター**です。

# %%
import numpy as np
import matplotlib.pyplot as plt
import eispac

from workshop import EIS_FILE, ensure_eis

ensure_eis()
cube = eispac.read_cube(EIS_FILE, 195.119)      # Fe XII 195.119 Å
print("shape (y, x, wavelength) =", cube.data.shape)
print("単位 =", cube.unit)

# %% [markdown]
# **`(512, 60, 24)`**
#
# | 軸 | 数 | 正体 |
# |---|---|---|
# | 0 | 512 | **スリットに沿った空間**（1″/画素） |
# | 1 | 60 | **ラスターのステップ**（2″/画素） |
# | 2 | 24 | **波長**（0.0223 Å/波長画素） |
#
# `read_cube` は、指定した波長を含む**スペクトルウィンドウ**を丸ごと読みます。

# %% [markdown]
# ## 1-2. スペクトルを 1 本見る

# %%
y, x = 250, 35          # 活動領域コアの中の 1 画素
plt.figure(figsize=(7, 4))
plt.plot(cube.wavelength[y, x, :], cube.data[y, x, :], "o-", ms=4)
plt.axvline(195.119, color="r", ls="--", lw=1, label="Fe XII 195.119")
plt.xlabel("wavelength [Å]")
plt.ylabel(f"intensity [{cube.unit}]")
plt.title(f"EIS spectrum, single pixel (y={y}, x={x})")
plt.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# - 山が 1 つ。これが Fe XII 195.119 Å。**幅は 0.066 Å（FWHM）**で、
#   波長画素（0.0223 Å）にすると 3 つぶんしかありません
# - 山の裾が 0 まで落ちていません。スペクトル全体が一定のレベルの上に乗っており、
#   これが**背景**（連続光と散乱光）です。フィットではガウシアンと一緒に差し引きます
# - **この山の面積が輻射強度**、**中心のずれがドップラー速度**、**幅が非熱的速度**(を含む) —— 第 2〜4 章でやります

# %% [markdown]
# ## 1-3. 輝線を変えると、まったく別の太陽が見える
#
# 波長方向に積めば輻射強度のマップになります（フィットせずに数秒でできる簡易版）。
# **形成温度の順**に 8 枚並べます。

# %%
def raster_image(datafile, wvl):
    """波長方向に積んで輻射強度のマップにする。"""
    c = eispac.read_cube(datafile, wvl)
    # 欠損サンプルは平均に入れない（理由は付録 A）
    d = np.where(np.asarray(c.mask, dtype=bool), np.nan, c.data)
    return np.nanmean(d, axis=2) * d.shape[2]


PANELS = [
    (275.368, "Si VII 275.4",  "0.6 MK  moss"),
    (184.536, "Fe X 184.5",    "1.1 MK"),
    (195.119, "Fe XII 195.1",  "1.6 MK"),
    (202.044, "Fe XIII 202.0", "1.8 MK"),
    (262.984, "Fe XVI 263.0",  "2.8 MK"),
    (193.874, "Ca XIV 193.9",  "3.5 MK"),
    (200.972, "Ca XV 201.0",   "4.5 MK"),
    (192.858, "Ca XVII 192.9", "5.6 MK"),
]

ext = cube.meta["extent_arcsec"]        # [x0, x1, y0, y1]（太陽面座標, arcsec）
fig, axes = plt.subplots(1, len(PANELS), figsize=(2.3 * len(PANELS), 9))
for ax, (wvl, label, temp) in zip(axes, PANELS):
    v = np.sqrt(np.clip(raster_image(EIS_FILE, wvl), 0, None))   # sqrt で暗部を持ち上げる
    lo, hi = np.nanpercentile(v, [1, 99.5])
    ax.imshow(v, origin="lower", extent=ext, aspect="equal",
              cmap="inferno", vmin=lo, vmax=hi)
    ax.set_title(f"{label}\n{temp}", fontsize=9)
    ax.set_xlabel("Solar X [″]")
    if ax is not axes[0]:
        ax.set_yticklabels([])
axes[0].set_ylabel("Solar Y [″]")
fig.suptitle("NOAA 1243   2011-07-02 03:07 UT   (same field of view)", fontsize=12)
fig.tight_layout(rect=[0, 0.01, 1, 0.965])
plt.show()

# %% [markdown]
# **同じ場所なのに、別物に見えます。**
#
# | 温度 | 見えるもの |
# |---|---|
# | Si VII (0.6 MK) | **まだら模様** = moss。高温ループの**足元**が遷移層で光っている |
# | Fe X–XIII (1–2 MK) | **細いループ**が何本も。周辺まで広がる |
# | Fe XVI 以上 (2.8 MK–) | ループが消え、**中心部の塊**だけ = 活動領域コア |
#
# コロナが単一温度なら、どの輝線でも同じ絵になるはずです。
# そうならないのは、**視線上に色々な温度のプラズマが混ざっている**から。
# その量を温度ごとに測るのが **DEM 解析**（第 5 章）です。
#
# **図のラベルが英語なのは**、Colab に日本語フォントが無いためです
# （日本語だと豆腐になります）。

# %% [markdown]
# ## 1-4. この「画像」は同時刻ではない
#
# スリットを 1 ステップずつ動かすので、**x 軸は空間であると同時に時間軸**です。

# %%
from astropy.time import Time

h = cube.meta["index"]
t0, t1 = Time(h["date_obs"]), Time(h["date_end"])
total = (t1 - t0).to_value("s")

print("観測プログラム :", h["stud_acr"], " (提案者:", h["st_auth"] + ")")
print("開始 / 終了    :", h["date_obs"], "/", h["date_end"])
print(f"所要時間       : {total/60:.1f} 分  ({h['nraster']} ステップ)")
print(f"1 ステップ     : {total/h['nraster']:.0f} 秒")
print(f"視野           : {h['fovx']:.0f}″ x {h['fovy']:.0f}″"
      f"  （x は {h['fovx']/h['nraster']:.1f}″/step、スリット幅は {h['slit_id']}）")

# %% [markdown]
# **62 分かかっています。** 左端と右端では 1 時間離れている。
#
# - 時間変化する現象（フレア、ジェット）には使えない
# - ドップラー速度のマップも、同時刻の速度場ではない
# - 一方、定常的な構造を測るなら問題なく、むしろ S/N の面で有利
#
# **Solar-C EUVST はカデンス 1 秒**なので、この制約は大きく緩みます。
# ただし「ラスターは掃いて作る」こと自体は変わりません。

# %% [markdown]
# ## 1-5. 演習
#
# 1. `PANELS` に自分で輝線を足して描く。使える波長は下の一覧から選ぶ
# 2. 1-2 のスペクトルを、**moss の上**（`y=250, x=10` 付近）と
#    **コアの中**（`y=250, x=35`）で描き比べる。Si VII 275.368 でやると差が大きい

# %%
wi = eispac.read_wininfo(EIS_FILE.replace(".data.h5", ".head.h5"))
print(f"この観測に入っているスペクトルウィンドウ（{len(wi)} 個）")
print(f"{'#':>3} {'line_id':<22} {'wvl_min':>9} {'wvl_max':>9}")
for w in wi:
    print(f"{w['iwin']:3d} {str(w['line_id']):<22} {w['wvl_min']:9.3f} {w['wvl_max']:9.3f}")

# %% [markdown]
# **どの輝線が使えるかは、観測プログラムの設計時に決まっています。**
# 全波長を降ろすとテレメトリが足りないので、必要な輝線の周りだけを切り出します。
# この観測は Ca XIV–XVII を含んでいるので、3 MK 以上を測れます。
