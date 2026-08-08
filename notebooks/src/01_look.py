# %% [markdown]
# # 第 1 章: EIS のデータを見る
#
# **20 分**
#
# EIS が何を観測しているのかを、実際のデータを触って確認します。

# %% [markdown]
# ## 1-1. スペクトル観測、画像生成
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

from workshop import EIS_FILE, ensure_eis   # データのパスと、無ければ落とす関数

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
# | 2 | 24 | **波長**（0.0223 Å 刻み） |
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
# - 山が 1 つ。これが Fe XII 195.119 Å。幅は **0.066 Å（FWHM）**で、
#   これはほぼ EIS の**波長分解能**（195 Å で 0.059 Å）です。
#   記録は 0.0223 Å 刻みなので、山は 3 点ほどしかありません
# - 山の裾が 0 まで落ちていません。スペクトル全体が一定のレベルの上に乗っており、
#   これが**背景**（連続光と散乱光）です。フィットではガウシアンと一緒に差し引きます
# - **この山の面積が輻射強度**、**中心のずれがドップラー速度**、**幅が非熱的速度**(を含む) —— 第 2〜4 章でやります

# %% [markdown]
# ## 1-3. 輝線ごとの観測データ
#
# 各画素のスペクトルを**波長方向に足し合わせて波長刻みを掛ける**と、
# そのスペクトルウィンドウに入っている放射を**全部含んだ輻射強度**になります。
# フィット（第 2 章）で得られるのは、そのうち**目的の輝線だけ**の輻射強度です。
#
# ここでは「どこに何があるか」が分かればよいので、連続光や隣の輝線が
# 混ざったままで構いません。フィットしないぶん数秒で済みます。
# **形成温度の順**に 8 枚並べます。

# %%
def raster_image(datafile, wvl):
    """スペクトルを波長方向に足し合わせて、輻射強度のマップにする。

    データは波長 1 Å あたりの値なので、足したあとに波長刻みを掛ける
    （そうして初めて erg cm^-2 s^-1 sr^-1 になる）。
    """
    c = eispac.read_cube(datafile, wvl)
    # 欠損サンプル（c.mask が True）を NaN に置き換える。
    # こうしておくと下の np.nanmean が無視してくれる。
    # 元データでは欠損は「大きな負の値」なので、そのままだと平均に入ってしまう（付録 A）
    d = np.where(np.asarray(c.mask, dtype=bool), np.nan, c.data)
    total = np.nanmean(d, axis=2) * d.shape[2]      # 欠損を除いた平均 × サンプル数
    dwave = float(np.median(np.diff(np.asarray(c.wavelength, float), axis=2)))
    return total * dwave


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

# %% [markdown]
# ## 1-4. 演習
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
