# %% [markdown]
# # モジュール 1: EIS のデータを見る
#
# **所要時間 40 分**
#
# **このノートで身につくこと**
#
# 1. EIS が撮っているのは**画像ではなくスペクトル**だと体で分かる
# 2. ラスター画像の x 軸が**空間であると同時に時間**であることを確認する
# 3. 輝線を変えると**同じ場所がまったく違う姿に見える**ことを見る
#    → これが講習会全体の動機付け
# 4. 「どの輝線が使えるか」は観測プログラムの設計時に決まっていることを知る
#
# 前提: モジュール 0（環境構築とデータ取得）が終わっていること。

# %% [markdown]
# ## 1-0. 準備
#
# Colab で**このノートから始めた人**もここで動くように、
# 足りないファイルだけ取ってくる（既にあれば何もしない）。

# %%
!pip install -q eispac

# %% [markdown]
# ### Colab のためのおまじない（ローカルで動かしている人は素通りします）
#
# **Colab はノートブック 1 冊ごとに新しい仮想マシンが立ち上がる。**
# GitHub から開いただけでは教材リポジトリもデータも無いので、ここで用意する。
# セッションが切れたときも、このセルをもう一度実行すれば復帰できる。

# %%
import os
import subprocess
import sys

REPO = "https://github.com/hottahd/EIS_practice.git"
if not os.path.exists("scripts/lines_warren2012.py"):      # リポジトリの外にいる
    if not os.path.exists("EIS_practice"):
        print("教材リポジトリを取得中 ...")
        subprocess.run(["git", "clone", "-q", REPO], check=True)
    os.chdir("EIS_practice")
sys.path.insert(0, "scripts")
print("作業ディレクトリ:", os.getcwd())

# %%
import os
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import eispac

EIS_FILE = "data/eis/eis_20110702_030712.data.h5"


def ensure(url, path):
    """path が無ければ url から落とす。あれば何もしない。"""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print("downloading", url)
    urllib.request.urlretrieve(url, path)
    return path


base = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
for ext in ("data", "head"):
    ensure(f"{base}/eis_20110702_030712.{ext}.h5",
           f"data/eis/eis_20110702_030712.{ext}.h5")
print("ok")

# %% [markdown]
# ## 1-1. データを 1 つ読む
#
# `eispac.read_cube(ファイル, 波長)` で、その波長を含む
# **スペクトルウィンドウ**を丸ごと読む。

# %%
cube = eispac.read_cube(EIS_FILE, 195.119)      # Fe XII 195.119 Å
print("shape (y, x, wavelength) =", cube.data.shape)
print("単位                     =", cube.unit)

# %% [markdown]
# **`(512, 60, 24)` が意味するもの**
#
# | 軸 | 数 | 正体 |
# |---|---|---|
# | 0 | 512 | **スリットに沿った空間**（1″/画素） |
# | 1 | 60 | **ラスターのステップ**（2″/画素） |
# | 2 | 24 | **波長**（0.0223 Å/画素） |
#
# EIS は細長いスリット（1″×512″）を太陽に当て、その 1 次元の像を
# 波長分散させて CCD に落とす。**1 回の露出で得られるのは
# (空間 512) × (波長) の 2 次元**であって、画像ではない。
#
# 2 次元の画像がほしければ**スリットを横に振る**。これがラスター。

# %% [markdown]
# ## 1-2. まずスペクトルを 1 本見る
#
# 「EIS はスペクトルを撮っている」を実感するために、
# 1 画素分のスペクトルをそのまま描く。

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
# **見えていること**
#
# - 山が 1 つ。これが Fe XII 195.119 Å。**幅は 3–4 画素しかない**。
# - 台が浮いている。これが背景（連続光＋散乱光）。
# - **輝線強度はこのガウシアンの面積**。次のモジュールでこれを測る。
#
# 山の幅 σ ≈ 0.030 Å は**ほとんどが装置の幅**であって、
# プラズマの熱運動や乱流速度はその上に乗るわずかな超過分。
# だから線幅から速度を出すのは繊細な仕事になる。

# %% [markdown]
# ## 1-3. ラスター画像を作る（フィット無しのクイックルック）
#
# 波長方向にただ足すだけ。連続光や隣の輝線も混ざるので**強度としては不正確**だが、
# 「どこに何があるか」を掴むには数秒で済むこの方法が便利。
#
# **8 つの輝線で同じ場所を見る。** 並べる順は形成温度の順。

# %% [markdown]
# ### ★ その前に: 欠損値の罠
#
# 単純に `np.nansum(cube.data, axis=2)` とやると**横に黒い縞**が出る。
# まずそれを見てから、原因を確かめる。

# %%
img_naive = np.nansum(cube.data, axis=2)        # 素直に足しただけ
plt.figure(figsize=(4, 8))
v = np.sqrt(np.clip(img_naive, 0, None))
plt.imshow(v, origin="lower", aspect="auto", cmap="inferno",
           vmin=0, vmax=np.nanpercentile(v, 99.5))
plt.title("np.nansum only → black stripes")
plt.xlabel("x [pix]"); plt.ylabel("y [pix]")
plt.tight_layout()
plt.show()

# 縞の正体を数字で確かめる
print("NaN の数        :", int(np.isnan(cube.data).sum()))      # → 0 個！
print("データの最小値  :", float(np.nanmin(cube.data)))         # → 大きな負の数
print("cube.mask で落とされるサンプル数:", int(np.asarray(cube.mask).sum()))

# %% [markdown]
# **★ 欠損値は NaN では入っていない。**
#
# EIS の level-1 では、不良画素・宇宙線ヒットで捨てられたサンプルは
# **大きな負のフラグ値**（level-0 の `-100` に較正係数を掛けたもの）として入っている。
# だから
#
# - `np.isnan` / `np.nansum` では**素通りする**
# - 足すと大きな負の数が入り、その行だけ暗くなる → 黒い縞
# - **エラーは一切出ない**
#
# eispac はこれを `cube.mask`（True = 使ってはいけない）に立ててくれるので、
# **必ずこれで落とす**。この先の全モジュールで効いてくる:
#
# - 箱の中で平均するとき、欠損を入れると強度が静かに下がる（モジュール 2）。
#   論文 §3 もわざわざ *"In computing these averaged profiles, missing data are
#   not included"* と書いている。
# - この教材の準備でも、マスクを忘れたまま解析していて、
#   弱い線で数 %（Ca XVII +6.9%、Ca XVI −4.6%、Ar XIV +3.2%）ずれていた。

# %%
def raster_image(datafile, wvl):
    """波長方向に積んで強度マップにする。欠損サンプルは平均に入れない。"""
    c = eispac.read_cube(datafile, wvl)
    d = np.where(np.asarray(c.mask, dtype=bool), np.nan, c.data)
    return np.nanmean(d, axis=2) * d.shape[2]      # 平均 × サンプル数 = 積分値

img = raster_image(EIS_FILE, 195.119)
plt.figure(figsize=(4, 8))
v = np.sqrt(np.clip(img, 0, None))
plt.imshow(v, origin="lower", aspect="auto", cmap="inferno",
           vmin=0, vmax=np.nanpercentile(v, 99.5))
plt.title("with cube.mask → stripes gone")
plt.xlabel("x [pix]"); plt.ylabel("y [pix]")
plt.tight_layout()
plt.show()

# %% [markdown]
# ### 8 つの輝線で並べる

# %%
PANELS = [
    (275.368, "Si VII 275.4",  "0.6 MK  moss"),
    (184.536, "Fe X 184.5",    "1.1 MK"),
    (195.119, "Fe XII 195.1",  "1.6 MK"),
    (202.044, "Fe XIII 202.0", "1.8 MK"),
    (262.984, "Fe XVI 263.0",  "2.8 MK"),
    (193.874, "Ca XIV 193.9",  "3.5 MK"),
    (200.972, "Ca XV 201.0",   "4.5 MK"),
    (192.858, "Ca XVII 192.9", "5.6 MK  (blended)"),
]

ext = cube.meta["extent_arcsec"]        # [x0, x1, y0, y1] （太陽面座標, arcsec）
fig, axes = plt.subplots(1, len(PANELS), figsize=(2.3 * len(PANELS), 9))
for ax, (wvl, label, temp) in zip(axes, PANELS):
    v = np.sqrt(np.clip(raster_image(EIS_FILE, wvl), 0, None))   # 平方根で暗部を持ち上げる
    lo, hi = np.nanpercentile(v, [1, 99.5])
    ax.imshow(v, origin="lower", extent=ext, aspect="equal",
              cmap="inferno", vmin=lo, vmax=hi)
    ax.set_title(f"{label}\n{temp}", fontsize=9)
    ax.set_xlabel("Solar X [″]")
    if ax is not axes[0]:
        ax.set_yticklabels([])
axes[0].set_ylabel("Solar Y [″]")
fig.suptitle("NOAA 1243   2011-07-02 03:07 UT   (same field of view, 8 spectral lines)",
             fontsize=12)
fig.tight_layout(rect=[0, 0.01, 1, 0.965])      # suptitle と重ならないように
plt.show()

# %% [markdown]
# **図の中の英語について**: Colab には日本語フォントが入っていないので、
# 図のラベルを日本語にすると豆腐（□）になる。
# この教材では**図は英語、説明は日本語**で統一している。
# どうしても図に日本語を入れたければ `!pip install japanize-matplotlib` を使う。

# %% [markdown]
# ## 1-4. ★ ここが講習会全体の動機
#
# **同じ場所なのに、輝線を変えると別物に見える。**
#
# | 温度 | 見えるもの |
# |---|---|
# | Si VII (0.6 MK) | **まだら模様** = moss。高温ループの**足元**が遷移層で光っている |
# | Fe X–XIII (1–2 MK) | **細いループ**が何本も見える。周辺部まで広がる |
# | Fe XVI 以上 (2.8 MK–) | ループが消え、**中心部の塊**だけが残る = 活動領域コア |
#
# コロナが単一温度なら、どの輝線で撮っても同じ絵になるはずである。
# そうならないということは、**視線上に色々な温度のプラズマが混ざっている**。
# その温度ごとの量を測るのが **DEM 解析**（モジュール 6, 7）。
#
# **★ Ca XVII のパネルをよく見る。** 5.6 MK の線なのに
# Fe XVI や Ca XV（もっと低温）よりも Fe XII (1.6 MK) に似ていないか？
# → **ブレンドしている**（Fe XI と O V が混ざっている）。モジュール 8 で解く。
#
# **★ 輝線によって残る筋の位置が違う**のにも注意。
# ウィンドウごとに CCD の別の領域を使うので、不良画素の位置も輝線ごとに違う。
# 「1 本の線でうまくいったから全部大丈夫」とはならない。

# %% [markdown]
# ## 1-5. ★ この「画像」は同時刻ではない
#
# ラスターは**スリットを 1 ステップずつ動かして**作る。
# つまり x 軸は空間であると同時に**時間軸**でもある。

# %%
h = cube.meta["index"]
dur = cube.meta["duration"]             # ステップごとの露出時間 [s]
print("観測プログラム :", h["stud_acr"], " (提案者:", h["st_auth"] + ")")
print("開始           :", h["date_obs"])
print("終了           :", h["date_end"])
print("ステップ数     :", h["nraster"])
print("露出時間 [s]   :", f"{dur[0]:.1f}  (1 ステップあたり)")
print()
from astropy.time import Time
t0, t1 = Time(h["date_obs"]), Time(h["date_end"])
total = (t1 - t0).to_value("s")
print(f"全体の所要時間 : {total/60:.1f} 分")
print(f"1 ステップ     : {total/h['nraster']:.1f} 秒  "
      f"(= 露出 {dur[0]:.0f} 秒 + 読み出し等)")
print(f"視野           : {h['fovx']:.1f}″ x {h['fovy']:.1f}″")
print(f"x 方向のサンプリング : {h['fovx']/h['nraster']:.2f}″/step  "
      f"(スリット幅は {h['slit_id']})")

# %% [markdown]
# **62 分かけて撮っている。** 左端と右端では 1 時間離れている。
#
# 帰結（受講者が必ず引っかかるところ）:
#
# - **時間変化する現象には使えない。** フレアやジェットがラスターの途中で
#   起きると、その列だけ別の状態が写る。
# - 「速度マップ」を作っても、それは同時刻の速度場ではない。
# - 一方、活動領域コアの**定常的な**構造を測る今回の目的には問題ない。
#   むしろ 62 分の平均になるので S/N の面では有利。
#
# **★ AIA (12 秒カデンス) と重ねるときは、EIS のどの列がどの時刻かを意識する。**
# 我々は AIA 03:38 UT（ラスターのほぼ中央の時刻）を使う。

# %% [markdown]
# ## 1-6. どの輝線が使えるかは、観測の設計時に決まっている
#
# EIS は 171–212 Å と 245–291 Å を観測できるが、
# **全波長を降ろすとテレメトリが足りない**。
# そこで必要な輝線の周りだけを切り出して降ろす → **スペクトルウィンドウ**。

# %%
wi = eispac.read_wininfo(EIS_FILE.replace(".data.h5", ".head.h5"))
print(f"{'#':>3} {'line_id':<22} {'wvl_min':>9} {'wvl_max':>9} {'幅[Å]':>7}")
for w in wi:
    print(f"{w['iwin']:3d} {str(w['line_id']):<22} "
          f"{w['wvl_min']:9.3f} {w['wvl_max']:9.3f} {w['wvl_max']-w['wvl_min']:7.3f}")
print(f"\n{len(wi)} ウィンドウ")

# %% [markdown]
# **読み方**
#
# - ウィンドウ名（`line_id`）は**目安**であって、そこに入っている輝線が
#   1 本とは限らない。`CA XVII 192.470` のウィンドウは 1.2 Å と広く、
#   Fe XI・O V・Ca XVII が全部入っている（だからブレンドが解ける）。
# - **Ca XIV / Ca XV / Ca XVI / Ca XVII が入っているのが決定的に重要**。
#   これらが 3 MK 以上を拘束する。入っていない観測では
#   「コアの高温プラズマ」の議論ができない。
# - `FE XXIV 255.000` はフレア用。今回は使わない。
#
# 論文 Table 1 の 15 活動領域のうち **region 1 (2010-06-19) には Ca XIV–XVI が無い**。
# 同じ解析ができない例として、モジュール 10 で扱う。
#
# → **自分の科学のために EIS を使うときは、まず
#   「その観測プログラムに必要な線が入っているか」を確認する。**

# %% [markdown]
# ## 1-7. 演習
#
# 1. `PANELS` に自分で輝線を足して描いてみる。
#    使える波長は 1-6 のウィンドウ一覧から選ぶ
#    （例: `186.750` Fe XII、`276.300` Mg V、`278.400` Mg VII）。
# 2. 1-2 のスペクトルを、**moss の上**（例 `y=250, x=10` 付近）と
#    **コアの中**（`y=250, x=35`）で描き比べる。
#    Si VII 275.368 のウィンドウでやると差が大きい。
# 3. `Ca XVII 192.9` のパネルが Fe XII に似ている理由を、
#    1-6 のウィンドウ幅の表から説明してみる。
#
# ## まとめ
#
# - EIS が撮るのは **(空間 512) × (波長)**。画像はスリットを振って自分で作る
# - **x 軸は時間軸でもある**（この観測は 62 分）
# - **輝線ごとに違う温度が見える** → これを定量するのが DEM 解析
# - **使える輝線は観測プログラムで決まっている**
#
# 次（モジュール 2）では、この山にガウシアンを当てて
# **輝線強度という数値**を取り出す。
