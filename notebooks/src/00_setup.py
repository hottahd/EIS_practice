# %% [markdown]
# # モジュール 0: 環境構築とデータ取得
#
# Hinode/EIS データ解析講習会 — Warren, Winebarger & Brooks (2012) を再現する
#
# **このノートで何をするか**
#
# 1. 必要なパッケージを入れる（すべて pip、5 分程度）
# 2. 観測データを取る（EIS 94 MB + AIA 3 MB、**ユーザ登録は不要**）
# 3. 動作確認
#
# **題材**: 2011 年 7 月 2 日 03:07 UT、活動領域 NOAA 1243
# （論文 Table 1 の region 7）。**論文に観測値の表が載っている唯一の天体**なので、
# 自分の解析結果を 1 行ずつ答え合わせできる。

# %% [markdown]
# ## 0-1. パッケージを入れる
#
# | パッケージ | 用途 |
# |---|---|
# | `eispac` | EIS の level-1 HDF5 を読む、輝線フィット |
# | `sunpy` / `aiapy` | AIA/HMI の画像を扱う |
# | `fiasco` | CHIANTI の原子データから寄与関数 G(T) を作る |
# | `demregpy` | 正則化インバージョンで DEM を解く |

# %%
!pip install -q eispac aiapy fiasco demregpy

# %%
import eispac, sunpy, aiapy, fiasco, demregpy
import numpy as np, matplotlib.pyplot as plt
print("eispac  ", eispac.__version__)
print("sunpy   ", sunpy.__version__)
print("aiapy   ", aiapy.__version__)
print("fiasco  ", fiasco.__version__)
print("numpy   ", np.__version__)

# %% [markdown]
# ## 0-2. 教材リポジトリを取ってくる
#
# スクリプトと、IDL/SolarSoft 側で作った**参照データ**（合計 96 KB）が入っている。
#
# 参照データの中身:
#
# | ファイル | 中身 | 何のために |
# |---|---|---|
# | `work/gofnt_chianti901.txt` | 22 輝線の G(T)（0.1 dex） | CHIANTI を落とさなくても DEM が解ける |
# | `work/gofnt_chianti901_005.txt` | 同（0.05 dex） | demregpy 用 |
# | `work/aia94_fe18_response.txt` | AIA 94 の Fe XVIII 応答 | 公式応答は低温線込みで使えない |
# | `work/eis_calcurve_20110702.txt` | 打ち上げ後較正カーブ | 較正の効きを試す |
# | `work/idl_intensities_tied.csv` | IDL で出した輝線強度 | 自分の結果の答え合わせ |
# | `work/mcmc_dem_result.txt` | PINTofALE MCMC の DEM | 手法比較の相手 |

# %%
!git clone -q https://github.com/hottahd/EIS_practice.git
%cd EIS_practice
!ls work/ | head -20

# %% [markdown]
# ## 0-3. EIS のデータを取る
#
# **NRL のアーカイブ**から level-1 HDF5 を直接落とす。ユーザ登録は要らない。
#
#     https://eis.nrl.navy.mil/level1/hdf5/YYYY/MM/DD/eis_YYYYMMDD_HHMMSS.{data,head}.h5
#
# - `.data.h5` (94 MB): スペクトルの中身
# - `.head.h5` (421 KB): ヘッダ（波長軸、ポインティング、露出時間など）
#
# **level-1 とは**: CCD の生データ（level-0）に対して
# ペデスタル・暗電流を引き、不良画素と宇宙線を除き、
# 実効面積で割って物理単位 (erg cm⁻² s⁻¹ sr⁻¹ Å⁻¹) にしたもの。

# %%
import os
os.makedirs("data/eis", exist_ok=True)
base = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
for f in ["eis_20110702_030712.data.h5", "eis_20110702_030712.head.h5"]:
    if not os.path.exists(f"data/eis/{f}"):
        !wget -q -c -P data/eis {base}/{f}
!ls -lh data/eis/

# %% [markdown]
# ## 0-4. AIA の画像を取る
#
# **JSOC の synoptic アーカイブ**を使う。1024×1024（2.4″/画素）、1 枚 1 MB 弱、
# **登録不要**。
#
#     http://jsoc.stanford.edu/data/aia/synoptic/YYYY/MM/DD/HHHH/AIAyyyymmdd_hhmm_wwww.fits
#
# フルディスクの level-1（4096²、1 枚 65 MB）は VSO 経由で取れるが、
# サーバが遅くてタイムアウトしやすい。**Colab では synoptic を使うこと。**
# 15″×23″ の箱なら 2.4″/画素でも 6×10 画素あり、平均値を出すには十分。

# %%
os.makedirs("data/sdo/synoptic", exist_ok=True)
base = "http://jsoc.stanford.edu/data/aia/synoptic/2011/07/02/H0300"
for w in ["0094", "0171", "0193"]:
    f = f"AIA20110702_0338_{w}.fits"
    if not os.path.exists(f"data/sdo/synoptic/{f}"):
        !wget -q -c -P data/sdo/synoptic {base}/{f}
!ls -lh data/sdo/synoptic/

# %% [markdown]
# ## 0-5. 動作確認
#
# EIS のデータを 1 つ読んでみる。Fe XII 195.119 Å は活動領域で最も明るい輝線。

# %%
cube = eispac.read_cube("data/eis/eis_20110702_030712.data.h5", 195.119)
print("データの形 (ny, nx, nwvl) =", cube.data.shape)
print("波長 [Å]:", float(cube.wavelength[0,0,0]), "-", float(cube.wavelength[0,0,-1]))
print("単位:", cube.unit)

# %% [markdown]
# ### この観測がどういうものか
#
# - **512 × 60 × 24** = (スリット方向の画素) × (ラスターのステップ) × (波長の画素)
# - EIS は**スリット分光器**。細長いスリット（1″×512″）を太陽に当て、
#   その 1 次元の像を波長分散させて CCD に落とす。
# - 2 次元の画像がほしいので、**スリットを横に 60 回振る** → これがラスター。
# - **★ 重要**: 1 ステップ約 60 秒なので、**全体で約 62 分かかる**。
#   画像に見えるが**同時刻ではない**。左端と右端で 1 時間離れている。

# %%
h = cube.meta["index"]
print("観測プログラム :", h["stud_acr"])          # スタディの略称
print("提案者         :", h["st_auth"])           # study author
print("開始           :", h["date_obs"])
print("終了           :", h["date_end"])
print("ラスター step数 :", h["nraster"])
print("露出時間 [s]   :", cube.meta["duration"][0])   # ステップごとに入っている
print("スリット       :", h["slit_id"])

# %% [markdown]
# 観測プログラム名は `HPW021_VEL_120x512v1`、study author は `Harry Warren`。
# **HPW = Harry P. Warren** — 論文の著者本人が設計した観測である。
# 「使いたい輝線が入った観測を自分で設計する」ところまで含めて研究になる。
#
# **★ 露出時間はヘッダの `exptime` には入っていない**（`None` が返る）。
# EIS はステップごとに露出時間を持ちうるので、
# eispac は `cube.meta["duration"]`（長さ = ラスターのステップ数）に入れている。
#
# ---
#
# ## 0-6. スペクトルウィンドウを見る
#
# EIS は 171–212 Å と 245–291 Å の 2 バンドを観測できるが、
# 全部を降ろすとテレメトリが足りない。そこで
# **必要な輝線の周りだけを切り出して降ろす** → これが「スペクトルウィンドウ」。
#
# **どの輝線が使えるかは観測プログラムの設計時に決まっている。**
# この観測には 25 個のウィンドウがあり、論文が使う 22 輝線が全部入っている。

# %%
from eispac import read_wininfo
wi = read_wininfo("data/eis/eis_20110702_030712.head.h5")
print(f"{'#':>3} {'line_id':<22} {'wvl_min':>9} {'wvl_max':>9}")
for w in wi:
    print(f"{w['iwin']:3d} {str(w['line_id']):<22} {w['wvl_min']:9.3f} {w['wvl_max']:9.3f}")

# %% [markdown]
# ---
#
# ## まとめ
#
# - パッケージはすべて pip で入る
# - データは**登録不要**で取れる（EIS は NRL、AIA は JSOC synoptic）
# - この観測は **62 分かけて撮った 60 ステップのラスター**
# - **25 個のスペクトルウィンドウ**に論文の 22 輝線が全部入っている
#
# 次（モジュール 1）では、フィットせずに手早く強度マップを作り、
# 「どこを解析するか」を決める。
