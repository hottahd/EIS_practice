# %% [markdown]
# # Hinode/EIS データ解析講習会
#
# ## この講習会について
#
# Hinode/EISのデータを触り、慣れ親しむことで、SOLAR-C_EUVSTの打ち上げに備えることが目標です。
#
# EUVST は EIS と同じ**スリット走査型の EUV 分光器**で、
# **EUV バンド 170–215 Å は EIS の短波長帯 171–212 Å とほぼ同じ**です。
# ここでの経験はEUVSTの解析にも活かせるでしょう。
#
# | | EUVST | Hinode/EIS (2006–) |
# |---|---|---|
# | 波長帯 | **170–215 Å** + 460–1220 Å | 171–212 / 245–291 Å |
# | 温度被覆 | 2×10⁴ – 1.5×10⁷ K（シームレス） | 飛び飛び |
# | 空間分解能 | **0.4″** | ~2″ |
# | カデンス | **1 秒** | 数十秒 |
# | 実効面積 | **10–30 倍** | — |
#
# ## 今日の内容
#
# | 章 | 内容 |
# |---|---|
# | 1 | EIS のデータを見る |
# | 2 | フィットして**輻射強度**を出す |
# | 3 | **ドップラー速度**を出す |
# | 4 | 線幅から**非熱的速度**を出す |
# | 5 | **温度分布 (Differential Emission Measure: DEM)** を出す |
#
# **eispac（EIS を扱う Python パッケージ）の使い方は、
# [付録 K](#eispac) に早見表としてまとめてあります。**
# 手が止まったらそこを見てください（Colab では左の目次からも飛べます）。
#
# 題材は **2011 年 7 月 2 日 03:07 UT、活動領域 NOAA 1243**。
# Warren, Winebarger & Brooks (2012), ApJ 759, 141 が使ったデータで、
# **論文に観測値の表が載っている**ので、自分の結果と答え合わせができます。

# %% [markdown]
# ## 準備
#
# パッケージを入れて、教材リポジトリとデータを取ってきます。
# **観測データは必要になったところで自動的に取得**されます（既にあれば何もしません）。
#
# インストールで赤い `ERROR:` が出ることがありますが、Colab に元から入っている
# 別のパッケージとの食い違いの報告で、**教材で使うものは正しく入ります**。

# %%
!pip install -q eispac fiasco demregpy

# %%
import os
import subprocess
import sys
from importlib.metadata import version

# pip が numpy を入れ替えた場合だけランタイムを再起動する。
# 入れ替えの直後は、実行中のセッションが古い numpy を掴んだままになり、
# あとで ImportError が出るため（Colab で稀に起きる）。
if "numpy" in sys.modules and sys.modules["numpy"].__version__ != version("numpy"):
    print("numpy が入れ替わったのでランタイムを再起動します。"
          "再起動したら、もう一度先頭から実行してください。")
    import IPython
    if IPython.get_ipython() is not None:
        IPython.get_ipython().kernel.do_shutdown(True)

REPO = "https://github.com/hottahd/EIS_practice.git"
if not os.path.exists("scripts/lines_warren2012.py"):      # リポジトリの外にいる
    if not os.path.exists("EIS_practice"):
        print("教材リポジトリを取得中 ...")
        subprocess.run(["git", "clone", "-q", REPO], check=True)
    os.chdir("EIS_practice")
sys.path.insert(0, "scripts")
print("作業ディレクトリ:", os.getcwd())

# %% [markdown]
# ### 観測データを取る
#
# **NRL のアーカイブ**から level-1 の HDF5 を直接落とします。**ユーザ登録は不要**です。
# URL は日付でできています:
#
#     https://eis.nrl.navy.mil/level1/hdf5/YYYY/MM/DD/eis_YYYYMMDD_HHMMSS.{data,head}.h5
#
# - `.data.h5` (94 MB) … スペクトルの中身
# - `.head.h5` (421 KB) … ヘッダ（波長軸、ポインティング、露出時間など）
#
# **level-1 とは**: CCD の生データ（level-0）に対して、ペデスタルと暗電流を引き、
# 不良画素と宇宙線を除き、実効面積で割って物理単位にしたものです。

# %%
import eispac
import sunpy
import numpy as np

print("eispac", eispac.__version__, " sunpy", sunpy.__version__,
      " numpy", np.__version__)

import os
import urllib.request

BASE = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
EIS_FILE = "data/eis/eis_20110702_030712.data.h5"      # 今日使うデータ

os.makedirs("data/eis", exist_ok=True)
for ext in ("data", "head"):
    name = f"eis_20110702_030712.{ext}.h5"
    path = f"data/eis/{name}"
    if os.path.exists(path) and os.path.getsize(path) > 0:
        print("既にある:", path)
    else:
        print("取得中  :", name, "...")
        urllib.request.urlretrieve(f"{BASE}/{name}", path)

print("\n準備完了:", EIS_FILE)

# %% [markdown]
# **★ Colab の保存について**
#
# GitHub から開いたノートは**読み取り専用の一時セッション**です。
#
# - 編集や実行結果を残したいときは「ファイル → ドライブにコピーを保存」
# - 仮想マシンが切れるとダウンロードしたデータも消えますが、
#   **上から流し直せば復帰します**（既にあるファイルは取得し直しません）
