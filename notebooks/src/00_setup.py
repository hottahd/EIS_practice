# %% [markdown]
# # Hinode/EIS データ解析講習会
#
# ## この講習会について
#
# **狙いは Solar-C (EUVST) の準備**です。打ち上げたその日から解析できるように、
# いま手に入るデータで実力をつけておきます。
#
# EUVST は EIS と同じ**スリット走査型の EUV 分光器**で、
# **EUV バンド 170–215 Å は EIS の短波長帯 171–212 Å とほぼ同じ**です。
# つまり**同じ輝線を撮ります**。今日フィットする Fe XII や Ca XV は、
# そのまま EUVST の主力になります。
#
# | | EUVST (2028 打ち上げ予定) | Hinode/EIS (2006–) |
# |---|---|---|
# | 波長帯 | **170–215 Å** + 460–1220 Å | 171–212 / 245–291 Å |
# | 温度被覆 | 2×10⁴ – 1.5×10⁷ K（シームレス） | 飛び飛び |
# | 空間分解能 | **0.4″** | ~2″ |
# | カデンス | **1 秒** | 数十秒 |
# | 実効面積 | **10–30 倍** | — |
#
# **今日やることは、そのまま 2028 年に使えます。**
#
# ## 今日の内容
#
# | 章 | 内容 |
# |---|---|
# | 1 | EIS のデータを見る |
# | 2 | フィットして**強度**を出す |
# | 3 | **速度**を出す |
# | 4 | 線幅から**非熱的速度**を出す |
# | 5 | **温度分布 (DEM)** を出す |
#
# 題材は **2011 年 7 月 2 日 03:07 UT、活動領域 NOAA 1243**。
# Warren, Winebarger & Brooks (2012), ApJ 759, 141 が使ったデータで、
# **論文に観測値の表が載っている**ので、自分の結果と答え合わせができます。

# %% [markdown]
# ## 準備
#
# パッケージを入れて、教材リポジトリとデータを取ってきます。
# **観測データは必要になったところで自動的に取得**されます（既にあれば何もしません）。

# %%
!pip install -q eispac fiasco demregpy

# %% [markdown]
# ### 赤い `ERROR:` が出ても、たいていは無視してよい
#
# Colab では `google-colab 1.0.0 requires requests==2.32.4, but you have ...`
# のような行が出ることがあります。**インストールの失敗ではなく**、
# Colab に元から入っている別のパッケージとの食い違いの報告です。
# 教材で使うものは正しく入っています。
#
# 次のセルは、`numpy` が入れ替わった場合だけランタイムを再起動します
# （再起動したら、もう一度先頭から実行してください。2 回目は一瞬で終わります）。

# %%
import sys
from importlib.metadata import version

need_restart = False
try:
    loaded = sys.modules["numpy"].__version__ if "numpy" in sys.modules else None
    if loaded is not None and loaded != version("numpy"):
        need_restart = True
        print(f"numpy が {loaded} -> {version('numpy')} に入れ替わりました")
except Exception as e:
    need_restart = True
    print("numpy の状態を確認できませんでした:", e)

if need_restart:
    print("ランタイムを再起動します。"
          "再起動したら、もう一度このノートを先頭から実行してください。")
    try:
        import IPython
        ipy = IPython.get_ipython()
        if ipy is not None:
            ipy.kernel.do_shutdown(True)
    except Exception:
        import os
        os.kill(os.getpid(), 9)
else:
    print("numpy の入れ替えは起きていません。このまま先へ進んで大丈夫です。")

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
import eispac
import sunpy
import numpy as np

from workshop import ensure_eis

print("eispac", eispac.__version__, " sunpy", sunpy.__version__,
      " numpy", np.__version__)
ensure_eis()          # EIS の level-1 データ（94 MB）。既にあれば何もしない
print("準備完了")

# %% [markdown]
# **★ Colab の保存について**
#
# GitHub から開いたノートは**読み取り専用の一時セッション**です。
#
# - 編集や実行結果を残したいときは「ファイル → ドライブにコピーを保存」
# - 仮想マシンが切れるとダウンロードしたデータも消えますが、
#   **上から流し直せば復帰します**（既にあるファイルは取得し直しません）
