"""講習会用の Colab ノートブックを生成する。

.ipynb を直接手で書くと JSON が読みにくく、レビューも差分も辛い。
ここでは **セルの中身を素の Python / Markdown として書き**、
このスクリプトで .ipynb に変換する。

  - ノートブックの定義は notebooks/src/NN_*.py に**モジュールごとに**置く（後述の書式）
  - `python notebooks/build_notebooks.py` で、それらを繋いだ
    **1 冊の `EIS_workshop.ipynb`** を生成する
  - コードセルの中身は普通の Python なので、単体で実行して検証できる

★ 配るのは 1 冊だけ。Colab は**ノートブック 1 冊ごとに新しい VM** が立ち上がるので、
  分冊にすると pip install と 94 MB のデータ取得を冊数ぶん繰り返すことになる。
  編集の都合でソースはモジュールごとに分けたままにしてある（差分が読みやすいため）。

書式（notebooks/src/*.py）:

    # %% [markdown]
    # # 見出し
    # 本文...

    # %%
    print("コードセル")

生成される .ipynb は Colab で開ける（GitHub の URL を
https://colab.research.google.com/github/<owner>/<repo>/blob/main/... に貼る）。
"""
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
SRC = HERE / "src"


def parse(path):
    """# %% 区切りのファイルをセルの列に分解する。"""
    cells = []
    kind, buf = None, []

    def flush():
        if kind is None:
            return
        src = "\n".join(buf).strip("\n")
        if not src.strip():
            return
        if kind == "markdown":
            # 各行の先頭の "# " を剥がす
            lines = [l[2:] if l.startswith("# ") else (l[1:] if l == "#" else l)
                     for l in src.split("\n")]
            cells.append({"cell_type": "markdown", "metadata": {},
                          "source": "\n".join(lines)})
        else:
            cells.append({"cell_type": "code", "metadata": {},
                          "execution_count": None, "outputs": [],
                          "source": src})

    for line in path.read_text().split("\n"):
        if line.startswith("# %% [markdown]"):
            flush()
            kind, buf = "markdown", []
        elif line.startswith("# %%"):
            flush()
            kind, buf = "code", []
        else:
            buf.append(line)
    flush()
    return cells


def notebook(cells):
    return {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }


# --- 各モジュールを繋いで 1 冊にする ----------------------------------------
#
# 各モジュールの先頭にある「インストール」「ブートストラップ」のセルは
# **先頭に 1 組だけ**置き、2 回目以降は落とす。判定は下の関数による。
BOOTSTRAP_MD = "### Colab のためのおまじない"

PREAMBLE_MD = """## 準備（この 1 冊で 1 回だけ）

パッケージを入れて、教材リポジトリを取ってくる。
**観測データは、必要になったところで各モジュールが自分で取得する**
（既にあれば何もしないので、上から流し直しても無駄が無い）。"""

PREAMBLE_PIP = "!pip install -q eispac fiasco demregpy"

PREAMBLE_BOOT = '''import os
import subprocess
import sys

REPO = "https://github.com/hottahd/EIS_practice.git"
if not os.path.exists("scripts/lines_warren2012.py"):      # リポジトリの外にいる
    if not os.path.exists("EIS_practice"):
        print("教材リポジトリを取得中 ...")
        subprocess.run(["git", "clone", "-q", REPO], check=True)
    os.chdir("EIS_practice")
sys.path.insert(0, "scripts")
print("作業ディレクトリ:", os.getcwd())'''


def _is_install(src):
    return src.lstrip().startswith("!") and "pip install" in src


def _is_bootstrap(src):
    return "REPO = " in src


def _is_legacy_clone(src):
    # モジュール 0 の `!git clone` + `%cd`。通し版では上の bootstrap が代わりをする
    return "git clone" in src or src.lstrip().startswith("%cd")

COMBINED_HEADER = """# Hinode/EIS データ解析講習会（通し版）

Warren, Winebarger & Brooks (2012), ApJ 759, 141 と同じ解析を最後まで通します。

**この 1 冊にモジュール 0〜7 が全部入っています。** 上から順に実行してください。
インストールとデータ取得は**最初の 1 回だけ**で済みます。

| | 内容 | 目安 |
|---|---|---|
| 0 | 環境構築とデータ取得 | 5 分（EIS 94 MB のダウンロード込み） |
| 1 | EIS のデータを見る | 40 分 |
| 2 | スペクトル線フィット | 50 分 |
| 3 | AIA 94 → Fe XVIII | 50 分 |
| 4 | 座標合わせと箱の選択 | 40 分 |
| 5 | **論文 Table 2 と答え合わせ** | 30 分 |
| 6 | 寄与関数と EM loci | 40 分 |
| 7 | DEM インバージョン | 60 分 |

半日コースは 5 まで、1 日コースは 7 まで。
モジュールごとに分かれた版は
[`notebooks/`](https://github.com/hottahd/EIS_practice/tree/main/notebooks) にあります。

**★ Colab の保存について**: GitHub から開いたノートは読み取り専用の一時セッションです。
編集や実行結果を残したいときは「ファイル → ドライブにコピーを保存」。
仮想マシンが切れるとダウンロードしたデータも消えますが、
その場合は上から流し直せば復帰できます（既にあるファイルは取得し直しません）。
"""


def build_combined(paths, name="EIS_workshop"):
    """モジュールごとのソースを繋いで 1 冊のノートブックにする。"""
    md = lambda t: {"cell_type": "markdown", "metadata": {}, "source": t}
    code = lambda t: {"cell_type": "code", "metadata": {}, "execution_count": None,
                      "outputs": [], "source": t}

    cells = [md(COMBINED_HEADER), md(PREAMBLE_MD), code(PREAMBLE_PIP),
             code(PREAMBLE_BOOT)]
    for p in paths:
        for c in parse(p):
            src = c["source"]
            if c["cell_type"] == "markdown":
                if src.startswith(BOOTSTRAP_MD):
                    continue                     # 先頭に 1 つ置いたので不要
                cells.append(c)
                continue
            if _is_install(src) or _is_bootstrap(src) or _is_legacy_clone(src):
                continue                         # 準備は先頭の 1 組だけ
            cells.append(c)

    out = HERE / (name + ".ipynb")
    out.write_text(json.dumps(notebook(cells), ensure_ascii=False, indent=1))
    return out, len(cells)


def main():
    if not SRC.is_dir():
        sys.exit(f"{SRC} が無い")
    paths = sorted(SRC.glob("*.py"))
    out, n = build_combined(paths)
    print(f"{len(paths)} モジュール -> {out.name}  ({n} cells)")


if __name__ == "__main__":
    main()
