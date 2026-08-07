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
RESTART_MD_HEAD = "### インストール直後のランタイム再起動"

PREAMBLE_MD = """<!-- 準備は 00_setup 側に書いてある -->"""

PREAMBLE_PIP = "!pip install -q eispac fiasco demregpy"

PREAMBLE_RESTART_MD = """### インストール直後のランタイム再起動について
#
Colab では、pip が `numpy` などを入れ替えると、**実行中のセッションが
古いモジュールを掴んだまま**になり、あとで次のようなエラーが出ることがある:
#
```
ImportError: cannot import name '_center' from 'numpy._core.umath'
```
#
これはインストールの失敗ではなく、**再起動すれば直る**。
次のセルが入れ替えを検出して、必要なときだけ自動で再起動する。
#
**再起動が起きたら、もう一度このノートを先頭から実行すること。**
2 回目はインストールもダウンロードも済んでいるので一瞬で終わる。"""

PREAMBLE_RESTART = '''import sys
from importlib.metadata import version

need_restart = False
try:
    loaded = sys.modules["numpy"].__version__ if "numpy" in sys.modules else None
    if loaded is not None and loaded != version("numpy"):
        need_restart = True
        print(f"numpy が {loaded} -> {version('numpy')} に入れ替わりました")
except Exception as e:                      # 判定自体が失敗したら念のため再起動
    need_restart = True
    print("numpy の状態を確認できませんでした:", e)

if need_restart:
    print("ランタイムを再起動します。"
          "再起動したら、もう一度このノートを先頭から実行してください。")
    try:
        import IPython
        ipy = IPython.get_ipython()
        if ipy is not None:
            ipy.kernel.do_shutdown(True)    # Colab のランタイム再起動
    except Exception:
        import os
        os.kill(os.getpid(), 9)
else:
    print("numpy の入れ替えは起きていません。このまま先へ進んで大丈夫です。")'''

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


def _is_restart(src):
    return "need_restart" in src


def _is_legacy_clone(src):
    # モジュール 0 の `!git clone` + `%cd`。通し版では上の bootstrap が代わりをする
    return "git clone" in src or src.lstrip().startswith("%cd")

COMBINED_HEADER = """# Hinode/EIS データ解析講習会

**Solar-C (EUVST) の準備**として、いま手に入る Hinode/EIS のデータで
分光解析をひととおり通します。

EUVST の EUV バンド 170–215 Å は EIS の短波長帯とほぼ同じで、**同じ輝線を撮ります**。
**今日やることは、そのまま 2028 年に使えます。**

| 章 | 内容 | 目安 |
|---|---|---|
| 1 | EIS のデータを見る | 20 分 |
| 2 | フィットして**強度**を出す | 35 分 |
| 3 | **速度**を出す | 50 分 |
| 4 | 線幅から**非熱的速度**を出す | 30 分 |
| 5 | **温度分布 (DEM)** を出す | 30 分 |
| 付録 | 自分の研究で使うときに読む（当日は走らせません） | — |

上から順に実行してください。インストールとデータ取得は**最初の 1 回だけ**です。
"""


def build_combined(paths, name="EIS_workshop"):
    """モジュールごとのソースを繋いで 1 冊のノートブックにする。"""
    md = lambda t: {"cell_type": "markdown", "metadata": {}, "source": t}
    code = lambda t: {"cell_type": "code", "metadata": {}, "execution_count": None,
                      "outputs": [], "source": t}

    cells = [md(COMBINED_HEADER), md(PREAMBLE_MD), code(PREAMBLE_PIP),
             md(PREAMBLE_RESTART_MD), code(PREAMBLE_RESTART), code(PREAMBLE_BOOT)]
    for p in paths:
        for c in parse(p):
            src = c["source"]
            if c["cell_type"] == "markdown":
                if src.startswith(BOOTSTRAP_MD) or src.startswith(RESTART_MD_HEAD):
                    continue                     # 先頭に 1 つ置いたので不要
                cells.append(c)
                continue
            if (_is_install(src) or _is_bootstrap(src) or _is_restart(src)
                    or _is_legacy_clone(src)):
                continue                         # 準備は先頭の 1 組だけ
            cells.append(c)

    out = HERE / (name + ".ipynb")

    # ★ .ipynb を直接編集していた場合の保険。
    #   生成物がソースより新しければ、上書きする前に控えを取って警告する
    #   （Colab や Jupyter で直接いじった内容を黙って消さないため）。
    if out.exists():
        newest_src = max(p.stat().st_mtime for p in paths)
        if out.stat().st_mtime > newest_src + 1:
            bak = out.with_suffix(".ipynb.bak")
            bak.write_text(out.read_text())
            print(f"  ! {out.name} がソースより新しい（直接編集された可能性）。"
                  f"控えを {bak.name} に取って上書きします")

    out.write_text(json.dumps(notebook(cells), ensure_ascii=False, indent=1))
    return out, len(cells)


ANSWERS = HERE / "answers"


def build_answers():
    """演習の答えを別のノートにする（本編には答えを置かない）。

    ★ 答えノートは **実行結果を埋め込んだ状態**で配る。
      `python notebooks/run_answers.py` で実行して出力を入れる。
      Colab で開くだけで答えと結果が読めるので、VM を立ち上げ直さずに済む。
    """
    if not ANSWERS.is_dir():
        return None, 0
    cells = []
    for p in sorted(ANSWERS.glob("*.py")):
        cells += parse(p)
    out = HERE / "EIS_workshop_answers.ipynb"
    # 既に実行結果が入っていれば、それを残したまま中身だけ差し替える
    old = {}
    if out.exists():
        try:
            prev = json.loads(out.read_text())
            old = {c["source"]: c for c in prev["cells"] if c["cell_type"] == "code"}
        except Exception:
            pass
    for c in cells:
        if c["cell_type"] == "code" and c["source"] in old:
            c["outputs"] = old[c["source"]].get("outputs", [])
            c["execution_count"] = old[c["source"]].get("execution_count")
    out.write_text(json.dumps(notebook(cells), ensure_ascii=False, indent=1))
    return out, len(cells)


def main():
    if not SRC.is_dir():
        sys.exit(f"{SRC} が無い")
    paths = sorted(SRC.glob("*.py"))
    out, n = build_combined(paths)
    print(f"{len(paths)} 章 -> {out.name}  ({n} cells)")
    out_a, n_a = build_answers()
    if out_a:
        print(f"答え      -> {out_a.name}  ({n_a} cells)")


if __name__ == "__main__":
    main()
