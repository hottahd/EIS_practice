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


# --- 章ごとのソースを繋いで 1 冊にする ------------------------------------
#
# ★ ノートに出る文章は**すべて notebooks/src/*.py の中**にある。
#   このスクリプトは文章を持たない（持たせると編集する人が探せなくなる）。
#
# 準備（インストール・リポジトリ取得）は先頭の 00_setup にだけ書く。
# 2 つ目以降の章に同じセルがあれば落とす。


def _is_install(src):
    return src.lstrip().startswith("!") and "pip install" in src


def _is_bootstrap(src):
    return "REPO = " in src


def build_combined(paths, name="EIS_workshop"):
    """章ごとのソースを順に繋ぐ。準備のセルは最初の 1 組だけ残す。"""
    cells = []
    seen_install = seen_boot = False
    for p in paths:
        for c in parse(p):
            src = c["source"]
            if c["cell_type"] == "code":
                if _is_install(src):
                    if seen_install:
                        continue
                    seen_install = True
                elif _is_bootstrap(src):
                    if seen_boot:
                        continue
                    seen_boot = True
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
