"""講習会用の Colab ノートブックを生成する。

.ipynb を直接手で書くと JSON が読みにくく、レビューも差分も辛い。
ここでは **セルの中身を素の Python / Markdown として書き**、
このスクリプトで .ipynb に変換する。

  - ノートブックの定義は notebooks/src/NN_*.py に置く（後述の書式）
  - `python notebooks/build_notebooks.py` で notebooks/*.ipynb を生成
  - コードセルの中身は普通の Python なので、単体で実行して検証できる

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


def build(path):
    cells = parse(path)
    nb = {
        "cells": cells,
        "metadata": {
            "colab": {"provenance": [], "toc_visible": True},
            "kernelspec": {"display_name": "Python 3", "name": "python3"},
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 0,
    }
    out = HERE / (path.stem + ".ipynb")
    out.write_text(json.dumps(nb, ensure_ascii=False, indent=1))
    return out, len(cells)


def main():
    if not SRC.is_dir():
        sys.exit(f"{SRC} が無い")
    for p in sorted(SRC.glob("*.py")):
        out, n = build(p)
        print(f"{p.name:<34} -> {out.name:<34} ({n} cells)")


if __name__ == "__main__":
    main()
