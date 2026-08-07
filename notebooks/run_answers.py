"""答えノートを実行して、出力を埋め込んだ状態で保存する。

**なぜ要るか**: 答えは本編と別のノートに置いてある。Colab はノート 1 冊ごとに
VM が変わるので、答えを見るためだけにデータを取り直すのは無駄。
**実行結果を埋め込んで配れば、開くだけで答えと結果が読める。**

    python notebooks/build_notebooks.py   # まずソースからノートを作る
    python notebooks/run_answers.py       # 実行して出力を埋め込む（数分）

`!pip` などの Colab 専用セルは飛ばす。
"""
import json
import os
import pathlib
import sys
import tempfile

import nbformat
from nbclient import NotebookClient

HERE = pathlib.Path(__file__).parent
ROOT = HERE.parent
NB = HERE / "EIS_workshop_answers.ipynb"


# フィッターやファイル読み込みの進捗表示は、配る答えノートには要らない
NOISE = ("+ working on exposure", "+ computing fits", "+ running mpfit",
         "作業ディレクトリ:",
         "Found a wavelength", "Data file,", "Header file,", "Finished computing fits",
         "runtime :", "spectra fit without issues", "spectra have <",
         "spectra have bad or invalid")


def clean_outputs(nb):
    """進捗ログを落として、答えノートの出力を読みやすくする。"""
    for c in nb.cells:
        if c.cell_type != "code" or not c.get("outputs"):
            continue
        keep = []
        for o in c.outputs:
            if o.get("output_type") != "stream":
                keep.append(o)
                continue
            text = o["text"] if isinstance(o["text"], str) else "".join(o["text"])
            lines = [l for l in text.split("\n")
                     if l.strip() and not any(n in l for n in NOISE)
                     and not l.strip().startswith("/")]
            if lines:
                o["text"] = "\n".join(lines) + "\n"
                keep.append(o)
        c.outputs = keep


def _make_kernelspec():
    """いま動いている python を指すカーネルを一時的に用意する。

    登録済みのカーネルは環境構築時の絶対パスを持っており、
    マウントが変わると `FileNotFoundError` になる（実際になった）。
    """
    tmp = pathlib.Path(tempfile.mkdtemp())
    d = tmp / "kernels" / "eis-run"
    d.mkdir(parents=True)
    (d / "kernel.json").write_text(json.dumps({
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "eis-run",
        "language": "python",
    }))
    os.environ["JUPYTER_PATH"] = str(tmp)
    return "eis-run"


def main():
    if not NB.exists():
        sys.exit(f"{NB} が無い。先に build_notebooks.py を実行すること")

    nb = nbformat.read(NB, as_version=4)

    # Colab 専用セル（!pip など）は中身を差し替えて実行させない
    skipped, saved = [], {}
    for k, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        if any(l.strip().startswith(("!", "%")) for l in c.source.split("\n")):
            saved[k] = c.source
            c.source = "pass"
            skipped.append(k)

    kernel = _make_kernelspec()
    print(f"実行中（{len(nb.cells)} セル、Colab 専用の {len(skipped)} 個は skip）...")
    NotebookClient(nb, timeout=1800, kernel_name=kernel,
                   resources={"metadata": {"path": str(ROOT)}}).execute()

    for k, src in saved.items():
        nb.cells[k].source = src
        nb.cells[k].outputs = []
    clean_outputs(nb)
    nbformat.write(nb, NB)

    n_out = sum(1 for c in nb.cells if c.cell_type == "code" and c.get("outputs"))
    print(f"出力を埋め込んで保存: {NB.name}  "
          f"（{n_out} セルに出力あり、{NB.stat().st_size/1e6:.1f} MB）")


if __name__ == "__main__":
    main()
