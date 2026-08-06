"""ノートブックのコードセルをローカルで通しで実行して壊れていないか確かめる。

**なぜ要るか**: 講習会当日に「動かない」が一番まずい。
Colab で手で流して確認するのは再現性が無いので、コードセルを
**上から順に exec** する。セル間で変数を共有するので、
ノートを頭から実行したのと同じ状態になる。

Colab 専用のセル（`!pip`, `%cd` など shell/magic を含むもの）はローカルでは
実行できないので飛ばす。飛ばしたセルは番号を表示するので、
「何を検証していないか」が分かる。

使い方（リポジトリのルートで）:

    python notebooks/verify_notebooks.py            # 配布物 EIS_workshop.ipynb を通しで
    python notebooks/verify_notebooks.py 03 04      # 編集中のモジュールだけ（src/*.py）

★ 通し版は**全モジュールが 1 つの名前空間**で動く。モジュール間で変数名が
  衝突していないかは、こちらを実行しないと分からない。編集したら必ず通すこと。

前提: eispac などが入った環境。データは data/ に取得済みであること
      （`python scripts/fetch_data.py`。無ければノートが自分で取りに行く）。
"""
import json
import pathlib
import sys
import time
import traceback

import matplotlib
matplotlib.use("Agg")

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))
from build_notebooks import parse  # noqa: E402

ROOT = HERE.parent
COMBINED = HERE / "EIS_workshop.ipynb"


def is_colab_only(src):
    """shell (!) や magic (%) を含むセルは Colab 専用とみなす。"""
    return any(l.strip().startswith(("!", "%")) for l in src.split("\n"))


def run(label, sources):
    """コードセルの列を 1 つの名前空間で順に実行する。"""
    ns = {"__name__": "__main__"}
    skipped, t0 = [], time.time()
    print(f"\n=== {label}  ({len(sources)} code cells) ===")
    for k, src in enumerate(sources):
        if is_colab_only(src):
            skipped.append(k)
            continue
        t = time.time()
        try:
            exec(compile(src, f"{label}[cell {k}]", "exec"), ns)
        except Exception:
            print(f"  cell {k:2d}  ✗ FAILED")
            print("  " + "-" * 60)
            traceback.print_exc()
            print("  " + "-" * 60)
            print("  該当セルの中身:")
            for line in src.split("\n"):
                print("    | " + line)
            return False
        dt = time.time() - t
        print(f"  cell {k:2d}  ok  {dt:6.1f}s" if dt > 1 else f"  cell {k:2d}  ok")
    print(f"  --> 全 {len(sources)} セル中 {len(sources)-len(skipped)} 個を実行, "
          f"{len(skipped)} 個は Colab 専用のため skip {skipped}  "
          f"({time.time()-t0:.0f}s)")
    return True


def main(args):
    import os
    os.chdir(ROOT)          # ノートは Colab でリポジトリのルートに cd している前提

    if args:                # モジュールを番号で絞る（編集中の確認用）
        srcs = [p for p in sorted((HERE / "src").glob("*.py"))
                if any(p.name.startswith(a) for a in args)]
        if not srcs:
            sys.exit("対象のモジュールが無い")
        ng = [p.name for p in srcs
              if not run(p.name, [c["source"] for c in parse(p)
                                  if c["cell_type"] == "code"])]
    else:                   # 配布物そのものを通しで
        if not COMBINED.exists():
            sys.exit(f"{COMBINED} が無い。先に build_notebooks.py を実行すること")
        nb = json.loads(COMBINED.read_text())
        ng = [] if run(COMBINED.name,
                       [c["source"] for c in nb["cells"]
                        if c["cell_type"] == "code"]) else [COMBINED.name]

    print()
    if ng:
        sys.exit(f"失敗: {', '.join(ng)}")
    print("すべて通った")


if __name__ == "__main__":
    main(sys.argv[1:])
