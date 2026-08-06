"""ノートブックのコードセルをローカルで通しで実行して壊れていないか確かめる。

**なぜ要るか**: 講習会当日に「動かない」が一番まずい。
.ipynb を Colab で手で流して確認するのは再現性が無いので、
`notebooks/src/*.py` のコードセルをこのスクリプトで**上から順に exec** する。
セル間の変数は共有されるので、ノートを頭から実行したのと同じ状態になる。

Colab 専用のセル（`!pip`, `!wget`, `%cd` など shell/magic を含むもの）は
ローカルでは実行できないので飛ばす。飛ばしたセルは番号を表示するので、
「何を検証していないか」が分かる。

使い方（リポジトリのルートで）:
    python notebooks/verify_notebooks.py              # 全部
    python notebooks/verify_notebooks.py 01 02        # 番号で絞る

前提: eispac などが入った環境。データは data/ に取得済みであること
      （`python scripts/fetch_data.py`）。
"""
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


def is_colab_only(src):
    """shell (!) や magic (%) を含むセルは Colab 専用とみなす。"""
    for line in src.split("\n"):
        s = line.strip()
        if s.startswith("!") or s.startswith("%"):
            return True
    return False


def run(path):
    cells = [c for c in parse(path) if c["cell_type"] == "code"]
    ns = {"__name__": "__main__"}
    skipped, t0 = [], time.time()
    print(f"\n=== {path.name}  ({len(cells)} code cells) ===")
    for k, c in enumerate(cells):
        if is_colab_only(c["source"]):
            skipped.append(k)
            continue
        t = time.time()
        try:
            exec(compile(c["source"], f"{path.name}[cell {k}]", "exec"), ns)
        except Exception:
            print(f"  cell {k:2d}  ✗ FAILED")
            print("  " + "-" * 60)
            traceback.print_exc()
            print("  " + "-" * 60)
            print("  該当セルの中身:")
            for line in c["source"].split("\n"):
                print("    | " + line)
            return False
        dt = time.time() - t
        print(f"  cell {k:2d}  ok  {dt:6.1f}s" if dt > 1 else f"  cell {k:2d}  ok")
    print(f"  --> 全 {len(cells)} セル中 {len(cells)-len(skipped)} 個を実行, "
          f"{len(skipped)} 個は Colab 専用のため skip {skipped}  "
          f"({time.time()-t0:.0f}s)")
    return True


def main(args):
    import os
    os.chdir(ROOT)          # ノートは Colab でリポジトリのルートに cd している前提
    srcs = sorted((HERE / "src").glob("*.py"))
    if args:
        srcs = [p for p in srcs if any(p.name.startswith(a) for a in args)]
    if not srcs:
        sys.exit("対象のノートブックが無い")
    ng = [p.name for p in srcs if not run(p)]
    print()
    if ng:
        sys.exit(f"失敗: {', '.join(ng)}")
    print("すべて通った")


if __name__ == "__main__":
    main(sys.argv[1:])
