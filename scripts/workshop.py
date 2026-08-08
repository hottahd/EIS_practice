"""講習会ノートブックの共通ヘルパ（データ取得と定数だけ）。

★ **ここには解析を書かない。** フィットの仕方、平均の取り方、速度の出し方は
  すべてノートブック側に書く。ここに畳むと受講者から見えなくなり、
  「何をやっているか分からないまま結果だけ出る」ことになる。

ここに置いてよいのは:
  - どのノートからも同じ値を使いたい定数（データのパス、採用する箱）
  - 教える中身ではない雑務（ファイルのダウンロード）
"""
import os
import urllib.request

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EIS_FILE = "data/eis/eis_20110702_030712.data.h5"
EIS_BASE = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
AIA_DIR = "data/sdo/synoptic"
AIA_BASE = "http://jsoc.stanford.edu/data/aia/synoptic/2011/07/02/H0300"

# 論文 Table 2 の region 7 で採用した inter-moss 箱（モジュール 4 で選ぶ）
BOX = dict(y0=244, y1=274, x0=32, x1=40)

# AIA 94 A の Fe XVIII 分離（論文 Appendix。指数の誤植は修正済み。モジュール 3）
_A = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]


def download(url, path):
    """path が無ければ url から落とす。あれば何もしない。"""
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"  downloading {os.path.basename(path)} ...")
    urllib.request.urlretrieve(url, path)
    return path


def ensure_eis():
    """EIS level-1 HDF5（94 MB + 421 KB）を用意する。"""
    for ext in ("data", "head"):
        download(f"{EIS_BASE}/eis_20110702_030712.{ext}.h5",
                 f"data/eis/eis_20110702_030712.{ext}.h5")
    return EIS_FILE


