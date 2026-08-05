"""EIS ラスターの中点時刻に合わせて SDO/AIA と SDO/HMI のデータを取得する。

Warren et al. (2012) は JSOC からラスター中点 +-300 s のデータを取っている。
ここでは中点に最も近い 1 枚ずつを取る。

取得経路について（講習会での注意点）:
  * VSO (Virtual Solar Observatory) 経由 → **ユーザ登録不要**。
    ただしフルディスク画像（1枚 65 MB）になる。
  * JSOC 経由 (a.jsoc.Series + a.jsoc.Notify) → 切り出し(cutout)が使えて軽いが、
    JSOC に登録したメールアドレスが必要。
  講習会では登録の手間が事故のもとなので VSO を既定にする。
  Colab でも問題なく落とせる（合計 ~260 MB、数分）。

使い方:
    python download_sdo.py <YYYY-MM-DDTHH:MM:SS> [出力ディレクトリ]
"""
import sys
import os
import astropy.units as u
from astropy.time import Time
from sunpy.net import Fido, attrs as a

AIA_WAVES = [94, 171, 193] * u.AA   # Fe XVIII 分離に必要な 3 波長


def _fetch_one(query_row, outdir, label, max_try=3):
    """VSO は時々 staging に失敗するのでリトライする（講習会中の事故を減らす）。"""
    for attempt in range(1, max_try + 1):
        got = Fido.fetch(query_row, path=os.path.join(outdir, "{file}"), progress=False)
        if len(got) > 0 and len(got.errors) == 0:
            print(f"{label}: {got[0]}")
            return list(got)
        print(f"{label}: 取得失敗 (試行 {attempt}/{max_try}) errors={got.errors}")
    print(f"!! {label} : {max_try} 回試して取得できませんでした")
    return []


def download(midtime, outdir="data/sdo", window=15 * u.s):
    t0 = Time(midtime)
    os.makedirs(outdir, exist_ok=True)
    files = []

    for w in AIA_WAVES:
        tr = a.Time(t0 - window, t0 + window)
        q = Fido.search(tr, a.Instrument.aia, a.Wavelength(w))
        if len(q[0]) == 0:
            print(f"!! AIA {w} : 見つかりません")
            continue
        # 中点に最も近い 1 枚だけ
        i = _nearest(q[0]["Start Time"], t0)
        files += _fetch_one(q[0][i:i + 1], outdir, f"AIA {w.value:.0f}")

    # HMI 磁場（line-of-sight magnetogram）は 45 s ケイデンス
    tr = a.Time(t0 - 60 * u.s, t0 + 60 * u.s)
    q = Fido.search(tr, a.Instrument.hmi, a.Physobs.los_magnetic_field)
    if len(q[0]) == 0:
        print("!! HMI : 見つかりません")
    else:
        i = _nearest(q[0]["Start Time"], t0)
        files += _fetch_one(q[0][i:i + 1], outdir, "HMI magnetogram")

    return files


def _nearest(times, t0):
    import numpy as np
    return int(np.argmin(abs((Time(times) - t0).sec)))


if __name__ == "__main__":
    mid = sys.argv[1] if len(sys.argv) > 1 else "2011-07-02T03:38:08"
    out = sys.argv[2] if len(sys.argv) > 2 else "data/sdo"
    download(mid, out)
