"""別のマシンにクローンしたあと、解析に必要なデータを一括で取ってくる。

観測データは巨大なのでリポジトリには入れていない（.gitignore）。
これを 1 回走らせれば `data/` が再構成され、他のスクリプトがそのまま動く。

    python scripts/fetch_data.py            # 既定 = region 7 (2011-07-02)
    python scripts/fetch_data.py --region 8 # 別の天体

必要なもの: eispac, sunpy, aiapy（`docs/03_environment.md` 参照）
所要: EIS 94 MB + SDO 46 MB、回線次第で数分
"""
import os
import sys
import argparse
import urllib.request

# docs/01_paper_analysis.md の Table 1 より（論文が使った 15 天体）
# region 番号 -> (EIS ファイル名, ラスター中点時刻 UT)
REGIONS = {
    1:  ("eis_20100619_014433", "2010-06-19T01:57:44"),
    2:  ("eis_20110212_143019", "2011-02-12T15:32:13"),
    3:  ("eis_20100621_011541", "2010-06-21T01:46:37"),
    4:  ("eis_20110725_090513", "2011-07-25T09:36:09"),
    5:  ("eis_20110131_102326", "2011-01-31T11:25:19"),
    6:  ("eis_20110121_133954", "2011-01-21T14:10:50"),
    7:  ("eis_20110702_030712", "2011-07-02T03:38:08"),   # 論文 Table 2 の天体
    8:  ("eis_20100723_143210", "2010-07-23T15:03:07"),
    9:  ("eis_20100929_223226", "2010-09-29T23:51:36"),
    10: ("eis_20110419_123027", "2011-04-19T13:32:20"),
    11: ("eis_20110411_105848", "2011-04-11T12:00:42"),
    12: ("eis_20110821_105251", "2011-08-21T12:25:42"),
    13: ("eis_20110415_001526", "2011-04-15T01:17:19"),
    14: ("eis_20111108_181234", "2011-11-08T19:14:27"),
    15: ("eis_20111110_100028", "2011-11-10T11:33:19"),
}

EIS_BASE = "https://eis.nrl.navy.mil/level1/hdf5"


def fetch_eis(tag, outdir="data/eis"):
    """EIS Level-1 HDF5 を NRL のアーカイブから取る。"""
    os.makedirs(outdir, exist_ok=True)
    ymd = tag.split("_")[1]                      # eis_20110702_030712 -> 20110702
    y, m, d = ymd[:4], ymd[4:6], ymd[6:8]
    got = []
    for suffix in ("data.h5", "head.h5"):
        fname = f"{tag}.{suffix}"
        dest = os.path.join(outdir, fname)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"  既にある: {fname}")
            got.append(dest)
            continue
        url = f"{EIS_BASE}/{y}/{m}/{d}/{fname}"
        print(f"  取得中: {url}")
        urllib.request.urlretrieve(url, dest)
        print(f"    -> {dest} ({os.path.getsize(dest)/1e6:.1f} MB)")
        got.append(dest)
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", type=int, default=7,
                    help="論文 Table 1 の region 番号 (既定: 7)")
    ap.add_argument("--eis-only", action="store_true", help="SDO は取らない")
    args = ap.parse_args()

    if args.region not in REGIONS:
        sys.exit(f"region は 1-15 のいずれか。指定: {args.region}")
    tag, midtime = REGIONS[args.region]

    print(f"=== region {args.region}: {tag}  (ラスター中点 {midtime}) ===")
    print("[1/2] EIS Level-1 HDF5")
    fetch_eis(tag)

    if args.eis_only:
        print("SDO はスキップ")
        return

    print("[2/2] SDO/AIA 94,171,193 + HMI （VSO 経由・登録不要）")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from download_sdo import download
    download(midtime, "data/sdo")

    print("\n完了。次はこれが動くはず:")
    print(f"  python scripts/quicklook_raster.py data/eis/{tag}.data.h5 "
          f"figures/quicklook.png")


if __name__ == "__main__":
    main()
