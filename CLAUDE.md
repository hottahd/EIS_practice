# EIS_practice — Hinode/EIS データ解析講習会の準備

## このリポジトリの目的

日本の太陽物理研究者が Hinode/EIS のデータを自力で扱えるようになるための
**講習会（ハンズオン）教材**を作る。

到達目標は明確で、**Warren, Winebarger & Brooks (2012), ApJ 759, 141**
"A Systematic Survey of High-Temperature Emission in Solar Active Regions"
と同じ解析を受講者が自分で実行できるようになること。

論文本体: `papers/Warren_2012_*.pdf`

## セッションをまたぐための記録ルール（重要）

Claude の会話は切れる前提で作業する。**作業内容は必ず `docs/` に記録すること。**

| ファイル | 役割 |
|---|---|
| `CLAUDE.md` | このファイル。プロジェクト全体像と現在地 |
| `docs/00_log.md` | 作業ログ（日付順に追記。何をやったか／次に何をやるか） |
| `docs/01_paper_analysis.md` | Warren+2012 の手法を分解したもの。教材設計の元ネタ |
| `docs/02_workshop_plan.md` | 講習会のカリキュラム設計 |
| `docs/03_environment.md` | 環境構築手順（受講者に配る想定） |

新しいセッションを始めたら、まず `docs/00_log.md` の末尾を読むこと。

## 検証済みスクリプト（`scripts/`）

論文の解析を通しで再現できることを実データで確認済み。実行順:

| スクリプト | 役割 |
|---|---|
| `lines_warren2012.py` | 論文 Table 2 の 22 輝線 + eispac テンプレート対応表。成分の自動判定 `pick_component()` |
| `quicklook_raster.py` | フィット無しのラスタークイックルック（数秒）。箱選び用 |
| `fit_box_spectra.py` | 箱内平均 → 22 輝線フィット（8 秒） |
| `download_sdo.py` | AIA 94/171/193 + HMI を VSO から取得（登録不要・リトライ付き） |
| `aia_fe18.py` | AIA 94 → Fe XVIII 分離（論文 Appendix。**誤植を修正済み**） |
| `make_fe18_map.py` | 論文 Figure 1–3 相当の 5 枚並び図 |
| `coalign_eis_aia.py` | EIS↔AIA 相互相関で座標合わせ、AIA を EIS 格子へ |
| `select_intermoss.py` | inter-moss 箱の候補を機械的に洗い出す |
| `compare_table2.py` | 論文 Table 2 と数値照合。「箱の選び方」を傾きで診断 |

## 到達済みの成果（2026-08-06）

**論文 Table 2 を再現できた。** 箱 y=[244:274], x=[32:40] で
21 輝線中 14 本が論文の 15% 以内、median ratio 0.89、ばらつき 0.10 dex。

準備の過程で 2 つの重要な発見:
1. **論文 Eq.(A1) の指数は誤植**。正しくは定数項つき 3 次式（実データで確定）
2. **eispac に Ca XVII 192.858 のブレンド分離テンプレートが無い**（未解決・要自作）

詳細は `docs/00_log.md`。

## 別のマシンにクローンしたときの再開手順

観測データはリポジトリに入っていない（`.gitignore` で `data/` を除外）。
下記 3 ステップで完全に再現できる。

```bash
git clone git@github.com:hottahd/EIS_practice.git
cd EIS_practice

# 1. 環境（詳細は docs/03_environment.md）
mamba create -n eis -c conda-forge -y python=3.12 numpy scipy matplotlib \
      astropy sunpy ndcube h5py pandas jupyterlab tqdm zeep drms parfive
mamba activate eis
pip install eispac demregpy aiapy fiasco

# 2. データ取得（EIS 94 MB + SDO 46 MB、数分）
python scripts/fetch_data.py            # 既定 = region 7 (論文 Table 2 の天体)

# 3. 動作確認：ここまでが再現済みの到達点
python scripts/quicklook_raster.py data/eis/eis_20110702_030712.data.h5 figures/q.png
python scripts/make_fe18_map.py       data/eis/eis_20110702_030712.data.h5
python scripts/coalign_eis_aia.py     data/eis/eis_20110702_030712.data.h5
python scripts/compare_table2.py      data/eis/eis_20110702_030712.data.h5 "244:274,32:40"
#   -> 21 輝線中 14 本が論文 Table 2 の 15% 以内、median ratio 0.89 になれば成功
```

`.claude/settings.json` の承認ルールにはこのマシン固有の
miniforge 絶対パスが入っているが、`Bash(python *)` 等の一般パターンも
入れてあるので別マシンでもそのまま使える。

## ★ IDL / SSWIDL があるサーバーでやるべきこと

Python (eispac) だけでは解決できず、SSW が要る項目。
**目的は「IDL に移行すること」ではなく、小さな検証結果を持ち帰ること。**
講習会の本線は Python/Colab のまま変えない。

優先度順:

1. **`eis_auto_fit` と eispac のフィット結果を突き合わせる** ← 最優先
   未解決の食い違い（Fe XVI 262.984 が 0.62、S XIII 256.686 が 0.71、
   Si VII 275.368 が 0.39）が eispac 固有なのか実在するのかを決める。
   同じ箱 y=[244:274], x=[32:40] で同じ 22 輝線を出して比較する。
   → 持ち帰るもの: 輝線ごとの強度の表（小さな CSV）

2. **打ち上げ後の EIS 較正カーブを作る**
   eispac は**打ち上げ前較正しか持っていない**（`radcal/<win>_pre` のみ。
   Del Zanna 2013 / Warren et al. 2014 の実装は無い）。
   SSW の `eis_recalibrate_intensity.pro` で較正曲線を出し、
   波長ごとの配列として保存する。
   → 持ち帰るもの: 較正カーブの .npz（eispac の `read_cube(radcal=...)` に渡せる）
   ※ ただし今回の食い違いは SW/LW で系統差が無いので、
      較正が主犯である可能性は低い（docs/00_log.md の訂正を参照）。
      それでも講習会モジュール 9 の教材として価値がある。

3. **Ca XVII 192.858 のブレンド分離を SSW でやる**
   eispac にはこのブレンドを解くテンプレートが無い（論文値の 4.75 倍が出る）。
   Ko et al. (2009) の方法を SSW で実行し、正解を作る。
   → 持ち帰るもの: 正しい Ca XVII 強度。自作 Python テンプレートの検証に使う

4. **PINTofALE の MCMC_DEM（入手できれば）**
   PINTofALE は SSW には含まれない。v2.97 (2016) で、公式サイトの
   ダウンロードは廃止され、著者に連絡して MEGA の共有リンクをもらう必要がある。
   入手できたら論文と同じ手法の DEM を 1 つ作る。
   → 持ち帰るもの: 参照 DEM。Python 側 (demregpy) の結果を較正できる
   ※ IDL があれば GDL の互換性問題は全部回避できる

5. CHIANTI IDL で寄与関数を出し、fiasco の結果と突き合わせる（余力があれば）

## 環境（2026-08-06 時点で確認済み）

- Python: `/home/sc/c0234hotta/miniforge3/bin/python` (3.12.13, miniforge base)
  - 導入済み: numpy 2.5.1, scipy 1.18.0, astropy 8.0.1, sunpy 8.0.0, matplotlib 3.11.1, h5py 3.16.0
  - **未導入**: eispac, ndcube, aiapy, demregpy, ChiantiPy → 講習会用の専用 conda env を作る想定
- IDL / SolarSoft: **この計算機には無い**（`idl`, `gdl` ともに不在, `$SSW` 未設定）
- ネットワーク: 外部アクセス可
  - EIS Level-1 HDF5 アーカイブ (NRL) : https://eis.nrl.navy.mil/level1/hdf5/YYYY/MM/DD/ → 200 OK
  - JSOC (AIA/HMI) : http://jsoc.stanford.edu → 200 OK
  - PyPI : 到達可能（eispac 0.99.4, demregpy 1.0.0, aiapy 0.12.1, ChiantiPy 0.16.0 が入手可）
- ディスク: `/scr/a000/c0234hotta` に 102 TB 空き（データ置き場に十分）

## 再現性の確認（済）

論文 Table 2（DEM モデルの観測強度・計算強度の比較表）は Table 1 の region 7
＝ 2011-07-02 03:07:12, NOAA 1243 のデータ。対応する Level-1 HDF5 が
アーカイブに実在することを確認した:

```
https://eis.nrl.navy.mil/level1/hdf5/2011/07/02/eis_20110702_030712.data.h5
https://eis.nrl.navy.mil/level1/hdf5/2011/07/02/eis_20110702_030712.head.h5
```

→ **この 1 天体を「答え合わせのできる題材」として講習会の主教材にするのが自然。**

同様に region 8（2010-07-23 14:32:10, NOAA 1089 = Warren et al. 2011 と同じ AR）も
`eis_20100723_143210.{data,head}.h5` が存在することを確認済み。
