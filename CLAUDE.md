# EIS_practice — Hinode/EIS データ解析講習会の準備

## 公開状態とライセンス

**public**（2026-08-07 に公開）。**CC BY 4.0**（`LICENSE`）。

- **論文 PDF はリポジトリに入れない。** 再配布の許諾が無い（出版社版は © AAS/IOP、
  arXiv 版は nonexclusive-distrib 1.0）。`.gitignore` で除外済み。
  出典と取得方法は `papers/README.md` / `papers/fetch_papers.sh`。
  **PDF を追加するコミットを作らないこと。**
- 公開前に履歴から PDF を除去済み（`git filter-repo`。2026-08-07。全 hash が変わった）。

## このリポジトリの目的

日本の太陽物理研究者が Hinode/EIS のデータを自力で扱えるようになるための
**講習会（ハンズオン）教材**を作る。

到達目標は明確で、**Warren, Winebarger & Brooks (2012), ApJ 759, 141**
"A Systematic Survey of High-Temperature Emission in Solar Active Regions"
と同じ解析を受講者が自分で実行できるようになること。

論文本体: `papers/Warren_2012_*.pdf`
**PDF はリポジトリに入っていない**（再配布の許諾が無いため）。
`bash papers/fetch_papers.sh` で arXiv 版が手元に揃う。詳細は `papers/README.md`。

## セッションをまたぐための記録ルール（重要）

Claude の会話は切れる前提で作業する。**作業内容は必ず `docs/` に記録すること。**

| ファイル | 役割 |
|---|---|
| `CLAUDE.md` | このファイル。プロジェクト全体像と現在地 |
| `docs/00_log.md` | 作業ログ（日付順に追記。何をやったか／次に何をやるか） |
| `docs/01_paper_analysis.md` | Warren+2012 の手法を分解したもの。教材設計の元ネタ |
| `docs/02_workshop_plan.md` | 講習会のカリキュラム設計 |
| `docs/03_environment.md` | 環境構築手順（受講者に配る想定） |
| `docs/04_physics_primer.md` | **理論屋向けの物理解説**（なぜそうするのかを全ステップで） |
| `docs/05_teaching_design.md` | **教材の設計方針**（作り直しの提案。★まずこれを読む） |

新しいセッションを始めたら、まず `docs/00_log.md` の末尾を読むこと。

## ★ この講習会の位置づけ（2026-08-07 に判明）

**Solar-C (EUVST) に向けた実力づくり**が目的。打ち上げ時にすぐ解析できる人を用意しておく。
EIS は練習台であって目的ではない。

- **EUVST の EUV バンド 170–215 Å は EIS の短波長帯 171–212 Å とほぼ同じ**＝同じ輝線を撮る。
  この教材で扱う線（Fe XI–XVI, Ca XIV–XVII, Si X, S X, Ar XIV）はそのまま主力になる。
- 変わるのは 0.4″・1 秒・実効面積 10–30 倍・温度被覆 2×10⁴–1.5×10⁷ K シームレス。
- 到達目標: **A** フィッティングできる／**B** 物理量を出せる
  （強度・**ドップラー速度**・**非熱的幅**）／**C** DEM 解析ができる。
- **★ 現行ノートは B がほぼ空**（速度はマップ 1 枚、線幅は「使わない」）。作り直しが要る。

→ 設計方針は `docs/05_teaching_design.md`。**ノートを触る前にこれを読むこと。**

## 講習会ノートブック（`notebooks/`）

**配るのは `notebooks/EIS_workshop.ipynb` の 1 冊だけ。** モジュール 0〜7 が全部入っている。
Colab は**ノート 1 冊ごとに VM が変わる**ので、分冊にすると pip install と
EIS 94 MB の取得を冊数ぶん繰り返すことになる。半日/1 日コースは
「どこで止めるか」の違いでしかないので、ノートを分けない。

**`.ipynb` を手で書かない。** ソースは `notebooks/src/NN_*.py` に
モジュールごとに `# %%` 区切りの素の Python で書き、
`python notebooks/build_notebooks.py` で 1 冊に繋いで生成する
（各モジュール先頭の install / bootstrap セルは先頭の 1 組にまとめられる）。

| 章 | 内容 | 状態 |
|---|---|---|
| 1 | EIS のデータを見る | ✅ |
| 2 | フィットして強度を出す | ✅ |
| 3 | **速度**を出す（wave_corr = 53 km/s、ゼロ点は自分で決める） | ✅ |
| 4 | 線幅から**非熱的速度**（装置幅は位置依存、中央値 18 km/s） | ✅ |
| 5 | 温度分布 DEM（EM loci → 1 回解く。ピーク 4.0 MK） | ✅ |
| 付録 | A 欠損値 / B 誤差 / C DEM の信頼度 / D 1/(4π) / E 較正 /<br>F AIA Fe XVIII / G Ca XVII / H Table 2 全面照合 / I 他領域 | ✅ |

**演習の答えは別ノート** `notebooks/EIS_workshop_answers.ipynb`（ソースは `notebooks/answers/`）。
本編に置くと見えてしまうため。**実行結果を埋め込んでコミットする**ので、
受講者は開くだけで読める（`python notebooks/run_answers.py`）。

**採否の基準は「楽しいか」。** 重要だが楽しくないものは**付録**に回し、
正しさは**コード側で担保**して本編にはポインタだけ置く。
「私が苦労したから」は入れない。**壊してから直す演出もしない**（`docs/05_teaching_design.md`）。

受講者は **GitHub の URL を Colab に食わせるだけ**で開ける（要 public）:
`https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop.ipynb`
→ 設計上の約束は `notebooks/README.md`。

**★ Colab は VM が使い捨て。** ノートは先頭のブートストラップで clone/cd し、
データは必要になったところで取得し（既にあれば何もしない）、
モジュール間の受け渡しファイルが無ければ `scripts/workshop.py` が作り直す。
**途中から始めても、セッションが切れても動く**ように保つこと。

**編集したら必ず検証すること**（講習会当日に動かないのが最悪なので）:

```bash
python notebooks/verify_notebooks.py        # 配布物 EIS_workshop.ipynb を通しで（約 100 秒）
python notebooks/verify_notebooks.py 03 04  # 編集中のモジュールだけ (src/*.py)
```

Colab 専用セル（`!pip` などを含むもの）は自動で飛ばし、番号を表示する。
**★ 引数なしを必ず 1 回通すこと。** 全モジュールが 1 つの名前空間で動くので、
変数名の衝突はこれでしか見つからない。

**図のラベルは英語で書く。** Colab に日本語フォントが無いので
日本語ラベルは豆腐（□）になる。**説明は日本語、図は英語**。

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
| `scan_boxes.py` | 箱を総当たりして論文に最も合う場所を探す（集計スコア） |
| `scan_perline.py` | 同上、ただし**輝線ごとの ratio を全記録**（上位だけ見ると誤る） |
| `scan_ratios.py` | 窓積分の比だけで高速に候補箱を探す |
| `fit_perpixel_box.py` | 「平均→fit」と「fit→平均」を比較（論文 §3 は前者） |
| `extract_paper_boxes.py` | **論文 Figure 1-3 の緑枠から 15 領域の箱の座標を実測** |
| `fetch_pintofale.sh` | PINTofALE 一式を取得（公式は「廃止」と書いているが生きている） |

### IDL / SolarSoft 側（`scripts/idl/`）

| スクリプト | 役割 |
|---|---|
| `run_sswidl.sh` | SSWIDL をヘッドレスで回すラッパ。踏んだ罠は全部コメントにある |
| `01_eis_prep.pro` | DARTS の level-0 → `eis_prep` → level-1 + error |
| `02_explore_windata.pro` | windata の構造と座標系を調べる（eispac との対応確認） |
| `03/04_fit_box*.pro` | 22 輝線フィット（04 は線幅共有＋欠損マスク修正版） |
| `05b_ca17_ratios.pro` | Ca XVII 分離用の原子データ比を CHIANTI から取得 |
| `06_dump_spectra.pro` | 箱内平均スペクトルを書き出す（Python と突き合わせる用） |
| `08_fit_ca.pro` | **Ca 線を論文 p.6 の拘束どおりに解く**（Ko+2009 の Ca XVII 分離込み） |
| `09_gofnt.pro` | CHIANTI 9.0.1 で 22 輝線の G(T)（Feldman 1992 + chianti.ioneq） |
| `10_poa_check.pro` | PINTofALE が IDL 9.2 でコンパイルできるか確認（29/29 成功） |
| `11_calcurve.pro` | **打ち上げ後較正カーブ**（Warren+2014 NRL / Del Zanna 2013） |
| `13_mcmc_dem.pro` | **PINTofALE MCMC_DEM で DEM を出す**。単位のつじつまはコメント参照 |

## 到達済みの成果（2026-08-06）

### Python (eispac) 側
**論文 Table 2 を再現できた。** 箱 y=[244:274], x=[32:40] で
21 輝線中 13 本が論文の 15% 以内、median ratio 0.89、ばらつき 0.10 dex。

### IDL / SolarSoft 側（同日、環境が見つかったので実施）
**level-0 → eis_prep → 22 輝線フィット → MCMC_DEM まで通した。**

1. **IDL と eispac は 22 輝線すべてで一致**（小数第 2 位まで）。median 0.887。
   → 未解決だった食い違い（Fe XVI 0.62、S XIII 0.72、Si VII 0.39）は
     **eispac 固有の問題ではない**ことが確定。
2. **MCMC_DEM が論文 Table 2 の R 列を再現した。** Fe XIII だけ大きく外れ、
   Ar XIV / Ca XIV / Ca XV が揃って 1.3-1.5 になるという論文の特徴まで一致
   （論文 1.36/1.31/1.43、我々 1.37/1.35/1.46）。
   EM 分布のピークは logT=6.60 = **4.0 MK**（論文アブストラクトと一致）。
   → **絶対強度が 11% 低くても DEM 解析の結論は論文と同じ。**
3. **Ca XVII のブレンド分離ができた**（Ko et al. 2009 相当）。4.99 → 0.767。
4. **論文が書いていない inter-moss 箱の座標を Figure 1-3 から復元**（15 領域ぶん）。

### 準備の過程での重要な発見
1. **論文 Eq.(A1) の指数は誤植**。正しくは定数項つき 3 次式（実データで確定）
2. **eispac に Ca XVII 192.858 のブレンド分離テンプレートが無い**
   （SSW 側で正解値を作ったので、自作テンプレートの検証はできる）
3. **Si VII 275.368 が論文の 0.39 倍なのは未解決。** ただし範囲は絞れた。
   - 容疑を 10 個潰した（フィッター実装 / eispac 固有 / 箱の位置 /
     打ち上げ後較正 2 種 / 実効面積のバージョン / despike / 欠損値処理 /
     フィットの順番 / 未モデル化のブレンド / 「箱で説明できる」説）
   - **決め手は AIA Fe XVIII**。Si VII が合う箱は AIA Fe XVIII が
     論文の 0.161 倍（6 分の 1）で、論文 Table 2 の AIA 行と整合しない。
     逆に AIA が合う箱では Si VII が必ず 2.6 倍低い。
   - → **論文 Table 2 の Si VII 値は、同じ Table 2 の AIA Fe XVIII 値と
     整合しない**という形まで絞れた。著者に問い合わせる価値がある。
   - Fe XIII が DEM から 1.3-2.8 倍ずれるのは **論文側でも同じ**
     （Warren+2011 で R=1.87/1.90）。原子データの既知の問題。
   - 参照論文は `papers/refs/`（arXiv 版。`papers/fetch_papers.sh` で取得）。
4. **欠損値はサンプル単位でマスクすること（IDL / Python の両方）**。
   論文 §3 も "In computing these averaged profiles, missing data are not
   included" とわざわざ書いている。
   - IDL: `eis_getwindata` の欠損値は `-100`。怠ると最も明るい線
     （Fe XII 195.119）だけ 24% 小さくなる。
   - **Python: 欠損は NaN では入っていない。** level-0 の -100 に較正係数を
     掛けた**大きな負のフラグ値**（-1000〜-6000）なので、
     `np.isfinite` / `np.nansum` では**素通りする**。
     eispac が立てる **`cube.mask`（True = 使用禁止）で落とす**こと。
     エラーは出ず、平均が静かに下がるだけ。クイックルック図に
     **横方向の黒い縞**が出たらこれを疑う。
     （2026-08-06 に `fit_box_spectra.py` / `quicklook_raster.py` を修正。
       採用箱では median 0%、弱い線だけ数 % 動いた）

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
python scripts/fetch_data.py            # 既定 = region 7 (論文 Table 2 の活動領域)

# 3. 動作確認：ここまでが再現済みの到達点
python scripts/quicklook_raster.py data/eis/eis_20110702_030712.data.h5 figures/q.png
python scripts/make_fe18_map.py       data/eis/eis_20110702_030712.data.h5
python scripts/coalign_eis_aia.py     data/eis/eis_20110702_030712.data.h5
python scripts/compare_table2.py      data/eis/eis_20110702_030712.data.h5 "244:274,32:40"
#   -> 21 輝線中 13 本が論文 Table 2 の 15% 以内、median ratio 0.89 になれば成功
```

`.claude/settings.json` の承認ルールにはこのマシン固有の
miniforge 絶対パスが入っているが、`Bash(python *)` 等の一般パターンも
入れてあるので別マシンでもそのまま使える。

## ★ IDL / SSWIDL でやるべきこと（**このマシンで全部できる**）

**訂正 (2026-08-06)**: 「この計算機には IDL / SSW が無い」は**誤り**だった。
実際には全部揃っている。下の「環境」節を参照。

| | 場所 |
|---|---|
| IDL 9.2.0 (NV5, 名古屋大学ライセンス) | `/usr/local/nv5/idl92`、`idl` は PATH 上 |
| SolarSoft (2021-04 版) | `/opt/ssw` → `/lustre/sc/ssw` |
| SSWDB | `/opt/sswdb` → `/lustre/sc/sswdb` |
| Hinode/EIS（較正データ 2.1 GB 込み） | `/opt/ssw/hinode/eis` |
| **PINTofALE** | `/opt/ssw/packages/poa`（放射率 DB 401 MB 込み） |
| CHIANTI 9.0.1 | `/opt/ssw/packages/chianti/dbase` (1.6 GB) |
| Warren 本人のルーチン群 | `/opt/ssw/hinode/eis/idl/atest/hwarren/` |

実行は `bash scripts/idl/run_sswidl.sh prog.pro logfile`。
lustre が遅く SSW の起動だけで 1–3 分かかるので、必ずバックグラウンドで回すこと。
ssw_batch は `.run` で実行する = メインプログラム扱いなので、
**`.pro` の末尾に `exit` と `end` の両方が必須**
（`end` が無いと "End of file encountered before end of program" で落ちる）。
`.run` なので複数行の begin/endfor は普通に書ける。
未定義のシステム変数（例: CHIANTI 未ロードでの `!xuvtop`）は
**コンパイルエラー**になるので、必要な instrument を `SSW_INSTR` に入れておくこと。

Python (eispac) だけでは解決できず、SSW が要る項目。
**2026-08-06 に 1〜4 をすべて実施した。結果は `docs/00_log.md`。**

1. **`eis_auto_fit` と eispac のフィット結果を突き合わせる** … ✅ **完了**
   → **22 輝線すべてで一致した**（小数第 2 位まで）。median 0.887。
     未解決だった Fe XVI 0.62 / S XIII 0.72 / Si VII 0.39 は
     **eispac 固有ではなく実データにそう出ている**ことが確定。
   → 成果物: `work/idl_intensities_tied.csv`, `work/idl_intensities_box246.csv`

2. **打ち上げ後の EIS 較正カーブを作る** … ✅ **完了**
   `scripts/idl/11_calcurve.pro` → `work/eis_calcurve_20110702.txt`
   （波長, NRL(Warren+2014) 倍率, Del Zanna(2013) 倍率。SW 170-212 / LW 246-291）
   → **ただし食い違いの原因ではなかった**。どちらの較正を仮定しても
     観測された ratio のパターンと整合しない。

3. **Ca XVII 192.858 のブレンド分離を SSW でやる** … ✅ **完了**
   `scripts/idl/08_fit_ca.pro`（Ko et al. 2009 相当）→ **4.99 → 0.767**。
   原子データは `scripts/idl/05b_ca17_ratios.pro` で CHIANTI 9.0.1 から取得:
   O V 多重線の分岐比、Fe XI 192.813/188.216 = 0.20896（エネルギー比）。
   → eispac 用の自作テンプレートを検証する「正解値」ができた。

4. **PINTofALE の MCMC_DEM** … ✅ **完了**
   自前ダウンロードの 2015 年版が IDL 9.2 でそのまま動く（依存 29 本が 29/29 通る）。
   `scripts/idl/13_mcmc_dem.pro` で DEM を計算。
   → **論文 Table 2 の R 列のパターンまで再現**（Fe XIII だけ大きく外れ、
     Ar XIV / Ca XIV / Ca XV が揃って 1.3-1.5）。
   → EM ピーク **logT 6.60 = 3.98 MK**（論文「4 MK 付近」）、α = 3.2-3.5（論文 2.9）。
   ※ SSW 同梱の `/opt/ssw/packages/poa` は `fitting/` のパーミッションが
     `drwxr--r--` で入れず、しかも 2004 年版。**使えない。**
     `/scr/a000/c0234hotta/PINTofALE`（2015 年版）を使うこと。

5. **CHIANTI IDL と fiasco の寄与関数を突き合わせる** … ✅ **完了**
   `scripts/gofnt_fiasco.py`。fiasco の `ascii_dbase_root` を SSW 同梱の
   CHIANTI 9.0.1 に向けて**同じ原子データファイル**を読ませた。
   → **22 輝線すべてで median 1.000、最大でも 3.5% 差、形成温度は完全一致。**
   → **講習会は Python (fiasco) だけで寄与関数を出せる。**
   ※ 最大の差が Fe XIII の 2 本（+3.5%/-2.8%）に出るのが示唆的。
     Fe XIII は DEM で R=1.3-2.8 と外れる唯一のイオンでもある。

**→ ★リストは全項目完了。**

### 追加でできたこと（当初の★リストに無かったもの）

- **AIA 94 Å の Fe XVIII 専用応答関数**を作った（論文 p.6 が要求するもの）。
  `scripts/idl/15_aia94_fe18_resp.pro` + `scripts/aia94_fe18_response.py`。
  R(T) ピーク 2.73e-27 DN cm^5 s^-1 pix^-1 at logT 6.90。
- **論文が書いていない inter-moss 箱の座標を Figure 1-3 から復元**（15 領域ぶん）。
  `scripts/extract_paper_boxes.py`。region 13 だけ箱が 2 つ出るなど、
  論文の記述と一致することで正しさを確認済み。
- **AIA Fe XVIII で箱を独立に検証**（`scripts/check_fe18_box.py`）。
  JSOC の synoptic アーカイブ（登録不要、1 枚 1 MB）を使う。
  → EIS 0.89 と AIA 0.93 が揃って低い = 装置でも処理でもなく箱の位置の差。

## 環境（2026-08-06 時点で確認済み）

- **解析に使う env（これを使うこと）**:
  `/scr/a000/c0234hotta/home/miniforge3/envs/eis/bin/python`
  - eispac 0.99.4, fiasco 0.8.2, sunpy 8.0.0, aiapy 0.12.1, demregpy 1.0.0,
    ChiantiPy 0.16.0, numpy 2.5.1
  - ★ 2026-08-06 訂正: 以前書いていた `/home/sc/c0234hotta/miniforge3/envs/eis`
    は**現在は解決しない**（`/home/sc/...` が見えない）。上のパスが実体。
- base: `/cidashome/sc/c0234hotta/miniforge3/bin/python` (3.12.13)
  - numpy/scipy/astropy/sunpy/matplotlib/h5py はあるが **eispac 等は無い**
- IDL / SolarSoft: **ある**（2026-08-06 に判明。以前の「無い」は誤り）
  - IDL 9.2.0 `/usr/local/nv5/idl92`（名古屋大学ライセンス、`/usr/local/bin/idl`）
  - SSW `/opt/ssw`、SSWDB `/opt/sswdb`。`$SSW` は未設定なので自分で export する
  - `scripts/idl/run_sswidl.sh` がラッパ
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

→ **この活動領域を「答え合わせのできる題材」として講習会の主教材にするのが自然。**

同様に region 8（2010-07-23 14:32:10, NOAA 1089 = Warren et al. 2011 と同じ AR）も
`eis_20100723_143210.{data,head}.h5` が存在することを確認済み。
