# Hinode/EIS データ解析講習会

**Warren, Winebarger & Brooks (2012), ApJ 759, 141**
*"A Systematic Survey of High-Temperature Emission in Solar Active Regions"*
と同じ解析を、**受講者が自分の手で最後まで通せるようになる**ための教材です。

Python のみ（IDL / SolarSoft 不要）。**Google Colab で完結**します。

## なぜこの論文なのか

論文 Table 2 に **region 7（2011-07-02, NOAA 1243）の 22 輝線 + AIA 94 Å の観測値**が
そのまま載っています。**論文に数値表がある唯一の天体**なので、
自分の解析結果を 1 行ずつ答え合わせできます。
「できたつもり」で終わらない教材になります。

## はじめ方

Colab のリンクを踏むだけです（インストール不要）:

| # | ノート | 内容 | |
|---|---|---|---|
| 0 | 環境構築とデータ取得 | pip とデータ取得（登録不要） | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/00_setup.ipynb) |
| 1 | EIS のデータを見る | 輝線ごとに別物に見える | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/01_raster.ipynb) |
| 2 | スペクトル線フィット | 強度はフィットの産物 | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/02_fitting.ipynb) |
| 3 | AIA 94 → Fe XVIII | 論文 Appendix の誤植を自分で確かめる | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/03_aia_fe18.ipynb) |
| 4 | 座標合わせと箱の選択 | どこを測るかを決めるのが解析の本体 | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/04_coalign.ipynb) |
| 5 | **論文 Table 2 と答え合わせ** | 講習会の山場 | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/05_table2.ipynb) |
| 6 | 寄与関数と EM loci | 1/(4π) の罠 | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/06_gofnt_emloci.ipynb) |
| 7 | DEM インバージョン | ill-posed と手法依存性 | [Colab](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/07_dem.ipynb) |

- **半日コース**: 0 → 5 ／ **1 日コース**: 0 → 7
- 各ノートは**単独で動きます**。途中から始めても、セッションが切れても大丈夫です。

ローカルで動かす場合:

```bash
git clone https://github.com/hottahd/EIS_practice.git
cd EIS_practice
mamba create -n eis -c conda-forge -y python=3.12 numpy scipy matplotlib \
      astropy sunpy ndcube h5py pandas jupyterlab tqdm zeep drms parfive
mamba activate eis
pip install eispac demregpy aiapy fiasco
python scripts/fetch_data.py          # 観測データ（EIS 94 MB + SDO）
```

## この教材が大事にしていること

**一致したことより、合わない理由を説明できることが実力。**
準備の過程で実際に踏んだ罠を、そのまま教材にしています:

- 論文 **Eq.(A1) の指数は誤植**（実データを通すと分かる）
- **欠損値は NaN ではない**。大きな負のフラグ値なので `np.nansum` を素通りする
- 同じ CHIANTI を読んでも **実装によって寄与関数が 12.6 倍ずれる**（1/4π の約束の違い）。エラーは出ない
- **DEM の傾き α は手法で 1.5–2.3 と動く**。一方で **EM ピーク 4 MK はどの手法でも動かない**
- 論文 Table 2 の **Si VII の値は、同じ表の AIA Fe XVIII の値と整合しない**（未解決として正直に扱う）

## リポジトリの構成

| | |
|---|---|
| `notebooks/` | 講習会ノートブック（`src/*.py` から生成。[README](notebooks/README.md)） |
| `scripts/` | 検証済みの解析スクリプト（Python）と `idl/`（SSW/IDL 版） |
| `docs/` | [作業ログ](docs/00_log.md) ・ [論文の手法分解](docs/01_paper_analysis.md) ・ [カリキュラム](docs/02_workshop_plan.md) ・ [環境構築](docs/03_environment.md) ・ [理論屋向けの物理解説](docs/04_physics_primer.md) |
| `work/` | 事前計算した寄与関数・AIA 応答・参照値（Colab で CHIANTI を落とさずに済ませるため） |
| `papers/` | 論文の出典と取得スクリプト（**PDF 本体は含みません**。[README](papers/README.md)） |

観測データ（`data/`）はリポジトリに含みません。`scripts/fetch_data.py` か
各ノートが自動で取得します。

## ライセンスと帰属

この教材（ノートブック、`docs/`、`scripts/`、`figures/`、`work/` の計算結果）は
**[CC BY 4.0](LICENSE)** です。**出典を示せば**自由に利用・改変・再配布できます。
講習会や講義でそのまま使っていただいて構いません。

> Hotta, H., *Hinode/EIS データ解析講習会 教材*, https://github.com/hottahd/EIS_practice (CC BY 4.0)

**この license が及ばないもの**:

- **論文 PDF は含まれていません**（再配布の許諾が無いため）。
  出典と取得方法は [`papers/README.md`](papers/README.md)。
- `docs/01_paper_analysis.md` などに転記した**論文の数値表は各論文に帰属**します
  （事実・データの、出典明記つきの引用として掲載）。
- `work/gofnt_chianti901*.txt` などは **CHIANTI 原子データベース**から計算したものです。
  利用の際は CHIANTI を引用してください
  (Dere et al. 1997, A&AS 125, 149 / Dere et al. 2019, ApJS 241, 22)。
- Hinode/EIS と SDO/AIA の観測データは各ミッションのデータポリシーに従います。

なお Creative Commons はソフトウェアへの CC ライセンスの適用を推奨していません。
`scripts/` を別ライセンス（MIT など）にしたい場合はご連絡ください。
