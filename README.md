# Hinode/EIS データ解析講習会

**Warren, Winebarger & Brooks (2012), ApJ 759, 141**
*"A Systematic Survey of High-Temperature Emission in Solar Active Regions"*
と同じ解析を、**受講者が自分の手で最後まで通せるようになる**ための教材です。

Python のみ（IDL / SolarSoft 不要）。**Google Colab で完結**します。

## なぜこの論文なのか

論文 Table 2 に **region 7（2011-07-02, NOAA 1243）の 22 輝線 + AIA 94 Å の観測値**が
そのまま載っています。**論文に数値表がある唯一の活動領域**なので、
自分の解析結果を 1 行ずつ答え合わせできます。
「できたつもり」で終わらない教材になります。

## はじめ方

**Colab のリンクを踏むだけです**（インストール不要、Google アカウントだけあれば OK）:

### → [Colab で開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop.ipynb)

演習の答えは[こちら](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop_answers.ipynb)（実行結果つき。開くだけで読めます）。

**5 章が 1 冊に入っています。** 上から順に実行してください（約 2 時間 45 分）。

| 章 | 内容 |
|---|---|
| 1 | EIS のデータを見る |
| 2 | フィットして**強度**を出す |
| 3 | **速度**を出す |
| 4 | 線幅から**非熱的速度**を出す |
| 5 | **温度分布 (DEM)** を出す |
| 付録 | 自分の研究で使うときに読む |

**狙いは Solar-C (EUVST) の準備です。** EUVST の EUV バンド 170–215 Å は
EIS の短波長帯とほぼ同じで、**同じ輝線を撮ります**。
今日やることは、そのまま 2028 年に使えます。

セッションが切れても、上から流し直せば復帰します
（既に取得したファイルは落とし直しません）。

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

- **速度ゼロは自分で決めるもの**。EIS に絶対的な波長基準は無い
- 軌道変動とスリット傾きで **53 km/s ぶん**波長が動く（測りたい信号より大きい）
- **装置幅を引かないと非熱的速度は出ない**。装置幅はスリット上の位置で変わる
- **欠損値は NaN ではない**。大きな負のフラグ値なので `np.nansum` を素通りする（付録 A）
- 同じ CHIANTI を読んでも **実装によって寄与関数が 12.6 倍ずれる**（付録 D）
- **DEM の傾きは手法で動く**が、**EM ピーク 4 MK は動かない**（付録 C）

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
