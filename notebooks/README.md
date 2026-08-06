# 講習会ノートブック

受講者は **GitHub の URL を Colab に食わせるだけ**で開けます。アップロード不要:

```
https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/01_raster.ipynb
                                        └─ owner/repo ─┘        └─ branch ─┘ └─ path ─┘
```

| # | ノート | 内容 | Colab で開く |
|---|---|---|---|
| 0 | `00_setup.ipynb` | 環境構築とデータ取得 | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/00_setup.ipynb) |
| 1 | `01_raster.ipynb` | EIS のデータを見る | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/01_raster.ipynb) |
| 2 | `02_fitting.ipynb` | スペクトル線フィット | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/02_fitting.ipynb) |
| 3 | `03_aia_fe18.ipynb` | AIA 94 → Fe XVIII | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/03_aia_fe18.ipynb) |
| 4 | `04_coalign.ipynb` | 座標合わせと箱の選択 | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/04_coalign.ipynb) |
| 5 | `05_table2.ipynb` | **論文 Table 2 と答え合わせ** | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/05_table2.ipynb) |
| 6 | `06_gofnt_emloci.ipynb` | 寄与関数と EM loci | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/06_gofnt_emloci.ipynb) |
| 7 | `07_dem.ipynb` | DEM インバージョン | [開く](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/07_dem.ipynb) |

- 半日コース: 0 → 5 ／ 1 日コース: 0 → 7
- **リポジトリが public でないと Colab の GitHub ローダーは開けません。**
  private のままにするなら、`.ipynb` を各自 Colab にアップロードしてもらうことになります。

## 受講者に最初に伝えること（Colab の保存の挙動）

GitHub から開いたノートは**読み取り専用の一時セッション**として扱われる。

- **編集も実行結果も保存されない。** 残すには「ファイル → ドライブにコピーを保存」
  （`MyDrive/Colab Notebooks/` にコピーができる。以降は教材側の更新は反映されない）
- **仮想マシンが消えるとダウンロードしたデータも消える**
  （放置で数十分、使っていても無料枠では最長で半日程度）
- → **EIS の 94 MB は毎セッション落とし直し**。だから各ノートが自分で取得する設計にしてある
- 落とし直しを避けるなら `drive.mount('/content/drive')` で Drive に置く手もあるが、
  許可ダイアログを踏ませることになるので既定にはしていない

## 設計上の約束（編集する人向け）

### 1 冊ずつ独立して動く

**Colab はノートブック 1 冊ごとに新しい VM が立ち上がる。**
前のノートが作ったファイルは残らないし、セッションが切れても消える。そこで:

- どのノートも先頭に**ブートストラップのセル**がある
  （リポジトリが無ければ `git clone` して `cd`、`sys.path` に `scripts` を追加）
- 観測データは**そのノートで必要な分だけ**その場で取得する（既にあれば何もしない）
- モジュール間の受け渡し（`work/box_intensities.csv`、
  `data/cache/aia_on_eis_grid.npz`）は、無ければ
  **`scripts/workshop.py` がその場で作り直す**（10〜20 秒）

→ 受講者が途中のモジュールから始めても、遅れて追いついても動く。

### `.ipynb` を手で書かない

`notebooks/src/NN_*.py` に `# %%` 区切りの**素の Python** として書き、

```bash
python notebooks/build_notebooks.py     # src/*.py -> *.ipynb
```

で生成する。JSON を直接いじると差分もレビューも辛いため。

### 編集したら必ず検証する

```bash
python notebooks/verify_notebooks.py        # 全ノートのコードセルを上から順に exec
python notebooks/verify_notebooks.py 03 04  # 番号で絞る
```

Colab 専用セル（`!pip` などの shell/magic を含むもの）は自動で飛ばし、
飛ばした番号を表示する（何を検証していないかが分かる）。

### 図のラベルは英語

Colab に日本語フォントが入っていないため、日本語ラベルは豆腐（□）になる。
**説明は日本語、図は英語**で統一している。

## 所要時間の目安（Colab の標準 VM、実測ベース）

| ノート | 重い処理 |
|---|---|
| 00 | EIS 94 MB のダウンロード（数分） |
| 02 | 部分ラスターのフィット（約 45 秒）＋ 22 輝線の箱フィット（約 5 秒） |
| 03, 04 | AIA 3 MB のダウンロード＋再投影（十数秒） |
| 05 | 箱を 3 通り比較（約 20 秒） |
| 06, 07 | 事前計算した G(T) を読むので数秒 |

## ライセンス

この教材は **CC BY 4.0**（リポジトリの [LICENSE](../LICENSE)）。
出典を示せば講義・講習会でそのまま使えます。
