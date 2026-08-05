# 環境構築（受講者に配る想定）

Hinode/EIS の解析に **IDL / SolarSoft は不要**です。Python だけで完結します。

---

## A. Google Colab（推奨・3 分）

ブラウザだけで動きます。何もインストールされていない状態から始められるので、
講習会当日に環境で詰まる事故がありません。

ノートブックの先頭で:

```python
!pip -q install eispac demregpy aiapy
# sunpy, ndcube, astropy, numpy, scipy, matplotlib, h5py は
# eispac / aiapy の依存として一緒に入る
```

インストール後に **ランタイムの再起動**を求められたら従ってください。

確認:
```python
import eispac, sunpy, aiapy, demregpy
print(eispac.__version__, sunpy.__version__, aiapy.__version__, demregpy.__version__)
```

### Colab を使うときの注意

| 事項 | 内容 |
|---|---|
| データ量 | EIS Level-1 は 1 ラスター約 **94 MB**、AIA/HMI 4 枚で約 **46 MB**。合計 140 MB 程度。数分で落ちる |
| セッション | 切れると `/content` の中身は消える。Google Drive をマウントして保存すると安全 |
| メモリ | 標準ランタイム（約 12 GB）で足りる。フルディスク AIA を 3 波長同時に扱うときだけ注意 |
| 速度 | **全ラスターの輝線フィットは避ける**（1 輝線で 95 秒）。論文どおり「箱の中で平均してからフィット」なら 22 輝線で 8 秒 |
| CHIANTI | フルデータベースは 1–2 GB あるので毎回は落とせない。**事前計算した寄与関数ファイルを配布**する（モジュール 6） |
| 日本語 | matplotlib の既定フォントに CJK が無い。**図のラベルは英語**にする（太陽物理の慣習にも合う） |

Google Drive に保存する場合:
```python
from google.colab import drive
drive.mount('/content/drive')
DATA = '/content/drive/MyDrive/eis_workshop'
```

---

## B. ローカル（conda）

```bash
mamba create -n eis -c conda-forge -y \
    python=3.12 numpy scipy matplotlib astropy sunpy ndcube h5py pandas \
    jupyterlab ipywidgets tqdm zeep drms parfive
mamba activate eis
pip install eispac demregpy aiapy
# モジュール 6/7（寄与関数を自分で作る場合）のみ
pip install fiasco
```

### この計算機での実績（2026-08-06 確認）

- 環境: `/home/sc/c0234hotta/miniforge3/envs/eis`
- eispac 0.99.4 / sunpy 8.0.0 / aiapy 0.12.1 / demregpy 1.0.0 / ndcube 2.4.0 /
  astropy 8.0.1 / numpy 2.5.1 / scipy 1.18.0 / fiasco 0.8.2
- 実行例: `/home/sc/c0234hotta/miniforge3/envs/eis/bin/python scripts/....py`

### pip だけで入れる場合（venv）

```bash
python3 -m venv eis-venv
source eis-venv/bin/activate
pip install eispac demregpy aiapy jupyterlab
```
Colab と同じ構成になるので、教材の互換性を確かめるのに使えます。

---

## C. データの取得

### EIS Level-1 (HDF5)

NRL が公開している Level-1 HDF5 を直接使います。
自分で `eis_prep` 相当をやる必要はありません。

```python
import eispac
eispac.download_hdf5_data(filename='eis_20110702_030712', local_top='data/eis')
```

または直接:
```
https://eis.nrl.navy.mil/level1/hdf5/2011/07/02/eis_20110702_030712.data.h5   (94 MB)
https://eis.nrl.navy.mil/level1/hdf5/2011/07/02/eis_20110702_030712.head.h5   (421 KB)
```

観測を探すときは EIS Science Nuggets や
[EIS の観測データベース](https://eis.nrl.navy.mil/) を使う。
**DEM 解析をしたいなら Ca XIV–XVII が含まれる観測プログラムを選ぶこと。**

### SDO/AIA, HMI

**VSO 経由（ユーザ登録不要）**を既定にします。

```python
import astropy.units as u
from sunpy.net import Fido, attrs as a
q = Fido.search(a.Time('2011-07-02T03:37:53','2011-07-02T03:38:23'),
                a.Instrument.aia, a.Wavelength(94*u.AA))
files = Fido.fetch(q[0][0:1], path='data/sdo/{file}')
```

- JSOC 経由 (`a.jsoc.Series` + `a.jsoc.Notify`) なら切り出し (cutout) が
  使えて軽いが、**JSOC に登録したメールアドレスが必要**。
  講習会では登録の手間が事故のもとなので使わない。
- **VSO は時々 staging に失敗する**。`Fido.fetch` の戻り値の `.errors` を
  確認してリトライすること（`scripts/download_sdo.py` に実装済み）。
- 検索結果には「1 枚 65 MB」と出るが、実体は Rice 圧縮 FITS で 7–15 MB。

---

## D. よくあるトラブル

| 症状 | 原因と対処 |
|---|---|
| `Fido.fetch` が空リストを返す | VSO の一時的な失敗。リトライする |
| AIA の引き算で shape が合わない | level-1.5 化しても波長ごとに配列サイズが違う（94/193 は 4096²、171 は 4094²）。`reproject_to` で共通グリッドに揃える |
| 輝線強度が論文と全然違う | テンプレートの成分番号を間違えている。`template['line_ids']` を確認する |
| フィットが終わらない | 全ラスターをフィットしている。箱の中で平均してからフィットする |
| `no name guard was found` 警告 | multiprocessing 用。スクリプトを `if __name__ == '__main__':` で囲むか `ncpu=1` にする |
| 図の日本語が □ になる | matplotlib のフォントに CJK が無い。ラベルは英語にする |
| Colab でセッションが切れた | Google Drive をマウントして中間ファイルを保存しておく |
