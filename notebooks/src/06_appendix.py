# %% [markdown]
# # 付録
#
# **講習会では走らせません。** 自分の研究で EIS / EUVST を使うときに読んでください。
#
# 本編で使っているコードは、ここに書いてあることを**既に正しく処理**しています。
# 「なぜそうしているのか」を知りたくなったときの参照先です。
#
# | | 内容 |
# |---|---|
# | A | 欠損値は NaN ではない |
# | B | 誤差の入れ方 |
# | C | DEM をどこまで信じてよいか |
# | D | G(T) と 1/(4π) の罠 |
# | E | 較正 |
# | F | AIA 94 Å から Fe XVIII を分離する |
# | G | Ca XVII のブレンド分離 |
# | H | 論文 Table 2 との全面照合 |
# | I | 他の活動領域・15 領域の統計 |
# | J | ラスターは同時刻の画像ではない |
# | **K** | **eispac の使い方（早見表）** |

# %% [markdown]
# ## 付録 A: 欠損値は NaN ではない
#
# **EIS の level-1 では、不良画素・宇宙線で捨てられたサンプルは
# 大きな負のフラグ値**（level-0 の `-100` に較正係数を掛けたもの。
# 実際には −1000 〜 −6000 程度）として入っています。NaN ではありません。
#
# したがって:
#
# - `np.isnan` / `np.nansum` では**素通りします**
# - 足すと大きな負の数が入り、その行だけ暗くなります（強度マップに**横縞**が出る）
# - **エラーは一切出ません**
#
# eispac は `cube.mask`（True = 使ってはいけない）を立ててくれるので、これで落とします:
#
# ```python
# d = np.where(np.asarray(cube.mask, dtype=bool), np.nan, cube.data)
# img = np.nanmean(d, axis=2) * d.shape[2]
# ```
#
# 論文 §3 も *"In computing these averaged profiles, missing data are not included"*
# とわざわざ書いています。
#
# **影響の実測**（本編の箱、22 輝線）: median は 0%、大半の線で 1% 未満。
# ただし弱い線では効きます（Ca XVII +6.9%、Ca XVI −4.6%、Ar XIV +3.2%）。
# IDL 側では最も明るい Fe XII 195.119 が 24% 小さくなった例もあります。
#
# **EUVST でも同じ発想が要ります。** 「欠損をどう表現しているか」は
# データ形式ごとに違うので、**必ず確認してから平均を取る**こと。

# %% [markdown]
# ## 付録 B: 誤差の入れ方
#
# 箱の中で数百画素を平均すると、**統計誤差は 0.2% 程度**まで落ちます。
# 一方、論文が使っている誤差は **22%** です。
#
# 差は**絶対較正の不確かさ**（系統誤差）で、これは平均しても減りません。
#
# 統計誤差だけで DEM を解こうとすると、モデル（滑らかな DEM）の
# わずかな不完全さが全部 χ² に化けて**発散**します。
#
# ```python
# edn = np.maximum(stat_err, 0.22 * intensity)   # 較正誤差を「床」として入れる
# ```
#
# **「誤差が小さい」ことは良いことではありません。**
# 何の誤差を見積もっているのかを意識してください。

# %% [markdown]
# ## 付録 C: DEM をどこまで信じてよいか
#
# DEM は**第一種 Fredholm 積分方程式**の逆問題で、本質的に ill-posed です。
#
# **実質的な自由度は輝線の本数より少ない。**
# 22 輝線 × 33 温度ビンの応答行列を特異値分解すると、
# 最大の 1/1000 以上ある特異値は **12 本**しかありません
# （G(T) が幅広く重なっているため）。
#
# ```python
# sv = np.linalg.svd(G / (4*np.pi), compute_uv=False)
# (sv > 1e-3 * sv[0]).sum()      # → 12
# ```
#
# **手法・設定で答えがどれだけ動くか**（同じ強度・同じ G(T) で実測）:
#
# | 設定 | chi2 | EM ピーク | 傾き α (logT 6.0–6.6) |
# |---|---:|---|---:|
# | demregpy 既定 | 3.43 | 6.65 (4.5 MK) | 1.86 |
# | demregpy reg_tweak=2 | 4.48 | 6.65 | 1.55 |
# | demregpy（MCMC を事前分布に） | 1.49 | 6.60 (4.0 MK) | 2.22 |
# | MCMC_DEM (PINTofALE) | — | 6.60 (4.0 MK) | 2.30 |
# | 論文 Table 1 | — | ~4 MK | 2.9 |
#
# → **ピーク温度 4 MK はどの手法でも動かない**（頑健な結論）。
# **傾き α は 1.5–2.3 と動く**（慎重に扱うべき量）。
# 加熱の頻度を α で議論するなら、この系統誤差を押さえる必要があります。
#
# **単位の罠**: DEM の単位は実装で違います。
#
# | | DEM の単位 | ビンあたりの EM |
# |---|---|---|
# | PINTofALE | cm⁻⁵ **/ logK** | DEM × ΔlogT |
# | demregpy | cm⁻⁵ **/ K** | DEM × T ln10 ΔlogT |
#
# 揃えないと 10⁷ ずれます。**傾きも +1 ずれます**（実際に間違えました）。

# %% [markdown]
# ## 付録 D: G(T) と 1/(4π) の罠
#
# 同じ CHIANTI 9.0.1 のファイルを 3 つの実装に読ませて G(T) を比べた結果:
#
# | 実装 | 単位の約束 | CHIANTI IDL との比 |
# |---|---|---|
# | CHIANTI IDL `emiss_calc` | 4π で割らない | 1.000 |
# | fiasco `contribution_function` | 4π で割らない | 0.972 – 1.035 |
# | **ChiantiPy `ion.emiss()`** | **sr⁻¹（4π で割ってある）** | **0.077 – 0.082** |
#
# **1/(4π) = 0.0796。12.6 倍ずれます。しかもエラーは出ません。**
# 全輝線が一律にずれるので、線ごとの比を見ている限り気づけません。
#
# **気づく方法はオーダーの検算だけ**です:
#
# | 量 | 覚えておく値 |
# |---|---|
# | Fe XII 195.119 の G ピーク | ~1.4×10⁻²³ erg cm³ s⁻¹ |
# | 活動領域の EM（視線積分） | 10²⁷ – 10²⁹ cm⁻⁵ |
# | コロナの電子密度 | ~10⁹ cm⁻³ |
# | → 視線長 L = EM/n_e² | ~100 Mm（妥当） |
#
# 残る 3% の差は**準安定準位の占有数**の扱いによるもので、
# **密度敏感線（Fe XIII 202/204 など）に集中**しています。
# Fe XIII が DEM で最も外れる線であることの独立な傍証になっています。
#
# 自分で計算するなら **fiasco** が扱いやすい（CHIANTI IDL と 3% 以内で一致）。
# 実装は `scripts/gofnt_fiasco.py`。
#
# **同じ話が EIS のデータ自体にもあります。** `cube.unit` は
# `erg / (s sr cm2)` と表示されますが、中身は**波長 1 Å あたり**の値です。
# 波長方向に足しただけでは輻射強度にならず、**波長刻み（0.0223 Å）を掛けて**
# 初めて erg cm⁻² s⁻¹ sr⁻¹ になります
# （フィットの `int` は $A\sigma\sqrt{2\pi}$ なので、こちらは掛かっています）。

# %% [markdown]
# ## 付録 E: 較正
#
# EIS は 2006 年打ち上げで、**有機物の付着などで感度が落ちています**。
# しかも**波長によって落ち方が違います**。
#
# 打ち上げ後較正が 2 つ提案されていますが（Del Zanna 2013、Warren et al. 2014）、
# **両者は食い違います**。どちらを使ったかで強度が数十 % 変わります。
#
# → **論文には必ず「どの較正を使ったか」を書く。**
# 他人の値と比べるときは、まず較正を揃える。
#
# 実装は `scripts/idl/11_calcurve.pro`、結果は `work/eis_calcurve_20110702.txt`。
#
# **絶対較正は 20% 程度ずれるもの**、という感覚を持っておくとよいです。
# EUVST でも同じ問題は必ず起きます。

# %% [markdown]
# ## 付録 F: AIA 94 Å から Fe XVIII を分離する
#
# EIS で観測できる最高温の強い輝線は Ca XVII 192.858（~5 MK）です。
# **それより上を拘束するものが無い**と、EM 分布の高温側が決まりません
# （第 5 章の EM loci で、log T > 6.9 の上限が跳ね上がっていたのがこれ）。
#
# そこで **AIA 94 Å の Fe XVIII 93.932 Å（~7 MK）**を使います。
# ただし 94 Å チャンネルは低温の線に汚染されているので、
# 171 Å と 193 Å から「低温成分」を経験的に見積もって引きます（論文 Appendix）。
#
# ```
# x        = (0.31 I_171 + 0.69 I_193) / 116.54
# I_94warm = 0.39 (a1 + a2 x + a3 x^2 + a4 x^3)
# a        = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]
# I_FeXVIII = I_94 - I_94warm
# ```
#
# **★ 論文に印刷されている式の指数は誤植です。** 字面どおり
# $\sum a_i x^i$ と読むと warm 成分が観測値を桁違いに超えます。
# 正しくは**定数項つきの 3 次式**（実データで確認）。
#
# 適用限界: フレア中は不可（Fe XXIV が 193 Å に入る）、
# 明るい moss でも破綻、AIA の感度劣化補正を掛けてはいけない。
#
# DEM に入れるときは、**公式の 94 Å 応答は使えません**（低温線込みのため）。
# Fe XVIII だけの応答関数を作ってあります（`work/aia94_fe18_response.txt`、
# ピーク 2.73×10⁻²⁷ DN cm⁵ s⁻¹ pix⁻¹ at log T 6.90）。
#
# 実装は `scripts/aia_fe18.py`, `scripts/aia94_fe18_response.py`。

# %% [markdown]
# ## 付録 G: Ca XVII のブレンド分離
#
# 第 2 章の演習 3 で Ca XVII 192.858 が論文の 5 倍になったのは、
# **Fe XI 192.813 と O V の複合線に埋もれている**ためです。
# eispac 同梱の `ca_17_192_858.1c` は、この波長域を単一ガウシアンで塗るだけで
# ブレンドを分離しません。
#
# 論文は Ko et al. (2009) の方法で分離しています。
# **自作テンプレート**で同じことができます（`scripts/ca17_template.py`）:
#
# | | Ca XVII 192.858 | 論文比 |
# |---|---:|---:|
# | eispac 同梱（1 成分） | 198.5 | 4.75 |
# | **自作 5 成分テンプレート** | **31.3** | **0.75** |
# | SSW/IDL 版 | 32.0 | 0.77 |
#
# 設計の要点:
#
# - **分離できない成分は統合する。** EIS の波長分解能（0.059 Å）に対し、
#   O V の 192.797/192.801 は分離不能 → 強度重み付き波長で 1 本にまとめる
# - 自由パラメータは 5 つだけ。他は**原子データと他の輝線から固定**
#   （O V の分岐比、Fe XI 192.813 = Fe XI 188.216 × 0.20896、
#     Ca XVII の線幅 = Ca XIV 193.874 の線幅）
#
# **成分を増やせば良くなるわけではない**、というのがこの作業の教訓です。

# %% [markdown]
# ## 付録 H: 論文 Table 2 との全面照合
#
# 22 輝線すべてを論文と比べると、**median 0.89、21 本中 13 本が 15% 以内**
# （箱 y=[244:274] x=[32:40]）。論文の誤差 ±22% の中に収まります。
#
# **★ 箱の選び方の診断法**: ratio を**形成温度に対して**並べ、傾きを見ます。
#
# | 傾き | 意味 |
# |---|---|
# | ≈ 0 | 論文と同じ温度組成の場所を見ている |
# | < 0 | 暖かいループ寄りを選んでいる |
# | > 0 | 高温コア寄りを選んでいる |
#
# 較正は波長の関数であって温度の関数ではないので、
# **温度に沿ったパターンが出たら、それは場所の違い**です。
#
# **★ 要約統計 1 つで判断しない**（実測）:
#
# | 箱 | median | 傾き | 15% 以内 | ばらつき |
# |---|---:|---:|---:|---:|
# | inter-moss（採用） | 0.89 | +0.21 | 13/21 | 0.10 dex |
# | 論文の箱サイズ | 0.93 | +0.27 | 15/21 | 0.12 dex |
# | **適当に明るいところ** | **0.96** | **−0.34** | **5/21** | **0.26 dex** |
#
# 一番下は **median が最も 1 に近いのに、15% 以内は 5 本しかありません**。
# 個々の線が上下に外れて打ち消し合っているだけです。
#
# 実装は `scripts/compare_table2.py`。

# %% [markdown]
# ## 付録 I: 他の活動領域・15 領域の統計
#
# 論文 Table 1 には 15 の活動領域があります（`docs/01_paper_analysis.md` に全リスト）。
# `scripts/fetch_data.py --region N` で別の領域のデータを取れます。
#
# - region 8 (2010-07-23, NOAA 1089) は Warren et al. 2011 と同じ活動領域
# - region 1 (2010-06-19) は **Ca XIV–XVI が無い**
#   → 「観測プログラムによって使える輝線が違う」実例。3 MK 以上を拘束できない
#
# 論文の主要な図:
#
# - Figure 4: 高温放射 I_hot と磁束 Φ_M の関係（べき指数 2.3）
# - Figure 9: 暖かい成分の EM と Φ_M は**逆相関**
#
# 自分が興味を持っている活動領域で同じ解析をしてみるのが、次の一歩です。

# %% [markdown]
# ## 付録 J: ラスターは同時刻の画像ではない
#
# スリットを 1 ステップずつ動かして作るので、**x 軸は空間であると同時に時間軸**です。
# 今日の観測がどれだけ時間をかけているか、ヘッダから確かめられます。

# %%
import eispac

from workshop import EIS_FILE
from astropy.time import Time

c = eispac.read_cube(EIS_FILE, 195.119)
h = c.meta["index"]
t0, t1 = Time(h["date_obs"]), Time(h["date_end"])
total = (t1 - t0).to_value("s")

print("観測プログラム :", h["stud_acr"], " (提案者:", h["st_auth"] + ")")
print("開始 / 終了    :", h["date_obs"], "/", h["date_end"])
print(f"所要時間       : {total/60:.1f} 分  ({h['nraster']} ステップ)")
print(f"1 ステップ     : {total/h['nraster']:.0f} 秒")
print(f"視野           : {h['fovx']:.0f}″ x {h['fovy']:.0f}″"
      f"  （x は {h['fovx']/h['nraster']:.1f}″/step、スリット幅は {h['slit_id']}）")

# %% [markdown]
# **62 分かかっています。** 左端と右端では 1 時間離れています。
#
# - 時間変化する現象（フレア、ジェット）には使えない
# - ドップラー速度のマップも、同時刻の速度場ではない
# - 一方、定常的な構造を測るなら問題なく、むしろ S/N の面では有利
#
# **Solar-C EUVST はカデンス 1 秒**なので、この制約は大きく緩みます。
# ただし「ラスターは掃いて作る」こと自体は変わりません。

# %% [markdown]
# <a id="eispac"></a>
# ## 付録 K: eispac の使い方（早見表）
#
# 演習で手が止まったときの参照用です。**主な関数の引数を全部説明**してあります。
# （Python では `help(eispac.read_cube)` でも同じものが読めます）
#
# 公式ドキュメント: https://eispac.readthedocs.io/

# %% [markdown]
# ### 1. データを取る
#
# ```python
# from eispac.download import download_hdf5_data
# download_hdf5_data(filename="eis_20110702_030712", local_top="data/eis")
# ```
#
# | 引数 | 既定値 | 意味 |
# |---|---|---|
# | `filename` | — | EIS のファイル名。`eis_YYYYMMDD_HHMMSS` の形。リストも可 |
# | `local_top` | `'data_eis'` | 保存先のディレクトリ |
# | `source` | `'nrl'` | 取得元。`'nrl'` / `'nasa'` / `'mssl'` |
# | `datetree` | `False` | `True` にすると `YYYY/MM/DD/` の階層を作る |
# | `nodata` / `nohead` | `False` | データ本体 / ヘッダを落とさない |
# | `overwrite` | `False` | 既にあっても取り直す |
# | `max_conn` | `2` | 同時接続数 |
#
# **ユーザ登録は不要**です。この教材では `workshop.ensure_eis()` が同じことをしています。

# %% [markdown]
# ### 2. 何が入っているか見る
#
# ```python
# wininfo = eispac.read_wininfo(head_file)   # 引数はヘッダファイルのパスだけ
# ```
#
# EIS は全波長を降ろせないので、**必要な輝線の周りだけ**を切り出して観測します。
# 使いたい輝線が入っているかを、まずこれで確認します。
# 返り値は構造化配列で、`iwin`（ウィンドウ番号）、`line_id`、`wvl_min`、`wvl_max` を持ちます。

# %%
import eispac
import numpy as np

from workshop import EIS_FILE

wi = eispac.read_wininfo(EIS_FILE.replace(".data.h5", ".head.h5"))
print(f"{len(wi)} ウィンドウ。最初の 3 つ:")
for w in wi[:3]:
    print(f"  {w['iwin']:2d} {str(w['line_id']):<20} {w['wvl_min']:.3f} – {w['wvl_max']:.3f} Å")

# %% [markdown]
# ### 3. 読む
#
# ```python
# cube = eispac.read_cube(datafile, 195.119)
# ```
#
# | 引数 | 既定値 | 意味 |
# |---|---|---|
# | `filename` | — | データファイルかヘッダファイルのパス（どちらでもよい） |
# | `window` | `0` | **ウィンドウ番号、またはその中の波長**。`195.119` のように波長で指定できる |
# | `exp_set` | `'sum'` | 1 つのラスター位置で複数露出がある観測のときだけ効く。`'sum'` は全部足す |
# | `apply_radcal` | `True` | 打ち上げ前の較正を当てて物理単位にする。`False` なら**光子カウントのまま** |
# | `radcal` | `None` | 自分で作った較正カーブを渡す（付録 E の較正を試すときに使う） |
# | `abs_errs` | `True` | カウントの絶対値から誤差を出す。暗電流引き算で負になった画素にも妥当な誤差が付く |
# | `count_offset` | `None` | カウントに定数を足してから較正する（処理の検証用） |
#
# 返ってくる `EISCube` の中身:
#
# | | |
# |---|---|
# | `cube.data` | (y, x, 波長) の配列。**波長 1 Å あたり**の値（付録 D） |
# | `cube.wavelength` | 同じ形の波長配列。**軌道変動とスリット傾きは補正済み**（第 3 章） |
# | `cube.uncertainty.array` | 誤差 |
# | `cube.mask` | True = 使ってはいけないサンプル（付録 A） |
# | `cube.meta` | ヘッダ類（下記） |
#
# `cube.meta` でよく使うもの:
#
# | キー | 中身 |
# |---|---|
# | `index` | FITS ヘッダ相当（`date_obs`, `nraster`, `fovx`, `stud_acr` …） |
# | `pointing` | `xcen`, `ycen`, `x_scale`, `y_scale` |
# | `extent_arcsec` | `[x0, x1, y0, y1]`。`imshow(extent=...)` にそのまま渡せる |
# | `wave_corr` | 波長補正（適用済み。第 3 章） |
# | `slit_width` | 波長分解能 FWHM [Å]、スリット位置ごと（第 4 章） |
# | `duration` | 露出時間 [s]、ステップごと |

# %%
cube = eispac.read_cube(EIS_FILE, 195.119)
print("data      :", cube.data.shape, cube.unit)
print("meta のキー:", ", ".join(list(cube.meta.keys())[:8]), "...")

# %% [markdown]
# ### 4. テンプレートを選ぶ
#
# ```python
# path  = eispac.data.get_fit_template_filepath("fe_12_195_119.2c.template.h5")
# tmplt = eispac.read_template(path)      # 引数はパスだけ
# ```
#
# - `.1c` / `.2c` / `.3c` はガウシアンの本数
# - `tmplt.template["line_ids"]` … **どの成分が何の輝線か。必ず確認する**（第 2 章）
# - `tmplt.template["n_gauss"]`, `["n_poly"]` … 成分の数と背景の次数
# - `tmplt.parinfo` … パラメータごとの辞書。中身は下の表
# - `tmplt.central_wave` … `read_cube` に渡せる波長
#
# `parinfo`（各パラメータの設定。mpfit の作法）:
#
# | キー | 意味 |
# |---|---|
# | `value` | 初期値 |
# | `fixed` | 1 なら動かさない |
# | `limited` / `limits` | `[下限を使うか, 上限を使うか]` と その値 |
# | `tied` | **他のパラメータに縛る式**。例 `p[2]`（第 1 成分と同じ幅）、`p[1]+0.06` |

# %%
names = eispac.data.fit_template_filenames()
print(f"同梱テンプレート {len(names)} 個。Fe XII のもの:")
for n in sorted(str(x).split("/")[-1] for x in names):
    if n.startswith("fe_12"):
        print("  ", n)
print("\n観測に合うものを探すには eispac.match_templates(eis_obs)")

# %% [markdown]
# ### 5. フィットする
#
# ```python
# fit = eispac.fit_spectra(cube, tmplt, ncpu=2, ignore_warnings=True)   # 全画素
# fit = eispac.fit_spectra(cube[180:340, :, :], tmplt, ncpu=2)          # 一部だけ
# fit = eispac.fit_spectra(inten, tmplt, wave=wave, errs=sig, ncpu=1)   # 1 本だけ
# ```
#
# | 引数 | 既定値 | 意味 |
# |---|---|---|
# | `inten` | — | `EISCube`、配列、ファイルパス。**3 次元ならラスター、1 次元なら 1 本**として扱う |
# | `template` | — | `EISFitTemplate`、辞書、テンプレートファイルのパス |
# | `parinfo` | `None` | テンプレートの設定を上書きする（自作テンプレートを作るとき。付録 G） |
# | `wave` | `None` | 波長配列。**`inten` を配列で渡すときは必須** |
# | `errs` | `None` | 誤差配列。同上 |
# | `min_points` | `7` | この数より有効点が少ないスペクトルは飛ばす。**パラメータ数以上が必要** |
# | `ncpu` | `'max'` | 並列数。`'max'` / `None` で全コア |
# | `ignore_warnings` | `False` | `True` にすると警告を出さない（大量に出るので実用上は `True`） |
# | `skip_fitting` | `False` | 初期値のまま返す（設定の確認用） |
#
# **★ `ncpu > 1` の注意**: スクリプトから使うときは
# `if __name__ == "__main__":` で囲まないと、単一プロセスに落ちます
# （ノートブックでは気にしなくて大丈夫です）。
#
# 結果 `EISFitResult` の中身（`fit.fit[...]`）:
#
# | キー | 中身 |
# |---|---|
# | `int` / `err_int` | **輻射強度**とその誤差（成分ごと） |
# | `params` / `perror` | 生のパラメータ。`[振幅, 中心, 幅] × 成分数 + 背景` |
# | `chi2` | フィットの χ² |
# | `line_ids` | 成分の並び |
# | `status` | mpfit の終了状態（負なら失敗） |

# %% [markdown]
# ### 6. 結果を取り出す
#
# ```python
# m = fit.get_map(component=0, measurement="intensity")
# ```
#
# | 引数 | 既定値 | 意味 |
# |---|---|---|
# | `component` | `0` | **何番目のガウシアンか**。`line_ids` を見て決める |
# | `measurement` | `'intensity'` | `'intensity'` / `'velocity'` / `'width'` |
#
# 返り値は `EISMap`（sunpy の `Map`）なので、`.data`、`.wcs`、`.plot()` が使えます。
#
# ```python
# wave, prof = fit.get_fit_profile()          # モデル曲線（観測点と同じ波長）
# wave, prof = fit.get_fit_profile(num_wavelengths=100)   # 滑らかに描きたいとき
# ```

# %% [markdown]
# ### 7. 速度に直す
#
# ```python
# from eispac.instr import calc_velocity
# v = calc_velocity(fit.fit["params"][..., 1], 195.119, corr_method="column")
# ```
#
# | 引数 | 既定値 | 意味 |
# |---|---|---|
# | `observed_wave` | — | 観測された中心波長の配列 |
# | `rest_wave` | — | 静止波長。`"Fe XII 195.119"` のような文字列でも可 |
# | `corr_method` | `'column'` | **ゼロ点の決め方**（第 3 章）。`'column'` = 列ごとの中央値を 0、`'image'` = 視野全体の中央値を 0、`None` = 何もしない |

# %% [markdown]
# ### 8. 保存する
#
# ```python
# eispac.save_fit(fit, save_dir="work")     # HDF5 で保存
# fit = eispac.read_fit("work/..._fit.h5")  # 読み戻す
# eispac.export_fits(fit, save_dir="work")  # FITS で書き出す
# ```
#
# `save_dir` を省くと元データと同じ場所に置きます。
# 全ラスターのフィットは数分かかるので、**結果は保存しておくのが実用的**です。

# %% [markdown]
# ### つまずきやすいところ
#
# | | |
# |---|---|
# | 成分の取り違え | `line_ids` を見ずに `component=0` とすると別のイオンが返る（第 2 章） |
# | 欠損値 | NaN ではなく大きな負の値。`cube.mask` で落とす（付録 A） |
# | 単位 | `cube.data` は Å⁻¹ あたり。波長刻みを掛けて輻射強度（付録 D） |
# | 波長のゼロ点 | 絶対基準が無い。自分で決める（第 3 章） |
# | 実行時間 | 全ラスター 1 輝線で数分。`ncpu` を上げるか範囲を絞る |
# | 警告の量 | `ignore_warnings=True` を付けないと画素ごとに出て読めなくなる |

# %% [markdown]
# ---
#
# ## この教材について
#
# - リポジトリ: https://github.com/hottahd/EIS_practice
# - ライセンス: CC BY 4.0（出典を示せば自由に使えます）
# - 検証済みの解析スクリプトは `scripts/`、準備の記録は `docs/`
#
# 元になった論文:
# Warren, H. P., Winebarger, A. R., & Brooks, D. H. 2012,
# *"A Systematic Survey of High-Temperature Emission in Solar Active Regions"*,
# ApJ, 759, 141 — https://doi.org/10.1088/0004-637X/759/2/141
