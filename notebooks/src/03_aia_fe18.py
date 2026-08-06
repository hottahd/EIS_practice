# %% [markdown]
# # モジュール 3: AIA 94 Å から Fe XVIII (7 MK) を取り出す
#
# **所要時間 50 分**
#
# **このノートで身につくこと**
#
# 1. **なぜ EIS だけでは足りないか**（EIS の最高温は Ca XVII の ~5 MK）
# 2. AIA 94 Å の Fe XVIII を経験式で分離する（論文 Appendix）
# 3. ★ **論文に印刷された式の指数が誤植である**ことを、実データで自分で確かめる
# 4. moss（171 Å）と高温ループ（Fe XVIII）の空間分布の違いを見る
#    → 次のモジュールで箱を選ぶための下準備
#
# 前提: モジュール 1, 2。EIS のデータは使わないので単独でも動く。

# %%
!pip install -q sunpy

# %% [markdown]
# ### Colab のためのおまじない（ローカルで動かしている人は素通りします）
#
# **Colab はノートブック 1 冊ごとに新しい仮想マシンが立ち上がる。**
# GitHub から開いただけでは教材リポジトリもデータも無いので、ここで用意する。
# セッションが切れたときも、このセルをもう一度実行すれば復帰できる。

# %%
import os
import subprocess
import sys

REPO = "https://github.com/hottahd/EIS_practice.git"
if not os.path.exists("scripts/lines_warren2012.py"):      # リポジトリの外にいる
    if not os.path.exists("EIS_practice"):
        print("教材リポジトリを取得中 ...")
        subprocess.run(["git", "clone", "-q", REPO], check=True)
    os.chdir("EIS_practice")
sys.path.insert(0, "scripts")
print("作業ディレクトリ:", os.getcwd())

# %% [markdown]
# ## 3-0. データを取る（登録不要）
#
# **JSOC の synoptic アーカイブ**を使う。1024×1024（2.4″/画素）、1 枚 1 MB 弱。
#
#     http://jsoc.stanford.edu/data/aia/synoptic/YYYY/MM/DD/HHHH/AIAyyyymmdd_hhmm_wwww.fits
#
# - フルディスク level-1（4096²、1 枚 65 MB）は VSO 経由で取れるが遅い。
#   **Colab では synoptic を使う。**
# - synoptic は既に **level-1.5**（`lvl_num=1.5`）なので `aiapy` の
#   `register` / `update_pointing` は不要。読んですぐ使える。
# - 時刻 **03:38 UT** は EIS ラスター（03:07–04:09）のほぼ中央。

# %%
import os
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
from astropy.coordinates import SkyCoord

AIA_DIR = "data/sdo/synoptic"
AIA_BASE = "http://jsoc.stanford.edu/data/aia/synoptic/2011/07/02/H0300"


def ensure(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print("downloading", os.path.basename(path))
    urllib.request.urlretrieve(url, path)
    return path


def load_aia(wave):
    """synoptic の AIA を読んで DN/s にする。"""
    f = ensure(f"{AIA_BASE}/AIA20110702_0338_{wave:04d}.fits",
               f"{AIA_DIR}/AIA20110702_0338_{wave:04d}.fits")
    m = sunpy.map.Map(f)
    exp = m.meta["exptime"]
    print(f"AIA {wave:4d}: {m.data.shape}  exptime={exp:.3f}s  "
          f"lvl={m.meta.get('lvl_num')}  {m.scale[0]:.2f}")
    return sunpy.map.Map(m.data / exp, m.meta)      # ★ 露光時間で割る


a94, a171, a193 = load_aia(94), load_aia(171), load_aia(193)

# %% [markdown]
# ## 3-1. なぜ AIA が要るのか
#
# EIS で観測できる**最高温の強い輝線は Ca XVII 192.858（~5 MK）**。
# それより上を拘束するものが無いと、EM 分布の高温側の裾が決まらない。
#
# しかも Ca XVII はブレンドしていて扱いが難しい（モジュール 8）。
# **7 MK に効く独立な測定**がほしい。
#
# そこで **AIA 94 Å の Fe XVIII 93.932 Å**（形成温度 log T ≈ 6.85 = 7 MK）を使う。
#
# **問題**: 94 Å チャンネルには Fe XVIII 以外に
# **Fe X 94.012 Å をはじめとする低温の線**が入っている。
# 静穏領域やループの足元では、94 Å の信号の**ほとんどが低温成分**。
# そのまま使うと「7 MK のプラズマが大量にある」ことになってしまう。

# %% [markdown]
# ## 3-2. 論文 Appendix の経験式
#
# 論文は 171 Å と 193 Å の混合から「warm 成分」を経験的に見積もって引く:
#
# $$ x = \frac{f\,I_{171} + (1-f)\,I_{193}}{116.54}, \qquad f = 0.31 $$
# $$ I_{94}^{\rm warm} = 0.39 \sum_i a_i x^{?}, \qquad
#    a = [-7.31\times10^{-2},\ 9.75\times10^{-1},\ 9.90\times10^{-2},\ -2.84\times10^{-3}] $$
# $$ I_{\rm FeXVIII} = I_{94} - I_{94}^{\rm warm} $$
#
# 較正の元データは 2010-03-22 12–13 UT の 1 時間平均（明るい輝点＋静穏領域）で、
# 規格化定数 116.54 と 0.39 は**そのときの中央値**。
# つまり「**中央値の場所では x = 1、warm 成分 = 0.39 DN/s**」という約束になっている。
#
# **★ 問題は指数 `?` である。** 論文の印刷は
#
# $$ I_{94}^{\rm warm} = 0.39 \sum_{i=1}^{4} a_i x^{i} $$
#
# だが、この字面どおりだと定数項が無い（x=0 で 0）。
# 一方、定数項つきの 3 次式（$a_i x^{i-1}$）という読み方もできる。
# **どちらが正しいかは、実データを通してみれば決まる。**

# %% [markdown]
# ## 3-3. ★ 誤植を自分で確かめる
#
# 2 つの読み方それぞれで warm 成分を計算し、**観測された 94 Å の強度**と比べる。
# warm 成分は 94 Å の**一部**なのだから、観測値を超えたらその読み方は誤り。

# %%
A = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]
F_MIX, NORM_COMPOSITE, NORM_94, X_MAX = 0.31, 116.54, 0.39, 30.0

x = (F_MIX * a171.data + (1 - F_MIX) * a193.data) / NORM_COMPOSITE
x = np.clip(x, 0.0, X_MAX)          # 論文の指示どおり x は 30 で頭打ち

warm_printed = NORM_94 * sum(a * x**(i + 1) for i, a in enumerate(A))   # a_i x^i
warm_fixed = NORM_94 * sum(a * x**i for i, a in enumerate(A))           # a_i x^(i-1)

print(f"{'x の範囲':>12} {'観測 I_94':>10} {'印刷どおり':>12} {'定数項つき':>12}  {'画素数':>8}")
for lo, hi in [(0.5, 1.5), (2, 3), (5, 8), (12, 18), (25, 30)]:
    sel = (x >= lo) & (x < hi)
    if sel.sum() == 0:
        continue
    print(f"{lo:5.1f} - {hi:4.1f} {np.median(a94.data[sel]):10.2f} "
          f"{np.median(warm_printed[sel]):12.2f} {np.median(warm_fixed[sel]):12.2f}"
          f" {sel.sum():9d}")

# %% [markdown]
# **結論は一目瞭然**: 印刷どおりの読み方だと、warm 成分が
# **観測された 94 Å の全強度を桁違いに超える**。物理的にありえない。
# 定数項つきの 3 次式なら観測値の少し下に収まる。
#
# → **論文 Eq.(A1) の指数は組版上の誤りで、正しくは $a_i x^{i-1}$**（定数項つき）。
#
# **★ 注意**: x = 1（較正の中央値）では**どちらの読み方も 0.389 になる**。
# 論文が明記している唯一の数値がここなので、**この点だけでは判別できない**。
# 実データを広い x の範囲で通して初めて決まる。
#
# **教訓**: 「論文に書いてある式」でも、**数値を通すまで信じない**。
# 特に指数・添字・単位は組版で落ちやすい。

# %%
# 差し引きの引き算そのもの
def fe18(i94, i171, i193, f=F_MIX):
    xx = np.clip((f * np.asarray(i171, float) + (1 - f) * np.asarray(i193, float))
                 / NORM_COMPOSITE, 0.0, X_MAX)
    warm = NORM_94 * sum(a * xx**i for i, a in enumerate(A))
    return np.asarray(i94, float) - warm


# 多項式が x=30 で頭打ちにされている理由も確かめておく
xs = np.linspace(0, 40, 4001)
poly = sum(a * xs**i for i, a in enumerate(A))
print(f"多項式の極大は x = {xs[np.argmax(poly)]:.1f}")
print("→ 論文が x を 30 で頭打ちにしているのは、この折り返しの手前で止めるため")
print(f"x=1 での warm = {NORM_94*sum(a*1.0**i for i, a in enumerate(A)):.3f} DN/s "
      f"（論文の中央値 0.39 と一致）")

# %% [markdown]
# ## 3-4. Fe XVIII マップを作る
#
# 活動領域の周りだけ切り出して 4 枚並べる。

# %%
XC, YC = -330, 200          # NOAA 1243 のだいたいの位置 [arcsec]
HALF = 200 * u.arcsec

bl = SkyCoord((XC * u.arcsec) - HALF, (YC * u.arcsec) - HALF, frame=a94.coordinate_frame)
tr = SkyCoord((XC * u.arcsec) + HALF, (YC * u.arcsec) + HALF, frame=a94.coordinate_frame)

# ★ 罠: 同じ SkyCoord で切り出しても、波長ごとに配列サイズが 1 画素ずれることがある
for m, w in [(a94, 94), (a171, 171), (a193, 193)]:
    print(f"AIA {w:4d} を SkyCoord で切り出すと {m.submap(bl, top_right=tr).data.shape}")

# %% [markdown]
# **サイズが揃わない。** 引き算しようとすると
# `operands could not be broadcast together` で落ちる。
#
# 原因は観測者距離 `dsun` がチャンネルごとに僅かに違い、
# arcsec → 画素の換算が 1 画素未満ずれて、切り上げ・切り捨ての境目に乗ること。
#
# **対処は 2 つ。今回のデータではどちらでもよいが、意味が違う。**
#
# 1. **同じ画素範囲で切る**（今回の synoptic は 3 波長とも
#    `crpix/crval/cdelt/crota` が完全に一致しているので、これで厳密に正しい）
# 2. **`reproject_to()` で揃える**（WCS が違う場合はこちらが必須。
#    フルディスク level-1 を自分で `aiapy.calibrate.register` した場合など)
#
# **★ 確認せずに引き算しないこと。** 揃っていない画像を引くと、
# 構造の縁に沿って偽の正負のパターンが出る。moss の縁が全部 Fe XVIII に化ける。

# %%
# 3 波長の WCS が同一であることを確認してから、同じ画素範囲で切る
for k in ("crpix1", "crpix2", "crval1", "crval2", "cdelt1", "crota2"):
    vals = {m.meta[k] for m in (a94, a171, a193)}
    assert len(vals) == 1, f"{k} が一致しない: {vals}"
print("3 波長の WCS は完全に一致 → 同じ画素範囲で切ってよい")

sub94 = a94.submap(bl, top_right=tr)
ny, nx = sub94.data.shape
px = a94.world_to_pixel(bl)
i0, j0 = int(round(px.y.value)), int(round(px.x.value))


def cut(m):
    return sunpy.map.Map(m.data[i0:i0 + ny, j0:j0 + nx], sub94.meta)


s94, s171, s193 = cut(a94), cut(a171), cut(a193)
s18 = sunpy.map.Map(fe18(s94.data, s171.data, s193.data), sub94.meta)
print("揃えた後:", s94.data.shape, s171.data.shape, s193.data.shape)

# %%
panels = [("AIA 171 (moss, 0.9 MK)", s171, "sdoaia171"),
          ("AIA 193 (1.6 MK)", s193, "sdoaia193"),
          ("AIA 94 (raw)", s94, "sdoaia94"),
          ("AIA 94 -> Fe XVIII (7 MK)", s18, "inferno")]

fig = plt.figure(figsize=(16, 4.6))
for k, (title, m, cmap) in enumerate(panels):
    ax = fig.add_subplot(1, 4, k + 1, projection=m)
    d = np.sqrt(np.clip(m.data, 0, None))
    ax.imshow(d, origin="lower", cmap=cmap,
              vmin=np.nanpercentile(d, 1), vmax=np.nanpercentile(d, 99.7))
    ax.set_title(title, fontsize=10)
    ax.grid(False)
    ax.set_xlabel("Solar X")
    ax.set_ylabel("Solar Y" if k == 0 else "")
    if k > 0:
        ax.coords[1].set_ticklabel_visible(False)
fig.tight_layout()
plt.show()

# %% [markdown]
# **見るべきこと**
#
# - **171 Å**: 活動領域コアの中に**まだらの明るい模様** = moss。
#   高温ループの足元が遷移層で光っている。
# - **94 Å（生）**: 171 の moss がそのまま透けて見える。**これが汚染**。
# - **Fe XVIII**: moss が消え、**滑らかで太いループ**だけが残る。
#   これが 7 MK のプラズマ。
#
# → **Fe XVIII で明るく、171 の moss を含まない場所**が、
#   論文の言う inter-moss 領域。次のモジュールで選ぶ。

# %% [markdown]
# ## 3-5. 論文 Table 2 の最終行と答え合わせ
#
# 論文 Table 2 の最終行は EIS の輝線ではなく
# **AIA 94 Å（Fe XVIII 分離後）の 7.20 ± 1.40 DN/s** で、
# EIS と同じ inter-moss 箱で測った値。
#
# **EIS とは完全に独立な測定**なので、箱の位置の検証に使える。
# 論文は箱の座標を書いていないが、Figure 2 の緑枠から実測できる
# （`scripts/extract_paper_boxes.py`）。

# %%
PAPER_BOX = dict(x0=-321.8, x1=-306.4, y0=202.6, y1=226.2)     # Fig.2 から実測
m18 = sunpy.map.Map(fe18(a94.data, a171.data, a193.data), a94.meta)

bl = SkyCoord(PAPER_BOX["x0"] * u.arcsec, PAPER_BOX["y0"] * u.arcsec,
              frame=m18.coordinate_frame)
tr = SkyCoord(PAPER_BOX["x1"] * u.arcsec, PAPER_BOX["y1"] * u.arcsec,
              frame=m18.coordinate_frame)
sub = m18.submap(bl, top_right=tr)
v = float(np.nanmean(sub.data))
print(f"論文の箱 X=[{PAPER_BOX['x0']}, {PAPER_BOX['x1']}] "
      f"Y=[{PAPER_BOX['y0']}, {PAPER_BOX['y1']}]  ({sub.data.shape[0]}x{sub.data.shape[1]} px)")
print(f"  Fe XVIII 平均 = {v:.3f} DN/s")
print(f"  論文 Table 2  = 7.200 ± 1.40 DN/s")
print(f"  ratio         = {v/7.20:.3f}")

# %% [markdown]
# **0.86 前後**になったはず。論文の誤差 ±19% の範囲に入っている。
#
# **★ ただし、測り方で 10% 動く**。同じ箱でも
#
# | 使ったデータ | Fe XVIII | 対論文 |
# |---|---:|---:|
# | synoptic 1024² (2.4″/px) をそのまま | 6.17 | 0.86 |
# | フルディスク level-1 (4096², 0.6″/px) を EIS 格子に再投影 | 6.81 | 0.95 |
#
# 15″×23″ の箱は 2.4″/画素だと **6×10 画素しかない**。
# 構造のある場所を粗い格子で測ると、平均は端の扱いで簡単に 10% 動く。
#
# **教訓**: 「小さい箱を粗い画像で測る」ときは、
# **解像度そのものが系統誤差になる**。論文と比べる前に、
# 自分の測定がどちらの向きに偏りうるかを見積もっておく。

# %% [markdown]
# ## 3-6. この経験式の適用限界（必ず知っておくこと）
#
# 1. **フレア中は使えない。** Fe XXIV 192.04 Å が 193 Å チャンネルに入るので、
#    warm 成分を過大評価して Fe XVIII が負になる。
# 2. **非常に明るい moss でも破綻する。** 較正データに無い強度域。
# 3. **AIA の感度劣化補正を掛けてはいけない。**
#    この経験式は劣化補正**なし**のデータ（2010 年 3 月）で導かれている。
#    `aiapy.calibrate.degradation()` を掛けると係数と整合しなくなる。
#    - ただし 2010 年の較正を 2011 年以降のデータに使うことの是非は
#      それ自体が議論の種。**論文どおりの再現をするなら掛けない。**
# 4. より丁寧にやるなら、**Fe XVIII だけの応答関数**を作って DEM に入れる
#    （公式の 94 Å 応答は低温線込みなので使えない）。
#    → `scripts/aia94_fe18_response.py`。モジュール 7 で使う。

# %% [markdown]
# ## 3-7. 演習
#
# 1. `f = 0.31` を 0.2 や 0.5 に変えて Fe XVIII マップがどう変わるか見る。
#    moss の消え方が変わるはず。**どの f が「正しい」と言えるか？**
# 2. Fe XVIII が**負になる画素**を数えて、どこに分布するか描いてみる。
#    （ヒント: 静穏領域と、非常に明るい moss）
# 3. 論文の経験式を**自分で導き直す**（発展）。
#    静穏領域を含む広い領域で `x` と `I_94` の散布図を作り、3 次式を当てる。
#    論文 Appendix の追体験になる。
#
# ## まとめ
#
# - EIS の上限（~5 MK）より上を拘束するために **AIA 94 Å の Fe XVIII** を使う
# - 94 Å は低温線に汚染されている → **171/193 から経験的に引く**
# - **論文 Eq.(A1) の指数は誤植**。実データを通せば自分で確かめられる
# - Fe XVIII で明るく 171 の moss が無い場所 = **inter-moss 領域**
# - 小さい箱を粗い画像で測ると、**解像度が 10% の系統誤差**になる
#
# 次（モジュール 4）では、EIS と AIA の座標を合わせて
# **実際に箱を選ぶ**。
