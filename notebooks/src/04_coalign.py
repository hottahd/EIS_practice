# %% [markdown]
# # モジュール 4: EIS と AIA の座標合わせ、inter-moss 領域の選択
#
# **所要時間 40 分**
#
# **このノートで身につくこと**
#
# 1. **EIS のポインティングは信用できない**ことを自分の目で確かめ、
#    相互相関で合わせる
# 2. AIA を EIS の画素格子に載せ替える（`reproject_to`）
# 3. **inter-moss 領域を選ぶ**。ここが解析で最も主観的な部分
# 4. 「どこを測るか」が結果を左右することを、次のモジュールへの伏線として持つ
#
# 前提: モジュール 1, 3。

# %%
!pip install -q eispac sunpy

# %%
import os
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import sunpy.map
import eispac
from astropy.coordinates import SkyCoord

EIS_FILE = "data/eis/eis_20110702_030712.data.h5"
AIA_DIR = "data/sdo/synoptic"
AIA_BASE = "http://jsoc.stanford.edu/data/aia/synoptic/2011/07/02/H0300"


def ensure(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print("downloading", os.path.basename(path))
    urllib.request.urlretrieve(url, path)
    return path


for ext in ("data", "head"):
    ensure(f"https://eis.nrl.navy.mil/level1/hdf5/2011/07/02/"
           f"eis_20110702_030712.{ext}.h5",
           f"data/eis/eis_20110702_030712.{ext}.h5")


def load_aia(wave):
    f = ensure(f"{AIA_BASE}/AIA20110702_0338_{wave:04d}.fits",
               f"{AIA_DIR}/AIA20110702_0338_{wave:04d}.fits")
    m = sunpy.map.Map(f)
    return sunpy.map.Map(m.data / m.meta["exptime"], m.meta)   # DN/s


a94, a171, a193 = load_aia(94), load_aia(171), load_aia(193)
print("AIA ok")

# %% [markdown]
# ## 4-1. EIS のラスターを sunpy Map にする
#
# 座標を扱うには WCS 付きの Map がほしい。
#
# **フィットは要らない**。相互相関に使うだけなので、
# モジュール 1 のクイックルック（波長方向の積分）で十分。
# 全ラスターをフィットすると 3 分かかるが、これなら 1 秒。

# %%
def eis_raster_map(datafile, wvl):
    """EIS のスペクトルウィンドウを波長方向に積んで sunpy Map にする。"""
    c = eispac.read_cube(datafile, wvl)
    d = np.where(np.asarray(c.mask, dtype=bool), np.nan, c.data)   # 欠損を除く
    img = np.nanmean(d, axis=2) * d.shape[2]

    h, p = c.meta["index"], c.meta["pointing"]
    ref = SkyCoord(p["xcen"] * u.arcsec, p["ycen"] * u.arcsec,
                   obstime=h["date_obs"], observer="earth", frame="helioprojective")
    hdr = sunpy.map.make_fitswcs_header(
        img, ref, scale=[p["x_scale"], p["y_scale"]] * u.arcsec / u.pix,
        instrument="EIS", wavelength=wvl * u.angstrom)
    hdr["measrmnt"] = "intensity"          # eispac の EISMap が要求するキー
    return sunpy.map.Map(img, hdr)


m_eis = eis_raster_map(EIS_FILE, 195.119)
print("EIS Fe XII map:", m_eis.data.shape)
print(f"  X = {m_eis.bottom_left_coord.Tx.value:.1f} .. "
      f"{m_eis.top_right_coord.Tx.value:.1f}\"")
print(f"  Y = {m_eis.bottom_left_coord.Ty.value:.1f} .. "
      f"{m_eis.top_right_coord.Ty.value:.1f}\"")
print(f"  画素 = {m_eis.scale[0]:.2f} x {m_eis.scale[1]:.2f}")

# %% [markdown]
# ## 4-2. AIA を EIS の格子に載せ替える
#
# `reproject_to` は WCS を見て座標変換してくれる。
# EIS の格子は x 2″ × y 1″ という**変な格子**だが、気にせず載せられる。

# %%
r193, r171, r94 = (m.reproject_to(m_eis.wcs) for m in (a193, a171, a94))
print("AIA on EIS grid:", r193.data.shape)

A = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]


def fe18(i94, i171, i193, f=0.31):
    x = np.clip((f * i171 + (1 - f) * i193) / 116.54, 0.0, 30.0)
    return i94 - 0.39 * sum(a * x**i for i, a in enumerate(A))


fe = fe18(r94.data, r171.data, r193.data)

# %% [markdown]
# ## 4-3. ★ ずれているのを見る
#
# EIS Fe XII 195.119（1.6 MK）と AIA 193（1.6 MK）は**ほぼ同じ温度**なので、
# 形態がよく似ているはず。並べて見る。

# %%
fig, axes = plt.subplots(1, 2, figsize=(7, 9))
for ax, (d, t) in zip(axes, [(m_eis.data, "EIS Fe XII 195.1"),
                             (r193.data, "AIA 193 (on EIS grid)")]):
    v = np.sqrt(np.clip(d, 0, None))
    lo, hi = np.nanpercentile(v, [1, 99.5])
    ax.imshow(v, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi)
    ax.set_title(t, fontsize=10)
    ax.set_xlabel("EIS x [pix]")
axes[0].set_ylabel("EIS y [pix]")
fig.tight_layout()
plt.show()

# %% [markdown]
# よく似ているが、**同じ場所に重なっていない**。
#
# **なぜずれるか**
#
# - EIS のポインティングには **数″〜十数″の系統誤差**がある
#   （熱変形、姿勢基準の違い、軌道中の指向ドリフト）
# - ヘッダの `xcen`/`ycen` を鵜呑みにすると、選んだ箱が別の場所を指す
# - 論文は箱を「AIA Fe XVIII で明るく AIA 171 の moss が無い場所」として
#   選んでいるので、**この座標合わせの精度がそのまま強度の精度になる**

# %% [markdown]
# ## 4-4. 相互相関でずれを測る
#
# 整数画素精度で十分（EIS の画素は 1″ × 2″）。
#
# **★ 相関係数の取り方に注意**: 重なり領域**ごとに**正規化した Pearson 相関を使う。
# 配列全体で一度だけ正規化すると、重なりが小さいシフトほど見かけの相関が
# 上がってしまい、ずれを過大評価する（実際にこれで一度間違えた）。

# %%
def cross_correlate_shift(ref, img, max_shift=25):
    """img を ref に合わせるための (dy, dx) を返す。"""
    def prep(a):
        a = np.array(a, float)
        a[~np.isfinite(a)] = np.nanmedian(a)
        return np.sqrt(np.clip(a, 0, None))     # 明るいコアだけで決まらないように

    r, i = prep(ref), prep(img)
    ny, nx = r.shape
    best, bdy, bdx = -np.inf, 0, 0
    for dy in range(-max_shift, max_shift + 1):
        for dx in range(-max_shift, max_shift + 1):
            rs = r[max(0, dy):ny + min(0, dy), max(0, dx):nx + min(0, dx)]
            is_ = i[max(0, -dy):ny + min(0, -dy), max(0, -dx):nx + min(0, -dx)]
            if rs.size < 0.5 * r.size:
                continue
            rc, ic = rs - rs.mean(), is_ - is_.mean()
            den = np.sqrt((rc**2).sum() * (ic**2).sum())
            c = float((rc * ic).sum() / den) if den > 0 else -np.inf
            if c > best:
                best, bdy, bdx = c, dy, dx
    return bdy, bdx, best


dy, dx, cc = cross_correlate_shift(m_eis.data, r193.data)
print(f"ずれ: dy = {dy} pix (= {dy*1.0:.0f}\"),  dx = {dx} pix "
      f"(= {dx*2.0:.0f}\")   相関 {cc:.3f}")

# %%
def shift(a, dy, dx):
    out = np.full_like(np.asarray(a, float), np.nan)
    ny, nx = a.shape
    out[max(0, dy):ny + min(0, dy), max(0, dx):nx + min(0, dx)] = \
        a[max(0, -dy):ny + min(0, -dy), max(0, -dx):nx + min(0, -dx)]
    return out


a193s, a171s, fes = (shift(v, dy, dx) for v in (r193.data, r171.data, fe))

fig, axes = plt.subplots(1, 4, figsize=(13, 9))
for ax, (d, t) in zip(axes, [(m_eis.data, "EIS Fe XII 195.1"),
                             (a193s, "AIA 193 (shifted)"),
                             (a171s, "AIA 171 = moss"),
                             (fes, "AIA Fe XVIII = 7 MK")]):
    v = np.sqrt(np.clip(d, 0, None))
    lo, hi = np.nanpercentile(v, [1, 99.5])
    ax.imshow(v, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi)
    ax.set_title(t, fontsize=10)
    ax.set_xlabel("EIS x [pix]")
axes[0].set_ylabel("EIS y [pix]")
fig.suptitle(f"co-aligned  (dy={dy}, dx={dx} pix)")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

# %% [markdown]
# ## 4-5. inter-moss 領域を選ぶ
#
# 論文の言い方:
#
# > the inter-moss region, that is, the region between the loop footpoints where
# > we are measuring the properties near the loop apex
#
# **なぜ moss を避けるのか**
#
# - moss = 高温ループの**足元**が遷移層で光っているもの（171 Å で明るい）
# - moss を含むと、視線上に「足元の 1 MK」と「ループ上部の 4 MK」が混ざる
# - 足元は熱伝導・彩層蒸発が絡む複雑な物理。**ループ上部だけを見たい**
#
# **数値化**: `score = median(Fe XVIII) / median(AIA 171)`
# （高温で明るく、moss が暗いほど大きい）
#
# **★ ただし最後は図を見て人間が決める。** 自動化しきらないのが正しい。
# 物理的な判断だからである。

# %%
def scan(fe, a171, ny, nx, step=4, min_fe=3.0):
    H, W = fe.shape
    out = []
    for y0 in range(0, H - ny, step):
        for x0 in range(0, W - nx, step):
            f, m = fe[y0:y0+ny, x0:x0+nx], a171[y0:y0+ny, x0:x0+nx]
            if not (np.isfinite(f).all() and np.isfinite(m).all()):
                continue
            fmed, mmed = np.median(f), np.median(m)
            if fmed < min_fe:            # Fe XVIII が暗い場所は候補にしない
                continue
            out.append((fmed / mmed, fmed, mmed, y0, y0+ny, x0, x0+nx))
    return sorted(out, reverse=True)


print(f"{'score':>8} {'FeXVIII':>9} {'AIA171':>8}   box")
cands = scan(fes, a171s, 30, 8)
for c in cands[:8]:
    s, f, m, y0, y1, x0, x1 = c
    print(f"{s:8.4f} {f:9.2f} {m:8.0f}   y=[{y0}:{y1}] x=[{x0}:{x1}]")

# %% [markdown]
# **★ 1 位に出てきた `y=[244:274] x=[32:40]` は、我々が採用した箱そのもの。**
#
# この箱は準備段階で「論文 Table 2 に最もよく合う場所」として
# 22 輝線の突き合わせから選んだものだが、
# **論文が書いている選択基準（Fe XVIII で明るく moss が無い）だけを
# 機械的に適用しても同じ場所に来る**。
#
# 選択基準が言葉どおりに再現できている、という確認になる。
# 逆に言えば、**この一致が無ければ「たまたま合う箱を探した」ことになる**。
# 答えを知っている問題では、ここを分けて考えるのが大事。
#
# ### 採用する箱
#
# 上位候補と、**論文が使った箱**（Figure 2 の緑枠から実測したもの、
# `scripts/extract_paper_boxes.py`）を重ねて見る。
#
# 論文は箱の座標を書いていないので、**図から読むしかない**。
# 論文の箱は 15.4″ × 23.5″。他の論文（Winebarger et al. 2011）も
# 5″ × 25″ と**縦長の細い箱**を使っており、それがこのグループの流儀。

# %%
BOX = dict(y0=244, y1=274, x0=32, x1=40)       # 準備段階で採用した箱

fig, axes = plt.subplots(1, 3, figsize=(11, 9))
for ax, (d, t) in zip(axes, [(fes, "AIA Fe XVIII"),
                             (a171s, "AIA 171 (moss)"),
                             (m_eis.data, "EIS Fe XII 195.1")]):
    v = np.sqrt(np.clip(d, 0, None))
    lo, hi = np.nanpercentile(v, [1, 99.5])
    ax.imshow(v, origin="lower", aspect="auto", cmap="inferno", vmin=lo, vmax=hi)
    for c in cands[:8]:                        # 候補（水色）
        _, _, _, y0, y1, x0, x1 = c
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0],
                color="cyan", lw=0.7, alpha=0.5)
    b = BOX                                    # 採用（白）
    ax.plot([b["x0"], b["x1"], b["x1"], b["x0"], b["x0"]],
            [b["y0"], b["y0"], b["y1"], b["y1"], b["y0"]], color="white", lw=2)
    ax.set_title(t, fontsize=10)
    ax.set_xlabel("EIS x [pix]")
axes[0].set_ylabel("EIS y [pix]")
fig.suptitle("inter-moss candidates (cyan) and adopted box (white)")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

# %%
b = BOX
print(f"採用箱 y=[{b['y0']}:{b['y1']}] x=[{b['x0']}:{b['x1']}]"
      f"  = {b['y1']-b['y0']}\" x {(b['x1']-b['x0'])*2}\"")
v18 = float(np.nanmean(fes[b["y0"]:b["y1"], b["x0"]:b["x1"]]))
print(f"  AIA Fe XVIII 平均 = {v18:6.2f} DN/s"
      f"   （論文 Table 2 = 7.20 ± 1.40 → 比 {v18/7.20:.2f}）")
print(f"  AIA 171     平均 = {np.nanmean(a171s[b['y0']:b['y1'], b['x0']:b['x1']]):6.0f} DN/s"
      f"   （視野の中央値 {np.nanmedian(a171s):.0f}）")

os.makedirs("data/cache", exist_ok=True)
np.savez("data/cache/aia_on_eis_grid.npz", aia171=a171s, aia193=a193s, fe18=fes,
         eis_fe12=m_eis.data, dy=dy, dx=dx)
print("\nwrote data/cache/aia_on_eis_grid.npz  （モジュール 5, 7 で使う）")

# %% [markdown]
# ## 4-6. ★ ここが解析で一番主観的なところ
#
# 論文は**目で見て手で箱を選んでおり、座標を書いていない**。
# 我々が Figure から復元した箱と、指標で選んだ箱は近いが同じではない。
#
# 準備段階で 230 通りの箱を総当たりした結果:
#
# | 箱の選び方 | 論文 Table 2 との median 比 |
# |---|---|
# | 適当に明るいところ | 0.3 〜 1.2 まで散らばる |
# | inter-moss の条件を満たす箱 | 0.83 〜 0.95 |
#
# **箱を動かすだけで結果は 2 割動く。** これは論文の誤差 ±22% と同程度。
#
# → **「どこを測るかを決めるのが解析の本体」**である。
#   装置やコードの議論より、まずここを疑う。
#
# 次のモジュールでは、この箱で出した 22 輝線の強度を論文と突き合わせ、
# **箱の選び方が正しいかを数値で診断する**方法を学ぶ。

# %% [markdown]
# ## 4-7. 演習
#
# 1. 相互相関の相手を **AIA 171** に変えるとどうなるか。
#    Fe XII (1.6 MK) と 171 (0.9 MK) では形態が違うので、
#    ずれの推定が悪化するはず。**温度の近い組を選ぶ**理由を確かめる。
# 2. `max_shift` を 5 にすると答えが変わるか。境界に張り付いていないか確認する。
# 3. 上位候補の箱をいくつか `BOX` に入れて、
#    モジュール 2 の 22 輝線フィットを回す。ratio がどう動くか。
# 4. `min_fe`（Fe XVIII の下限）を 1.0 に下げると、どんな場所が候補に入ってくるか。
#
# ## まとめ
#
# - **EIS のポインティングはずれている**。AIA との相互相関で合わせる
# - 相手は**温度の近い線**（EIS Fe XII ↔ AIA 193）を選ぶ
# - inter-moss = **Fe XVIII で明るく、171 の moss が無い**場所。ループ上部を見るため
# - **箱の選び方だけで結果は 2 割動く**。ここが解析で最も主観的
#
# 次（モジュール 5）は講習会の山場、**論文 Table 2 との答え合わせ**。
