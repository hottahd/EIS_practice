# %% [markdown]
# # 第 3 章: ドップラー速度を出す
#
# **50 分**
#
# 同じフィットから、中心波長も出ている。これを**ドップラー速度**に直す。
#
# 出てくるのは:
#
# - コロナの物質が**こちらに向かってくるのか、遠ざかるのか**の地図
# - **ドップラー速度のゼロ点**を**自分で決めなければならない**という、分光観測の宿命
#
# Solar-C EUVST は 0.4″・1 秒で同じことをする。**ここで身につけたことは
# そのまま使える**（ゼロ点の決め方は装置が変わっても付いて回る）。

# %% [markdown]
# ## 3-1. ラスターをフィットする
#
# 第 2 章では箱の中で平均してから 1 回フィットした。
# 速度**マップ**がほしいので、今度は**画素ごとにフィット**する。
#
# 全部（512×60）だと数分かかるので、活動領域が写っている
# y = 180–340 の 160 行だけにする（9600 スペクトル、1〜2 分）。

# %%
import numpy as np
import matplotlib.pyplot as plt

from workshop import fit_region

Y0, Y1 = 180, 340
fit, cube = fit_region(wvl=195.119, y0=Y0, y1=Y1, ncpu=2)   # Fe XII 195.119

cen = fit.fit["params"][..., 1]        # 中心波長 [Å]
sig = fit.fit["params"][..., 2]        # 線幅 σ [Å]  ← 第 4 章で使う
inten = fit.fit["int"][..., 0]         # 強度
print("マップの形:", cen.shape)
print(f"中心波長の範囲: {np.nanmin(cen):.4f} – {np.nanmax(cen):.4f} Å")

# %% [markdown]
# ## 3-2. 波長のずれをドップラー速度に直す
#
# 静止波長 $\lambda_0$ からのずれが、視線方向の速度になる:
#
# $$ v = c\,\frac{\lambda_{\rm obs} - \lambda_0}{\lambda_0} $$
#
# Fe XII の静止波長は $\lambda_0 = 195.119$ Å。
# 波長画素 1 つ（0.0223 Å）は 195 Å では **34 km/s** に相当する。
# つまり**波長画素の 1/10 以下のずれ**を測ることになる。

# %%
C_KMS = 2.998e5
LAM0 = 195.119

v_naive = C_KMS * (cen - LAM0) / LAM0

print(f"1 波長画素 = {C_KMS*0.0223/LAM0:.0f} km/s")
print(f"素直に計算したドップラー速度: 中央値 {np.nanmedian(v_naive):+.1f} km/s  "
      f"5–95% [{np.nanpercentile(v_naive,5):+.1f}, {np.nanpercentile(v_naive,95):+.1f}]")

# %% [markdown]
# **視野全体の中央値が 0 になっていない。**
#
# 活動領域全体が数 km/s で一様に動いているわけではないので、これは**装置側のずれ**。
# 理由は次の 2 つで、どちらも分光観測に共通する。

# %% [markdown]
# ### (a) EIS には絶対的な波長基準が無い
#
# 較正用の光源を積んでいないので、「この画素がちょうど 195.119 Å」という
# 基準が無い。**視野の中の何かをゼロと決めるしかない。**
#
# よく使われる決め方:
#
# | 決め方 | 意味 |
# |---|---|
# | 視野全体の中央値を 0 | 「平均的には静止している」と仮定する |
# | **列（露光）ごとの中央値を 0** | 上に加えて、露光ごとの装置の揺らぎも落とす |
# | 静穏領域の値を 0 | 静穏領域が静止していると仮定する |
#
# **どれを選んだかでドップラー速度のマップの意味が変わる。** 論文には必ず書く。
#
# ### (b) 軌道に伴う波長のずれ（eispac が補正済み）
#
# EIS は 98 分で地球を回るあいだに日陰・日照を繰り返し、**装置の温度が変わって
# 波長が動く**。さらにスリットは検出器に対してわずかに傾いている。
#
# eispac は `read_cube` の時点でこれを補正している（`cube.meta['wave_corr']`）。
# **自分で level-0 から処理するなら、この工程は必須。**

# %%
wc = np.asarray(cube.meta["wave_corr"], float)          # (y, x) [Å]
wct = np.asarray(cube.meta["wave_corr_t"], float)       # 時間依存（軌道）
wcx = np.asarray(cube.meta["wave_corr_tilt"], float)    # 位置依存（スリット傾き）

print(f"補正量 全体   : {wc.min():+.4f} 〜 {wc.max():+.4f} Å  "
      f"= {C_KMS*(wc.max()-wc.min())/LAM0:.0f} km/s ぶん")
print(f"  うち 軌道変動: {wct.min():+.4f} 〜 {wct.max():+.4f} Å  "
      f"({C_KMS*(wct.max()-wct.min())/LAM0:.0f} km/s)")
print(f"  うち 傾き    : {wcx.min():+.4f} 〜 {wcx.max():+.4f} Å  "
      f"({C_KMS*(wcx.max()-wcx.min())/LAM0:.0f} km/s)")
print("\n→ コロナの流れ（数十 km/s）と同じ大きさ。補正しなければ速度は測れない。")

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 3.4))
axes[0].plot(wct, lw=1.5)
axes[0].set_xlabel("raster step (= time)")
axes[0].set_ylabel("wavelength shift [Å]")
axes[0].set_title("orbital drift (time-dependent)")
axes[1].plot(wcx, np.arange(len(wcx)), lw=1.5)
axes[1].set_ylabel("y [pix] (along the slit)")
axes[1].set_xlabel("wavelength shift [Å]")
axes[1].set_title("slit tilt (position-dependent)")
for ax in axes:
    ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 3-3. ゼロ点を決めてドップラー速度のマップにする
#
# `eispac.instr.calc_velocity` がやってくれる。
# `corr_method` がゼロ点の決め方に対応する。

# %%
from eispac.instr import calc_velocity

v_image = calc_velocity(cen, LAM0, corr_method="image")    # 視野全体の中央値を 0
v_col = calc_velocity(cen, LAM0, corr_method="column")     # 列ごとの中央値を 0

for name, v in [("ゼロ点なし", v_naive), ("視野の中央値を 0", v_image),
                ("列ごとの中央値を 0", v_col)]:
    print(f"{name:20s} 中央値 {np.nanmedian(v):+6.2f}  "
          f"5–95% [{np.nanpercentile(v,5):+6.2f}, {np.nanpercentile(v,95):+6.2f}] km/s")

# %%
ext = [0, 60, Y0, Y1]
fig, axes = plt.subplots(1, 2, figsize=(9, 7))

d = np.sqrt(np.clip(inten, 0, None))
axes[0].imshow(d, origin="lower", aspect="auto", extent=ext, cmap="inferno",
               vmin=0, vmax=np.nanpercentile(d, 99.5))
axes[0].set_title("Fe XII 195.119  intensity")

im = axes[1].imshow(v_col, origin="lower", aspect="auto", extent=ext,
                    cmap="RdBu_r", vmin=-15, vmax=15)
axes[1].set_title("Doppler velocity [km/s]\n(blue = toward us)")
fig.colorbar(im, ax=axes[1], label="km/s")
for ax in axes:
    ax.set_xlabel("x [pix]")
axes[0].set_ylabel("y [pix]")
fig.tight_layout()
plt.show()

# %% [markdown]
# 輻射強度のマップとドップラー速度のマップは**似ていません**。明るい＝速い、ではない。
# 明るさとドップラー速度の関係を数字で見てみます。

# %%
q20, q80 = np.nanpercentile(inten, [20, 80])
print("輻射強度で 3 分割したときのドップラー速度 [km/s]（負 = こちらに向かう = 上昇流）")
for name, m in [("暗い 20%", inten < q20),
                ("中間", (inten >= q20) & (inten < q80)),
                ("明るい 20%", inten >= q80)]:
    print(f"  {name:10s} 中央値 {np.nanmedian(v_col[m]):+5.2f}   "
          f"平均 {np.nanmean(v_col[m]):+5.2f}")

# %% [markdown]
# **明るいところはわずかに赤方偏移（下降流）、暗いところはわずかに青方偏移。**
# ただし中央値の差は **1–2 km/s** しかありません。
#
# 一方、画素ごとの振れ幅は ±15 km/s 程度あります。
# つまり**平均的な傾向は小さく、場所ごとのばらつきの方が大きい**。
#
# ここで効いてくるのが 3-2 の話です。**装置由来のずれは 53 km/s ぶん**あり、
# 測りたい信号より大きい。補正とゼロ点の扱いを間違えれば、
# **簡単に嘘の流れが見えます。**

# %% [markdown]
# ## 3-4. 演習
#
# 1. **別の輝線でドップラー速度のマップを作る。** 温度が違えば速度も違うはず。
#    - `Fe XIII 202.044`（1.8 MK、`fe_13_202_044.1c.template.h5`）
#    - `Si VII 275.368`（0.6 MK、`si_07_275_368.1c.template.h5`、moss が見える）
#    ヒント: `fit_region(wvl=..., tmplt_name=..., y0=Y0, y1=Y1)`
#
# 2. **ゼロ点の決め方を変える**（`corr_method="image"` と `"column"`）。
#    マップのどこが変わるか。x 方向の縞が出たり消えたりするはず。
#
# 3. ドップラー速度の**誤差**を見る。`fit.fit["err_params"][..., 1]` が中心波長の誤差。
#    暗いところで速度がどれくらい信用できないか確かめる。

# %%
# 演習 1: ____ を埋めて実行してください（Fe XIII 202.044 のドップラー速度マップ）
#
# fit2, cube2 = fit_region(wvl=____, tmplt_name="____", y0=Y0, y1=Y1)
# cen2 = fit2.fit["params"][..., ____]
# v2 = calc_velocity(cen2, ____, corr_method="column")
# plt.imshow(v2, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-15, vmax=15)
# plt.colorbar(label="km/s"); plt.show()


# %% [markdown]
# **答えは別のノートにあります** →
# [演習の答え](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop_answers.ipynb)
#
# 実行結果も入れてあるので、開くだけで確認できます（走らせる必要はありません）。

# %% [markdown]
# ## まとめ
#
# - ドップラー速度は**中心波長のずれ**。波長画素 1 つが 34 km/s なので、その 1/10 以下を測る
# - **絶対的な波長基準は無い。ゼロ点は自分で決める**（＝仮定を置く）
# - 軌道変動とスリット傾きで **53 km/s ぶん**動く。eispac は補正済み
# - **Solar-C でも同じ。** 「ゼロ点を何と決めたか」を言えることが解析の一部
