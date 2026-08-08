# %% [markdown]
# # 第 4 章: 線幅から非熱的速度を出す
#
# **30 分**
#
# フィットから出る 3 つ目の量が**線幅**。ここから
# **非熱的速度**（熱運動では説明できない速度成分）が出る。
#
# 波・乱流・視線上に重なった細かい流れ —— コロナ加熱の議論に直結する量で、
# Solar-C の主戦場のひとつ。

# %% [markdown]
# ## 4-1. 線幅は 3 つの足し算
#
# 観測される線幅は、次の 3 つが**二乗和**で足さったもの:
#
# $$ \sigma_{\rm obs}^2 = \sigma_{\rm inst}^2 + \sigma_{\rm th}^2 + \sigma_{\rm nonth}^2 $$
#
# | | 中身 |
# |---|---|
# | $\sigma_{\rm inst}$ | **装置の幅**。EIS では最大の寄与 |
# | $\sigma_{\rm th}$ | **熱運動**。イオンが重いほど狭い。$\sigma_{\rm th} = \frac{\lambda_0}{c}\sqrt{kT/M}$ |
# | $\sigma_{\rm nonth}$ | **それ以外**。波・乱流・視線上の速度の重なり |
#
# 非熱的速度 $\xi$ は $\sigma_{\rm nonth} = \frac{\lambda_0}{c}\frac{\xi}{\sqrt{2}}$ で定義する
# （$\xi$ は最確速度）。つまり
#
# $$ \xi = \sqrt{2}\,\frac{c}{\lambda_0}
#    \sqrt{\sigma_{\rm obs}^2 - \sigma_{\rm inst}^2 - \sigma_{\rm th}^2} $$
#
# **引き算なので、装置幅を間違えると答えが大きく動く。**

# %%
import numpy as np
import matplotlib.pyplot as plt
import eispac

Y0, Y1 = 180, 340

# 第 3 章と同じフィット結果を使う。
# 上から順に実行していれば `fit` が残っているので、その場合は解き直さない
try:
    fit
except NameError:
    tmplt = eispac.read_template(
        eispac.data.get_fit_template_filepath("fe_12_195_119.2c.template.h5"))
    cube = eispac.read_cube(EIS_FILE, tmplt.central_wave)
    fit = eispac.fit_spectra(cube[Y0:Y1, :, :], tmplt, ncpu=2, ignore_warnings=True)

sig_obs = fit.fit["params"][..., 2]      # フィットで得た σ [Å]
inten = fit.fit["int"][..., 0]
print(f"観測された σ: 中央値 {np.nanmedian(sig_obs):.4f} Å")

# %% [markdown]
# ## 4-2. 装置幅はスリット上の位置で変わる
#
# 第 1 章で見た「山の幅はほぼ装置で決まっている」の中身がこれです。
# EIS の**波長分解能**は `cube.meta['slit_width']` に入っています（**FWHM**、単位 Å）。
# **スリットに沿って一定ではありません。**

# %%
fwhm_inst = np.asarray(cube.meta["slit_width"], float)[Y0:Y1]
sig_inst = fwhm_inst / (2 * np.sqrt(2 * np.log(2)))      # FWHM → σ

print(f"波長分解能 FWHM: {fwhm_inst.min():.4f} – {fwhm_inst.max():.4f} Å "
      f"（この範囲で {100*(fwhm_inst.max()/fwhm_inst.min()-1):.0f}% 変わる）")
print(f"σ に直すと : {sig_inst.min():.4f} – {sig_inst.max():.4f} Å")
print(f"観測の σ   : 中央値 {np.nanmedian(sig_obs):.4f} Å")
print("\n→ 観測した幅のほとんどが装置の幅。差の部分を取り出す作業になる。")

plt.figure(figsize=(4.5, 5))
plt.plot(fwhm_inst, np.arange(Y0, Y1), lw=1.5)
plt.xlabel("instrumental width FWHM [Å]")
plt.ylabel("y [pix] (along the slit)")
plt.title("instrumental width varies along the slit")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 4-3. 熱幅を引く
#
# 熱運動の幅は**イオンの質量と温度**で決まる。
# 温度は測れないので、**その輝線の形成温度を仮定する**（Fe XII なら log T = 6.2）。
#
# これは仮定なので、後で変えてみて影響を見る（演習 2）。

# %%
LAM0, C_KMS = 195.119, 2.998e5
M_AMU, LOGT = 56.0, 6.2                  # Fe XII: 鉄の質量数と形成温度
K_B, AMU = 1.380649e-16, 1.66054e-24     # cgs

v_th = np.sqrt(K_B * 10**LOGT / (M_AMU * AMU)) / 1e5      # km/s
sig_th = LAM0 * v_th / C_KMS

print(f"熱速度   sqrt(kT/M) = {v_th:.1f} km/s   （logT={LOGT}, 鉄）")
print(f"熱幅     σ_th       = {sig_th:.4f} Å")
print(f"装置幅   σ_inst     = {np.median(sig_inst):.4f} Å  ← こちらの方がずっと大きい")

# %% [markdown]
# ## 4-4. 非熱的速度のマップ

# %%
def nonthermal_velocity(sig_obs, sig_inst, sig_th, lam0=LAM0):
    """観測の σ から非熱的速度 ξ [km/s] を出す。引けない画素は NaN。"""
    excess = sig_obs**2 - sig_inst**2 - sig_th**2
    excess = np.where(excess > 0, excess, np.nan)
    return np.sqrt(2) * C_KMS / lam0 * np.sqrt(excess)


xi = nonthermal_velocity(sig_obs, sig_inst[:, None], sig_th)   # ← y ごとの装置幅

bright = inten > np.nanpercentile(inten, 60)      # 暗い場所は幅が信用できない
print(f"非熱的速度（明るい画素）: 中央値 {np.nanmedian(xi[bright]):.1f} km/s   "
      f"5–95% [{np.nanpercentile(xi[bright],5):.1f}, {np.nanpercentile(xi[bright],95):.1f}]")
print(f"引き算が成立しなかった画素: {100*np.isnan(xi[bright]).mean():.0f}%")

# %%
ext = [0, 60, Y0, Y1]
fig, axes = plt.subplots(1, 2, figsize=(9, 7))
d = np.sqrt(np.clip(inten, 0, None))
axes[0].imshow(d, origin="lower", aspect="auto", extent=ext, cmap="inferno",
               vmin=0, vmax=np.nanpercentile(d, 99.5))
axes[0].set_title("Fe XII 195.119  intensity")

im = axes[1].imshow(np.where(bright, xi, np.nan), origin="lower", aspect="auto",
                    extent=ext, cmap="viridis", vmin=5, vmax=35)
axes[1].set_title("non-thermal velocity [km/s]")
fig.colorbar(im, ax=axes[1], label="km/s")
for ax in axes:
    ax.set_xlabel("x [pix]")
axes[0].set_ylabel("y [pix]")
fig.tight_layout()
plt.show()

# %% [markdown]
# 輻射強度との関係を数字で見ます。

# %%
q20, q80 = np.nanpercentile(inten, [20, 80])
print("輻射強度で 3 分割したときの非熱的速度 [km/s]")
for name, m in [("暗い 20%", inten < q20),
                ("中間", (inten >= q20) & (inten < q80)),
                ("明るい 20%", inten >= q80)]:
    print(f"  {name:10s} 中央値 {np.nanmedian(xi[m]):5.1f}")
good = np.isfinite(xi) & (inten > 0)
print(f"\n相関係数 (log I, ξ) = {np.corrcoef(np.log10(inten[good]), xi[good])[0,1]:+.2f}")

# %% [markdown]
# - 全体の中央値は **18 km/s 前後**。活動領域として典型的な値
# - **明るいところほど非熱的速度は大きい**（暗 14 → 明 20 km/s）。
#   ただし相関は +0.1 程度で**弱い**
# - 暗い画素は線幅の測定誤差が大きく、引き算が破綻して NaN になる
#   （上の統計で暗い側の値が低めに出るのは、この効果も混ざっている）
#
# 「非熱的」の中身は 1 つではない。波かもしれないし、視線上に速度の違う
# 構造がいくつも重なっているだけかもしれない。
# **0.4″ の Solar-C で見ると、この一部は「分解できていなかっただけ」に変わるはず。**
# それを確かめるのが EUVST の仕事のひとつ。

# %% [markdown]
# ## 4-5. 演習
#
# 1. **装置幅を定数（平均値）で代用するとどうなるか。**
#    `nonthermal_velocity(sig_obs, sig_inst.mean(), sig_th)` に変えて引き算し、
#    差のマップを描く。**y 方向に系統的なパターン**が出るはず。
# 2. **形成温度の仮定を変える。** `LOGT` を 6.0 や 6.4 にすると ξ はどれだけ動くか。
#    「仮定が結果に効く」量なのかどうかを自分で確かめる。
# 3. 軽い元素の輝線（例: `Si X 258.375`、Si は鉄より軽い）で同じことをする。
#    熱幅の寄与が大きくなるので、仮定の効き方も変わる。

# %%
# 演習 1: ____ を埋めて実行してください（装置幅を定数で代用するとどうなるか）
#
# xi_const = nonthermal_velocity(sig_obs, ____, sig_th)
# diff = xi - xi_const
# plt.imshow(np.where(bright, diff, np.nan), origin="lower", aspect="auto",
#            extent=ext, cmap="coolwarm", vmin=-10, vmax=10)
# plt.colorbar(label="km/s"); plt.show()


# %% [markdown]
# **答えは別のノートにあります** →
# [演習の答え](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop_answers.ipynb)
#
# 実行結果も入れてあるので、開くだけで確認できます（走らせる必要はありません）。

# %% [markdown]
# ## まとめ
#
# - 線幅は **装置 + 熱運動 + それ以外** の二乗和。**引き算で非熱的速度を取り出す**
# - **装置幅はスリット上の位置で変わる**（この観測で 7%）。`slit_width` を使う
# - 熱幅には**温度の仮定**が入る
# - 活動領域の非熱的速度は **20 km/s 前後**
