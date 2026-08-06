# %% [markdown]
# # 第 2 章: フィットして強度を出す
#
# **35 分**
#
# 山にガウシアンを当てて、**輝線強度という数字**を取り出します。

# %% [markdown]
# ## 2-1. 何を測るのか
#
# $$ I(\lambda) = A\exp\left[-\frac{(\lambda-\lambda_0)^2}{2\sigma^2}\right] + b $$
#
# フィットで 3 つの量が出ます。
#
# | パラメータ | 物理量 | 使う章 |
# |---|---|---|
# | 面積 $A\sigma\sqrt{2\pi}$ | **輝線強度** [erg cm⁻² s⁻¹ sr⁻¹] | 第 2・5 章 |
# | 中心 $\lambda_0$ | **ドップラー速度** | 第 3 章 |
# | 幅 $\sigma$ | **熱運動 + 非熱的速度** | 第 4 章 |
#
# **強度は測定値ではなく、フィットの産物**です。背景をどこに引くか、
# 隣の線をどう扱うか、ガウシアンを何本当てるか —— すべてモデルの仮定です。

# %% [markdown]
# ## 2-2. eispac のテンプレート
#
# 輝線ごとのフィット設定（**テンプレート**）が同梱されています。

# %%
import os
import numpy as np
import matplotlib.pyplot as plt
import eispac

from workshop import EIS_FILE, ensure_eis

ensure_eis()
path = eispac.data.get_fit_template_filepath("fe_12_195_119.2c.template.h5")
tmplt = eispac.read_template(path)

print("ファイル :", os.path.basename(path))
print("line_ids :", tmplt.template["line_ids"])
print("ガウシアンの本数 :", tmplt.template["n_gauss"])
print()
print(f"{'#':>2} {'初期値':>12} {'下限〜上限':>24} {'tied（他に縛る）':>18}")
for i, p in enumerate(tmplt.parinfo):
    lim = f"{p['limits'][0]:.3f} 〜 {p['limits'][1]:.3f}" if p["limited"].any() else "—"
    print(f"{i:2d} {p['value']:12.4f} {lim:>24} {str(p['tied']):>18}")

# %% [markdown]
# パラメータは **[振幅, 中心, 幅] × ガウシアンの本数 + 背景**。
# `2c` は 2 成分なので 3×2 + 1 = 7 個です。
#
# Fe XII 195.119 に 2 成分あるのは、195.179 Å に弱い Fe XII があるためです
# （密度が高いと無視できない）。第 2 成分は `tied` で
# **位置と幅を第 1 成分に縛って**います（原子データで決め打ち）。
#
# ### ★ `line_ids` を必ず確認してから使う
#
# **目的の線が第 0 成分とは限りません。**

# %%
for name in ["fe_12_195_119.2c", "fe_13_203_826.2c", "fe_14_270_519.2c",
             "ar_14_194_396.2c", "ca_14_193_874.2c"]:
    t = eispac.read_template(eispac.data.get_fit_template_filepath(name + ".template.h5"))
    print(f"{name:<20} {[str(s) for s in t.template['line_ids']]}")

# %% [markdown]
# `fe_13_203_826.2c` の**第 0 成分は Fe XII 203.720** です。
# `component=0` と書くと、Fe XIII のつもりで**別のイオンの強度**を取ってしまいます。
# エラーは出ません。**波長で照合して成分番号を決める**のが安全です。

# %%
def pick_component(template, target_wvl):
    """line_ids を見て、目的波長に一番近い成分の番号を返す。"""
    ids = [str(s) for s in template.template["line_ids"]]
    best, bestd = 0, 1e9
    for i, s in enumerate(ids):
        try:
            w = float(s.split()[-1])
        except ValueError:
            continue
        if abs(w - target_wvl) < bestd:
            best, bestd = i, abs(w - target_wvl)
    return best, ids


t = eispac.read_template(eispac.data.get_fit_template_filepath(
    "fe_13_203_826.2c.template.h5"))
print("Fe XIII 203.826 は第", pick_component(t, 203.826)[0], "成分")

# %% [markdown]
# ## 2-3. 領域を決めて、平均してからフィットする
#
# 弱い線は 1 画素では埋もれているので、**まず空間平均して S/N を上げ**、
# それからフィットします（速いので Colab でも快適）。
#
# ここでは活動領域コアの中の **inter-moss 領域**（ループの足元ではなく
# 上部を見ている場所）を使います。

# %%
from workshop import BOX
from fit_box_spectra import average_spectrum      # 欠損値は落としてある（付録 A）

print("箱:", BOX)
wave, inten, sig, npix = average_spectrum(EIS_FILE, 195.119, **BOX)
print(f"平均に使った画素数: {npix}")

fit1 = eispac.fit_spectra(inten, tmplt, wave=wave, errs=sig, ncpu=1,
                          ignore_warnings=True)
wfit, pfit = fit1.get_fit_profile()

plt.figure(figsize=(7, 4.5))
plt.errorbar(wave, inten, yerr=sig, fmt="o", ms=4, label="observed (box average)")
plt.plot(np.ravel(wfit), np.ravel(pfit), "-", lw=2, label="fit")
plt.axvline(195.119, color="r", ls="--", lw=1)
plt.axvline(195.179, color="g", ls="--", lw=1)
plt.xlabel("wavelength [Å]")
plt.ylabel("intensity")
plt.title("Fe XII 195.119")
plt.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 2-4. 強度はガウシアンの面積

# %%
p = np.atleast_1d(fit1.fit["params"]).ravel()
A, lam0, sigma = p[0], p[1], p[2]
I_fit = float(np.atleast_1d(fit1.fit["int"][..., 0]).ravel()[0])

print(f"振幅 A     = {A:10.1f}")
print(f"中心 λ0    = {lam0:10.4f} Å")
print(f"幅   σ     = {sigma:10.4f} Å")
print(f"A σ √(2π)  = {A*sigma*np.sqrt(2*np.pi):10.2f}")
print(f"eispac の int = {I_fit:8.2f}   ← 一致する")
print()
print(f"論文 Table 2 の Fe XII 195.119 = 1147.35")
print(f"比 = {I_fit/1147.35:.2f}   （論文の誤差は ±22%）")

# %% [markdown]
# **論文の値と 1 割の一致。** 独立に処理した結果が合うのは気持ちがよいところです。
#
# 差の主な原因は**測った場所の違い**です。論文は箱の座標を書いていないので、
# 図から読み取るしかありません。
#
# ## 2-5. 演習
#
# 1. **別の輝線でやってみる。** テンプレート名は `scripts/lines_warren2012.py` の
#    一覧にあります。論文 Table 2 の値と比べてみましょう。
#    - `Fe XIII 202.044` → 論文 1076.80
#    - `Fe XV 284.160` → 論文 5931.55
#    - `Ca XV 200.972` → 論文 127.92
# 2. **箱を動かす。** `BOX` の y や x をずらすと強度はどれくらい変わるか。
# 3. `Ca XVII 192.858` をやると論文の 5 倍になります。なぜか考えてみてください
#    （ヒント: 第 1 章の Ca XVII のマップは Fe XII に似ていた。答えは付録 G）

# %%
# 演習 1（TODO を埋める）
#
# from lines_warren2012 import LINES
# for ion, wvl, tname, i_paper, sig_paper in LINES:
#     print(f"{ion:8s} {wvl:8.3f}  {tname}")
#
# tmplt2 = eispac.read_template(eispac.data.get_fit_template_filepath("____"))
# wave2, inten2, sig2, _ = average_spectrum(EIS_FILE, ____, **BOX)
# fit2 = eispac.fit_spectra(inten2, tmplt2, wave=wave2, errs=sig2, ncpu=1,
#                           ignore_warnings=True)
# comp2, ids2 = pick_component(tmplt2, ____)
# print(ids2, "→ 第", comp2, "成分")
# print(float(np.atleast_1d(fit2.fit["int"][..., comp2]).ravel()[0]))

# %% [markdown]
# ## まとめ
#
# - 輝線強度は**ガウシアンの面積** $A\sigma\sqrt{2\pi}$。**フィットの産物**
# - **`line_ids` を確認してから成分番号を決める**
# - 平均してからフィットすると速く、S/N も上がる
# - 論文 Table 2 と 1 割で一致した
