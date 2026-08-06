# %% [markdown]
# # モジュール 2: スペクトル線をフィットして「強度」を取り出す
#
# **所要時間 50 分**
#
# **このノートで身につくこと**
#
# 1. 輝線強度とは**ガウシアンの面積**であり、測定値ではなく**フィットの産物**だと分かる
# 2. eispac のテンプレート機構を使いこなす（`parinfo`, `line_ids`, `tied`）
# 3. **成分の順番の罠**を自分で踏んで確認する（受講者が最も高確率で間違えるところ）
# 4. 論文の手順 = **「箱の中で平均してからフィット」**の理由を、速度と S/N の両面で理解する
# 5. 誤差の入れ方を知る（統計誤差だけでは χ² が発散する）
#
# 前提: モジュール 1。

# %%
!pip install -q eispac

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

# %%
import os
import time
import urllib.request

import numpy as np
import matplotlib.pyplot as plt
import eispac

EIS_FILE = "data/eis/eis_20110702_030712.data.h5"


def ensure(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print("downloading", url)
    urllib.request.urlretrieve(url, path)
    return path


base = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
for ext in ("data", "head"):
    ensure(f"{base}/eis_20110702_030712.{ext}.h5",
           f"data/eis/eis_20110702_030712.{ext}.h5")
print("ok")

# %% [markdown]
# ## 2-1. 何を測るのか
#
# EIS が返すのは各波長画素の強度。そこから輝線強度を出すには
# **ガウシアン + 背景**を当てて、ガウシアンの面積を取る:
#
# $$ I(\lambda) = A \exp\left[-\frac{(\lambda-\lambda_0)^2}{2\sigma^2}\right] + b $$
#
# $$ I_{\rm line} = \int A\,e^{-(\lambda-\lambda_0)^2/2\sigma^2} d\lambda
#    = A\,\sigma\sqrt{2\pi} $$
#
# 得られる 3 つの量:
#
# | パラメータ | 物理量 |
# |---|---|
# | 面積 $A\sigma\sqrt{2\pi}$ | **輝線強度** [erg cm⁻² s⁻¹ sr⁻¹] → DEM に使う |
# | 中心 $\lambda_0$ | **ドップラー速度** |
# | 幅 $\sigma$ | 熱運動 + 非熱的速度（ただし装置幅が支配的） |
#
# **★ 「強度」は測定値ではない。** 背景をどこに引くか、隣の線をどう扱うか、
# ガウシアン何本を当てるか——**全部モデルの仮定**である。
# 論文の値と自分の値が違うとき、まずここを疑う。

# %% [markdown]
# ## 2-2. eispac のテンプレート
#
# eispac には輝線ごとのフィット設定（**テンプレート**）が同梱されている。

# %%
path = eispac.data.get_fit_template_filepath("fe_12_195_119.2c.template.h5")
tmplt = eispac.read_template(path)
print("ファイル:", os.path.basename(path))
print("line_ids:", tmplt.template["line_ids"])
print("n_gauss :", tmplt.template["n_gauss"], "  n_poly:", tmplt.template["n_poly"])
print()
print(f"{'#':>2} {'value':>12} {'fixed':>6} {'limited':>10} {'limits':>22} {'tied':>8}")
for i, p in enumerate(tmplt.parinfo):
    print(f"{i:2d} {p['value']:12.4f} {p['fixed']:6d} {str(p['limited']):>10} "
          f"{str(np.round(p['limits'], 3)):>22} {str(p['tied']):>8}")

# %% [markdown]
# **パラメータの並びは `[振幅, 中心, 幅] × ガウシアンの本数 + 背景の係数`**。
# `2c` テンプレートは 2 成分なので 3×2 + 1 = **7 パラメータ**。
#
# - `.1c` = 1 成分、`.2c` = 2 成分、`.3c` = 3 成分
# - `limited` / `limits` で範囲を縛る（振幅は正、中心は ±0.1 Å など）
# - **`tied`** は「他のパラメータに縛る」しくみ。上の出力では第 2 成分の
#   中心が `p[1]+0.06`（第 1 成分から 0.06 Å 離れた位置）に、
#   幅が `p[2]`（第 1 成分と同じ）に縛られている。
#   つまり **2 本目は原子データで位置を決め打ちしている**。
#   この機構がブレンドを解くときの主役になる（モジュール 8）。
#
# Fe XII 195.119 に**なぜ 2 成分**必要か: 195.179 Å に弱い Fe XII の線があり、
# 密度が高いと無視できない。1 成分で当てると強度が数 % 過大になる。

# %% [markdown]
# ## 2-3. ★ 罠: 成分の順番は波長順とは限らない
#
# 多成分テンプレートで**自分がほしい線が第 0 成分とは限らない**。
# 論文の 22 輝線のうち **4 本**がこれに該当する。必ず `line_ids` で確認すること。

# %%
for name in ["fe_12_195_119.2c", "fe_13_203_826.2c", "fe_14_270_519.2c",
             "ar_14_194_396.2c", "ca_14_193_874.2c"]:
    t = eispac.read_template(eispac.data.get_fit_template_filepath(name + ".template.h5"))
    ids = [str(s) for s in t.template["line_ids"]]
    print(f"{name:<20} {ids}")

# %% [markdown]
# `fe_13_203_826.2c` の**第 0 成分は Fe XII 203.720** であって、
# ほしい Fe XIII 203.826 は第 1 成分。`component=0` と書いたら
# **別のイオンの強度**を DEM に入れてしまう。しかもエラーは出ない。
#
# 手で書くと必ず間違えるので、**波長で自動照合する**関数を使う。

# %%
def pick_component(template, target_wvl):
    """line_ids を見て、目的波長に一番近いガウシアン成分の番号を返す。"""
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
print(pick_component(t, 203.826))      # → (1, [...]) になるはず

# %% [markdown]
# ## 2-4. まず素直に「全ラスターをフィット」してみる
#
# eispac の標準的な使い方。ただし **1 輝線 30720 スペクトルで 2-3 分**かかる。
# 22 輝線なら 1 時間。ここでは**活動領域の部分だけ**（y = 200–320）を試す。

# %%
tmplt = eispac.read_template(eispac.data.get_fit_template_filepath(
    "fe_12_195_119.2c.template.h5"))
cube = eispac.read_cube(EIS_FILE, tmplt.central_wave)

t0 = time.time()
fit = eispac.fit_spectra(cube[200:320, :, :], tmplt, ncpu=2, ignore_warnings=True)
dt = time.time() - t0
n = 120 * 60
print(f"\n{n} スペクトルを {dt:.0f} 秒でフィット "
      f"({1000*dt/n:.1f} ms/スペクトル)")
print(f"→ 全ラスター (512x60) なら約 {dt*512/120/60:.1f} 分、22 輝線なら約 "
      f"{dt*512/120*22/60:.0f} 分")

# %%
m_int = fit.get_map(component=0, measurement="intensity")
m_vel = fit.get_map(component=0, measurement="velocity")

fig, axes = plt.subplots(1, 2, figsize=(9, 7))
d = np.sqrt(np.clip(m_int.data, 0, None))
axes[0].imshow(d, origin="lower", aspect="auto", cmap="inferno",
               vmin=0, vmax=np.nanpercentile(d, 99.5))
axes[0].set_title("Fe XII 195.119  intensity")
axes[1].imshow(m_vel.data, origin="lower", aspect="auto", cmap="RdBu_r",
               vmin=-20, vmax=20)
axes[1].set_title("Doppler velocity [km/s]")
for ax in axes:
    ax.set_xlabel("x [pix]")
axes[0].set_ylabel("y [pix]  (200-320 の範囲)")
fig.tight_layout()
plt.show()

# %% [markdown]
# 速度マップに構造が見える（青 = 上昇流）。これはこれで面白いが、
# **論文がやりたいのは 22 輝線の強度**であって、この速度ではない。
# 22 輝線を全ラスターでフィットするのは時間の無駄になる。

# %% [markdown]
# ## 2-5. 論文の手順: **箱の中で平均してからフィットする**
#
# 論文 §3 はこう書いている:
#
# > we extract the EIS data from each spectral window in the selected field of view
# > and **average them together** (missing data are not included in the average)
# > to form high signal-to-noise line profiles, which are then fit with
# > single Gaussians
#
# **順番が「平均 → フィット」である**ことが重要。理由は 2 つ:
#
# 1. **速い**。240 画素を平均すれば、フィットするのは 1 本のプロファイルだけ。
#    22 輝線でも数秒で終わる（Colab で決定的）。
# 2. **S/N が上がる**。弱い線（Ca XVI は最も明るい線の 1/200）は
#    1 画素では埋もれている。平均して初めてフィットできる。
#
# 逆順（各画素をフィットしてから平均）とは**同じにならない**。
# フィットは非線形なので、平均と可換ではない。
# 論文がどちらをやったかを読み取ることが、再現の前提になる。
#
# 使う箱は **モジュール 4 で選ぶ**。ここでは結果を先取りして使う。

# %%
BOX = dict(y0=244, y1=274, x0=32, x1=40)     # モジュール 4 で決めた inter-moss 箱


def average_spectrum(datafile, wvl, y0, y1, x0, x1):
    """箱の中でスペクトルを平均する。欠損サンプルは平均に入れない。

    ★ 欠損は NaN ではなく**大きな負のフラグ値**（モジュール 1 参照）。
      eispac が立てる cube.mask で必ず落とす。
    """
    c = eispac.read_cube(datafile, wvl)
    data = c.data[y0:y1, x0:x1, :]
    errs = c.uncertainty.array[y0:y1, x0:x1, :]
    wave = c.wavelength[y0:y1, x0:x1, :]
    bad = np.asarray(c.mask[y0:y1, x0:x1, :], dtype=bool)

    good = np.isfinite(data) & ~bad
    n = good.sum(axis=(0, 1))                            # 波長ごとの有効画素数
    inten = np.nansum(np.where(good, data, 0), axis=(0, 1)) / np.maximum(n, 1)
    # 平均の誤差 = sqrt(Σσ²)/N
    sig = np.sqrt(np.nansum(np.where(good, errs**2, 0), axis=(0, 1))) / np.maximum(n, 1)
    inten[n == 0], sig[n == 0] = np.nan, np.nan
    return np.nanmean(wave, axis=(0, 1)), inten, sig, int(np.median(n))


wave, inten, sig, npix = average_spectrum(EIS_FILE, 195.119, **BOX)
print(f"箱 y=[{BOX['y0']}:{BOX['y1']}] x=[{BOX['x0']}:{BOX['x1']}] "
      f"= {(BOX['y1']-BOX['y0'])*(BOX['x1']-BOX['x0'])} 画素")
print(f"平均に使えた画素数（中央値）: {npix}")
print(f"ピーク強度 {np.nanmax(inten):.0f} ± {sig[np.nanargmax(inten)]:.1f}  "
      f"→ 相対誤差 {100*sig[np.nanargmax(inten)]/np.nanmax(inten):.2f}%")

# %%
t0 = time.time()
fit1 = eispac.fit_spectra(inten, tmplt, wave=wave, errs=sig, ncpu=1,
                          ignore_warnings=True)
print(f"1 本のプロファイルのフィット: {1000*(time.time()-t0):.0f} ms")

wfit, pfit = fit1.get_fit_profile()          # 細かい波長グリッドでのモデル曲線
plt.figure(figsize=(7, 4.5))
plt.errorbar(wave, inten, yerr=sig, fmt="o", ms=4, label="box-averaged data")
plt.plot(np.ravel(wfit), np.ravel(pfit), "-", lw=2, label="fit (2 Gaussians + bg)")
plt.axvline(195.119, color="r", ls="--", lw=1)
plt.axvline(195.179, color="g", ls="--", lw=1)
plt.xlabel("wavelength [Å]")
plt.ylabel("intensity")
plt.title("Fe XII 195.119, averaged over the box")
plt.legend()
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 2-6. 強度がガウシアンの面積であることを数値で確かめる

# %%
p = np.atleast_1d(fit1.fit["params"]).ravel()
A, lam0, sigma = p[0], p[1], p[2]
I_fit = float(np.atleast_1d(fit1.fit["int"][..., 0]).ravel()[0])
print(f"A     = {A:12.2f}")
print(f"λ0    = {lam0:12.4f} Å")
print(f"σ     = {sigma:12.4f} Å")
print(f"A σ √(2π) = {A*sigma*np.sqrt(2*np.pi):10.2f}")
print(f"eispac の int = {I_fit:10.2f}   ← 一致する")
print()
print(f"背景 b = {p[-1]:.1f}  （ピークの {100*p[-1]/A:.0f}%）")
print(f"論文 Table 2 の Fe XII 195.119 = 1147.35 → 比 {I_fit/1147.35:.3f}")

# %% [markdown]
# **σ = 0.028 Å の意味**: EIS の装置幅は σ ≈ 0.027 Å（1″ スリット）。
# つまり測っている幅の**ほとんどが装置由来**で、
# プラズマの熱運動 + 非熱的速度はその上に乗るわずかな超過分。
# 線幅から速度を出すには装置幅を正確に知る必要がある（今回は使わない）。

# %% [markdown]
# ## 2-7. ★ 誤差の入れ方（ここを間違えると DEM で破綻する）
#
# 上で見たとおり、240 画素を平均した後の**統計誤差は 0.2%** しかない。
# 一方、論文が使っている誤差は **22%**。

# %%
chi2 = float(np.atleast_1d(fit1.fit["chi2"]).ravel()[0])
print(f"このフィットの χ² = {chi2:.0f}   (データ点 {len(wave)}, パラメータ 7)")
print(f"χ²_red = {chi2/(len(wave)-7):.0f}   ← 1 のはずが桁違い")

# %% [markdown]
# **なぜ χ² が桁違いに大きいのか**
#
# 統計誤差だけを使うと、誤差が小さすぎて「ガウシアン + 直線背景」という
# **モデルのわずかな不完全さ**が全部 χ² に化ける。
# 実際のスペクトルには弱いブレンドや非対称性があり、
# 0.2% の精度でガウシアンに一致することはない。
#
# **論文の 22% はどこから来るか**
#
# - EIS の**絶対較正**の不確かさ（打ち上げ前較正で ~20%、劣化補正でさらに増える）
# - これは**系統誤差**なので、画素を平均しても減らない
#
# → DEM を解くときは **σ_I = max(統計誤差, 0.22 × I)** のように
#   **較正誤差を床として入れる**。これをやらないと χ² が発散して
#   DEM インバージョンが収束しない（モジュール 7 で実際に見る）。
#
# **★ 教訓**: 「誤差が小さい」ことは良いことではない。
# **何の誤差を見積もっているか**を意識する。

# %% [markdown]
# ## 2-8. 22 輝線をまとめてフィットする
#
# 論文 Table 2 の 22 輝線と、対応する eispac テンプレートの表は
# リポジトリの `scripts/lines_warren2012.py` にある。

# %%
import sys
sys.path.insert(0, "scripts")
from lines_warren2012 import LINES        # (ion, 波長, テンプレート, 論文値, 論文σ)

print(f"{'line':<16}{'template':<26}{'I_paper':>9}")
for ion, wvl, tname, ip, sp in LINES[:5]:
    print(f"{ion+' '+f'{wvl:.3f}':<16}{tname.replace('.template.h5',''):<26}{ip:9.2f}")
print(f"... 全 {len(LINES)} 本")

# %%
t0 = time.time()
rows = []
for ion, wvl, tname, i_paper, sig_paper in LINES:
    w, I, s, npix = average_spectrum(EIS_FILE, wvl, **BOX)
    t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
    comp, ids = pick_component(t, wvl)
    f = eispac.fit_spectra(I, t, wave=w, errs=s, ncpu=1, ignore_warnings=True)
    i_fit = float(np.atleast_1d(f.fit["int"][..., comp]).ravel()[0])
    rows.append((ion, wvl, i_fit, i_paper, i_fit / i_paper, comp, ids[comp]))
print(f"\n22 輝線を {time.time()-t0:.0f} 秒でフィットした\n")

print(f"{'line':<16}{'I_fit':>10}{'I_paper':>10}{'ratio':>7}  component")
for ion, wvl, i_fit, i_paper, r, comp, cid in rows:
    print(f"{ion+' '+f'{wvl:.3f}':<16}{i_fit:10.2f}{i_paper:10.2f}{r:7.2f}  [{comp}] {cid}")

# %%
import csv
os.makedirs("work", exist_ok=True)
with open("work/box_intensities.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["ion", "wvl", "I_fit", "I_paper", "ratio"])
    w.writerows([(r[0], r[1], r[2], r[3], r[4]) for r in rows])
print("wrote work/box_intensities.csv  （モジュール 5, 6, 7 で使う）")

# %% [markdown]
# ここまでで**論文 Table 2 と比べられる数字**が出た。
# 中身の議論（どれが合っていて、どれが合わないか）は**モジュール 5** で行う。
#
# ここで先に 1 つだけ言っておくと、**Ca XVII 192.858 の ratio が 5 前後**に
# なっているはずである。これはフィットの失敗ではなく
# **ブレンド**（Fe XI と O V が混ざっている）で、モジュール 8 で解く。
#
# 同じ処理は `scripts/fit_box_spectra.py` にまとめてある:
#
# ```bash
# python scripts/fit_box_spectra.py data/eis/eis_20110702_030712.data.h5 244 274 32 40
# ```

# %% [markdown]
# ## 2-9. 演習
#
# 1. **成分の順番を間違えるとどうなるか**を体験する。
#    上のループで `comp` を `0` に固定して Fe XIII 203.826 の値を見る。
#    論文値との比がどう変わるか。
# 2. **マスクを外すとどうなるか**。`average_spectrum` の `& ~bad` を消して
#    22 輝線を再実行し、どの線が何 % 変わるか調べる
#    （ヒント: 弱い線ほど効く。Ca XVII、Ca XVI、Ar XIV）。
# 3. **箱の大きさを変える**。`BOX` の y 範囲を 2 倍にして、
#    統計誤差と ratio がどう変わるか。誤差は減るが ratio は良くなるか？
# 4. 2-5 の逆順（各画素をフィットしてから平均）を Fe XII で試して、
#    「平均 → フィット」との差を測る（`scripts/fit_perpixel_box.py` に実装がある）。
#
# ## まとめ
#
# - 輝線強度は**ガウシアンの面積** $A\sigma\sqrt{2\pi}$。**フィットの産物**である
# - **成分の順番は波長順ではない**。必ず `line_ids` で照合する
# - 論文の手順は **「箱で平均 → フィット」**。速さと S/N の両方の理由がある
# - **統計誤差だけでは足りない**。較正の系統誤差 22% を床として入れる
#
# 次（モジュール 3）では、EIS では測れない 7 MK のプラズマを
# **AIA 94 Å から取り出す**。
