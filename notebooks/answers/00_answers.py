# %% [markdown]
# # 演習の答え
#
# [講習会ノートブック](https://colab.research.google.com/github/hottahd/EIS_practice/blob/main/notebooks/EIS_workshop.ipynb)
# の演習の答えです。
#
# **実行結果を入れてあるので、開くだけで確認できます。**
# 自分で動かしたい場合は、下の準備セルから順に実行してください
# （本編と同じデータを取り直すので数分かかります）。

# %%
!pip install -q eispac fiasco demregpy

# %%
import os
import subprocess
import sys
import urllib.request

REPO = "https://github.com/hottahd/EIS_practice.git"
if not os.path.exists("scripts/lines_warren2012.py"):
    if not os.path.exists("EIS_practice"):
        print("教材リポジトリを取得中 ...")
        subprocess.run(["git", "clone", "-q", REPO], check=True)
    os.chdir("EIS_practice")
sys.path.insert(0, "scripts")

import numpy as np
import matplotlib.pyplot as plt
# データを取る（本編と同じ。既にあれば何もしない）
BASE = "https://eis.nrl.navy.mil/level1/hdf5/2011/07/02"
EIS_FILE = "data/eis/eis_20110702_030712.data.h5"
REGION = dict(y0=244, y1=274, x0=32, x1=40)   # 本編と同じ領域

os.makedirs("data/eis", exist_ok=True)
for ext in ("data", "head"):
    name = f"eis_20110702_030712.{ext}.h5"
    if not os.path.exists(f"data/eis/{name}"):
        urllib.request.urlretrieve(f"{BASE}/{name}", f"data/eis/{name}")

print("作業ディレクトリ:", os.getcwd())

# %% [markdown]
# ---
# ## 第 3 章 演習 1: 別の輝線で速度マップを作る
#
# Fe XIII 202.044（1.8 MK）で同じことをして、Fe XII 195.119（1.6 MK）と比べます。

# %%
from eispac.instr import calc_velocity

Y0, Y1 = 180, 340
ext = [0, 60, Y0, Y1]

import eispac

tmplt = eispac.read_template(
    eispac.data.get_fit_template_filepath("fe_12_195_119.2c.template.h5"))
cube = eispac.read_cube(EIS_FILE, tmplt.central_wave)
fit = eispac.fit_spectra(cube[Y0:Y1, :, :], tmplt, ncpu=2, ignore_warnings=True)
v_col = calc_velocity(fit.fit["params"][..., 1], 195.119, corr_method="column")

# %%
tmplt2 = eispac.read_template(
    eispac.data.get_fit_template_filepath("fe_13_202_044.1c.template.h5"))
cube2 = eispac.read_cube(EIS_FILE, tmplt2.central_wave)
fit2 = eispac.fit_spectra(cube2[Y0:Y1, :, :], tmplt2, ncpu=2, ignore_warnings=True)
cen2 = fit2.fit["params"][..., 1]
v2 = calc_velocity(cen2, 202.044, corr_method="column")

fig, axes = plt.subplots(1, 2, figsize=(9, 7))
for ax, (d, t) in zip(axes, [(v_col, "Fe XII 195.119 (1.6 MK)"),
                             (v2, "Fe XIII 202.044 (1.8 MK)")]):
    im = ax.imshow(d, origin="lower", aspect="auto", extent=ext,
                   cmap="RdBu_r", vmin=-15, vmax=15)
    ax.set_title(t + "  velocity [km/s]", fontsize=10)
    ax.set_xlabel("x [pix]")
axes[0].set_ylabel("y [pix]")
fig.colorbar(im, ax=axes, label="km/s")
plt.show()

print(f"Fe XII  : 中央値 {np.nanmedian(v_col):+.2f}  "
      f"5–95% [{np.nanpercentile(v_col,5):+.1f}, {np.nanpercentile(v_col,95):+.1f}] km/s")
print(f"Fe XIII : 中央値 {np.nanmedian(v2):+.2f}  "
      f"5–95% [{np.nanpercentile(v2,5):+.1f}, {np.nanpercentile(v2,95):+.1f}] km/s")
print("→ 温度が近い 2 本なので、速度の分布もよく似る")

# %% [markdown]
# ---
# ## 第 4 章 演習 1: 装置幅を定数で代用するとどうなるか

# %%
LAM0, C_KMS = 195.119, 2.998e5
K_B, AMU = 1.380649e-16, 1.66054e-24

sig_obs = fit.fit["params"][..., 2]
inten = fit.fit["int"][..., 0]
fwhm_inst = np.asarray(cube.meta["slit_width"], float)[Y0:Y1]
sig_inst = fwhm_inst / (2 * np.sqrt(2 * np.log(2)))
sig_th = LAM0 * np.sqrt(K_B * 10**6.2 / (56.0 * AMU)) / 1e5 / C_KMS


def nonthermal_velocity(sig_obs, sig_inst, sig_th, lam0=LAM0):
    excess = sig_obs**2 - sig_inst**2 - sig_th**2
    excess = np.where(excess > 0, excess, np.nan)
    return np.sqrt(2) * C_KMS / lam0 * np.sqrt(excess)


xi = nonthermal_velocity(sig_obs, sig_inst[:, None], sig_th)
bright = inten > np.nanpercentile(inten, 60)

# %%
xi_const = nonthermal_velocity(sig_obs, sig_inst.mean(), sig_th)
diff = xi - xi_const

plt.figure(figsize=(4.5, 7))
plt.imshow(np.where(bright, diff, np.nan), origin="lower", aspect="auto",
           extent=ext, cmap="coolwarm", vmin=-10, vmax=10)
plt.colorbar(label="km/s")
plt.title("position-dependent minus constant")
plt.xlabel("x [pix]"); plt.ylabel("y [pix]")
plt.tight_layout()
plt.show()

print(f"差の範囲: {np.nanmin(diff[bright]):+.1f} 〜 {np.nanmax(diff[bright]):+.1f} km/s")
print("→ y 方向に系統的なパターンが出る。装置幅の y 依存がそのまま乗っている")

# %% [markdown]
# ---
# ## 第 5 章 演習 1: 高温側の輝線を落とすとどうなるか

# %%
from demregpy import dn2dem
import sys
sys.path.insert(0, "notebooks/src")


def read_gofnt(path):
    lines = open(path).readlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# nT nline"))
    nT, nline = (int(v) for v in lines[i + 1].split())
    k = i + 2

    def skip(tag):
        nonlocal k
        while not lines[k].startswith(tag):
            k += 1
        k += 1

    def take(n):
        nonlocal k
        out = []
        while len(out) < n:
            out += [float(x) for x in lines[k].split()]
            k += 1
        return np.array(out[:n])

    skip("# logT")
    logT = take(nT)
    skip("# ion")
    names, wvl = [], []
    for _ in range(nline):
        p = lines[k].split()
        names.append(" ".join(p[:-2]))
        wvl.append(float(p[-2]))
        k += 1
    skip("# G(T)")
    G = np.array([take(nT) for _ in range(nline)])
    return logT, names, np.array(wvl), G


logT, names, wvl, G = read_gofnt("work/gofnt_chianti901_005.txt")
iobs = np.zeros(len(names))
from lines_warren2012 import LINES, pick_component
from fit_box_spectra import average_spectrum

for ion, w, tname, i_paper, sig_paper in LINES:
    wave, inten, sig, _ = average_spectrum(EIS_FILE, w, **REGION)
    t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
    comp, _ = pick_component(t, w)
    f = eispac.fit_spectra(inten, t, wave=wave, errs=sig, ncpu=1, ignore_warnings=True)
    k = int(np.argmin(np.abs(wvl - w)))
    if abs(wvl[k] - w) < 0.01:
        iobs[k] = float(np.atleast_1d(f.fit["int"][..., comp]).ravel()[0])
iobs[int(np.argmin(np.abs(wvl - 192.858)))] = 0.0

keep = (logT >= 5.5) & (logT <= 7.1)
lt, Gk = logT[keep], G[:, keep]
dlt = lt[1] - lt[0]
T = 10**lt
tedges = 10 ** np.append(lt - dlt / 2, lt[-1] + dlt / 2)
sel = iobs > 0
dem, _, _, _, _ = dn2dem(iobs[sel], 0.22 * iobs[sel], (Gk[sel] / (4 * np.pi)).T,
                         lt, tedges, max_iter=30, warn=False)
em = np.atleast_1d(np.squeeze(dem)) * T * np.log(10) * dlt

# %%
sel2 = sel.copy()
for bad in ["Ca XIV", "Ca XV", "Ca XVI"]:
    sel2[[i for i, n in enumerate(names) if n == bad]] = False

dem2, _, _, chi2_2, _ = dn2dem(iobs[sel2], 0.22 * iobs[sel2],
                               (Gk[sel2] / (4 * np.pi)).T, lt, tedges,
                               max_iter=30, warn=False)
em2 = np.atleast_1d(np.squeeze(dem2)) * T * np.log(10) * dlt

plt.figure(figsize=(7.5, 5))
plt.plot(lt, em, "o-", lw=2, ms=4, label=f"all lines ({int(sel.sum())})")
plt.plot(lt, em2, "s--", lw=2, ms=4, label=f"without Ca ({int(sel2.sum())})")
plt.yscale("log"); plt.xlim(5.6, 7.1); plt.ylim(1e24, 1e28)
plt.xlabel("log T [K]"); plt.ylabel(r"EM per bin [cm$^{-5}$]")
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout(); plt.show()

print(f"全輝線     : ピーク logT {lt[int(np.nanargmax(em))]:.2f}"
      f"  ({10**lt[int(np.nanargmax(em))]/1e6:.1f} MK)")
print(f"Ca を除く  : ピーク logT {lt[int(np.nanargmax(em2))]:.2f}"
      f"  ({10**lt[int(np.nanargmax(em2))]/1e6:.1f} MK)")
print("→ 高温側を拘束する輝線を失うと、その温度域は決まらなくなる")
