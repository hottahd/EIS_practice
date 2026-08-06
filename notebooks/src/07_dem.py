# %% [markdown]
# # モジュール 7: DEM インバージョン
#
# **所要時間 60 分**
#
# **このノートで身につくこと**
#
# 1. DEM 逆問題が **ill-posed** であることを、特異値分解で**数値として**見る
# 2. 誤差の床（較正の系統誤差）を入れる理由を理解する
# 3. `demregpy`（正則化）で実際に解き、論文 Table 2 の R 列と比べる
# 4. ★ **手法・設定で答えがどれだけ動くか**を測る（これがこの章の核心）
# 5. 傾き α から加熱の物理を議論する
#
# 前提: モジュール 2（強度）、4（AIA）、6（G(T)）。

# %%
!pip install -q demregpy

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
import csv
import os

import numpy as np
import matplotlib.pyplot as plt

GOFNT_FINE = "work/gofnt_chianti901_005.txt"   # 0.05 dex（demregpy 用）
GOFNT = "work/gofnt_chianti901.txt"            # 0.10 dex（MCMC 用）
AIARESP = "work/aia94_fe18_response.txt"
MCMC = "work/mcmc_dem_result.txt"              # PINTofALE MCMC_DEM の結果（同梱）


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


logT, names, wvl, G = read_gofnt(GOFNT_FINE)
keep = (logT >= 5.5) & (logT <= 7.1)
logT, G = logT[keep], G[:, keep]
dlt = logT[1] - logT[0]
print(f"温度格子: logT {logT[0]:.2f}–{logT[-1]:.2f}, {len(logT)} ビン, {dlt:.2f} dex 刻み")

# %% [markdown]
# ## 7-1. ★ なぜ ill-posed なのか —— 特異値で見る
#
# 解きたいのは
#
# $$ I_\lambda = \frac{1}{4\pi}\int G_\lambda(T)\,\xi(T)\,dT $$
#
# で、$\xi(T)$ が未知（**第一種 Fredholm 積分方程式**）。
# 温度を離散化すると、ただの行列方程式 $\mathbf{I} = \mathsf{A}\,\boldsymbol{\xi}$ になる。
# $\mathsf{A}$ は「輝線 × 温度ビン」の行列で、中身は $G/(4\pi)$。
#
# **この行列がどれくらい"効いて"いるかは特異値を見れば分かる。**

# %%
A = G / (4 * np.pi)                 # (nline, nT)
sv = np.linalg.svd(A, compute_uv=False)
print(f"行列の形     : {A.shape}  (輝線 x 温度ビン)")
print(f"条件数       : {sv[0]/sv[-1]:.2e}")
print(f"特異値（上位）: {np.round(sv[:8]/sv[0], 4)}")
print(f"最大の 1/1000 以上ある特異値の本数 = {int((sv > 1e-3*sv[0]).sum())}")

plt.figure(figsize=(6, 4))
plt.semilogy(sv / sv[0], "o-")
plt.axhline(1e-3, color="r", ls="--", label="1/1000 of the largest")
plt.xlabel("index")
plt.ylabel("singular value / max")
plt.title("singular values of the response matrix")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# **22 本の輝線を測ったのに、実質的な自由度は 12 程度しかない。**
#
# 理由はモジュール 6 で見たとおり、**G(T) が幅広く重なっている**から。
# 隣り合う温度ビンの応答がほとんど同じなので、両者を区別できない。
#
# 数学的には、特異値が急速に 0 に近づくと逆作用素が非有界になる。
# 実際問題としては:
#
# - $\xi(T)$ に細かい構造を入れても、$G$ との積分で均されて $I_\lambda$ に現れない
# - 逆に、$I_\lambda$ をわずかに動かすと $\xi$ の細かい構造は激しく変わる
#   = **観測誤差が解に増幅される**
#
# → だから**何らかの追加の仮定（正則化・事前分布）が必ず要る**。
#   「DEM を解いた」と言うときは、**何を仮定したか**を必ずセットで言う。

# %% [markdown]
# ## 7-2. 観測強度と誤差
#
# モジュール 2 で出した 22 輝線の強度を読む。
#
# - **Ca XVII 192.858 は外す**（ブレンドしたまま。モジュール 8 で分離したら戻す）
# - 誤差は **22%**（論文と同じ）。統計誤差は 0.2% しかないが、
#   それを使うと χ² が発散する（モジュール 2 参照）

# %%
# モジュール 2 の出力を読む。Colab はノート 1 冊ごとに VM が変わるので、
# 無ければその場で作り直す（`scripts/workshop.py` にまとめてある。10 秒ほど）。
from workshop import box_intensities

iobs = np.zeros(len(names))
for ion, w, i_fit, i_paper, ratio in box_intensities():
    k = int(np.argmin(np.abs(wvl - w)))
    if abs(wvl[k] - w) < 0.01:
        iobs[k] = i_fit
iobs[int(np.argmin(np.abs(wvl - 192.858)))] = 0.0        # Ca XVII を外す

ok = iobs > 0
labels = [f"{names[i]} {wvl[i]:.3f}" for i in np.where(ok)[0]]
dn = iobs[ok]
edn = 0.22 * dn                                          # ★ 較正の系統誤差を床に
tresp = (A[ok]).T                                        # (nT, nf)
print(f"EIS から {len(dn)} 本")

# %% [markdown]
# ### AIA Fe XVIII を拘束に加える
#
# モジュール 6 で見たとおり、**log T > 6.9 は EIS だけでは決まらない**。
# 入れないと DEM が高温側に漏れて傾き β が出鱈目になる
# （PINTofALE の文書が **"toothpaste tube effect"** と呼ぶ現象）。
#
# **★ 公式の AIA 94 Å 応答は使えない。** 低温線の寄与が入っているため。
# **Fe XVIII だけの応答関数**を作ってある
# （`scripts/aia94_fe18_response.py`。CHIANTI の Fe XVIII 582 本を
# aiapy の波長応答に畳み込んだもの。実際は 93.932 Å の 1 本が寄与の 100%）。

# %%
d = np.loadtxt(AIARESP)
R = np.interp(logT, d[:, 0], d[:, 1])
print(f"R(T) のピーク = {R.max():.3e} DN cm^5 s^-1 pix^-1 "
      f"at logT {logT[int(np.argmax(R))]:.2f}")

# モジュール 4 の出力。無ければその場で作り直す（AIA 3 MB の取得込みで十数秒）
from workshop import aia_on_eis_grid, BOX

npz = aia_on_eis_grid()
aia = float(np.nanmean(npz["fe18"][BOX["y0"]:BOX["y1"], BOX["x0"]:BOX["x1"]]))
print(f"AIA Fe XVIII（同じ箱の平均） = {aia:.2f} DN/s  （論文 Table 2 = 7.20）")

tresp = np.column_stack([tresp, R])
dn = np.append(dn, aia)
edn = np.append(edn, 0.19 * aia)                          # 論文の 1.40/7.20 = 19%
labels.append("AIA 94 FeXVIII")
print(f"拘束の合計 = {len(dn)} 本")

# %% [markdown]
# **★ 単位のつじつま（3 回踏んだので明示する）**
#
# | | 単位 | demregpy に渡すもの |
# |---|---|---|
# | EIS 輝線 | $I$ = erg cm⁻² s⁻¹ sr⁻¹、$G$ = erg cm³ s⁻¹ | **$G/(4\pi)$** |
# | AIA | DN s⁻¹ pix⁻¹ | **$R(T)$ そのまま**（1/4π と画素立体角込み） |
#
# 返ってくる DEM は **[cm⁻⁵ K⁻¹]**。

# %% [markdown]
# ## 7-3. 解く

# %%
from demregpy import dn2dem

T = 10**logT
tedges = 10 ** np.append(logT - dlt / 2, logT[-1] + dlt / 2)


def solve(**kw):
    """解いて (ビンあたりの EM, chi2_red, モデル強度) を返す。"""
    de, _, _, ch, dr = dn2dem(dn, edn, tresp, logT, tedges, max_iter=30,
                              warn=False, **kw)
    de = np.atleast_1d(np.squeeze(de))
    # ★ demregpy の DEM は [cm^-5/K]。ビンあたりの EM にするには ΔT を掛ける
    return (de * T * np.log(10) * dlt, float(np.squeeze(ch)),
            np.atleast_1d(np.squeeze(dr)))


em, chi2, dn_reg = solve()
print(f"reduced chi2 = {chi2:.2f}")
print(f"EM ピーク    = logT {logT[int(np.nanargmax(em))]:.2f} "
      f"({T[int(np.nanargmax(em))]/1e6:.2f} MK)")
print(f"総 EM        = {em.sum():.2e} cm^-5")
print(f"  → n_e = 1e9 cm^-3 とすると視線長 L = {em.sum()/1e18/1e8:.0f} Mm")
print("  （EM = n_e^2 L。活動領域の視線長として妥当なオーダー）")

# %% [markdown]
# ## 7-4. R = I_obs / I_DEM のパターンを論文と比べる
#
# 論文 Table 2 の R 列がこれ。**個々の線がどれだけ再現できたか**を見る。

# %%
# 比較のため、事前分布に MCMC の解を入れた解（χ² が最も小さい）も出す
m0 = np.loadtxt(MCMC)
norm_mcmc = np.interp(logT, m0[:, 0], m0[:, 2])
norm_mcmc = norm_mcmc / norm_mcmc.max()
em_prior, chi2_prior, dn_prior = solve(dem_norm0=norm_mcmc)
print(f"既定の解 chi2 = {chi2:.2f}   MCMC を事前分布にした解 chi2 = {chi2_prior:.2f}\n")

# 論文 Table 2 の R 列（= I_obs / I_dem）をそのまま転記したもの
PAPER_R = {"Si VII 275.368": 1.10, "Fe IX 188.497": 1.01, "Fe IX 197.862": 0.93,
           "Fe X 184.536": 1.40, "Fe XI 180.401": 0.88, "Fe XI 188.216": 1.11,
           "S X 264.233": 1.00, "Si X 258.375": 0.78, "Fe XII 192.394": 1.01,
           "Fe XII 195.119": 1.04, "Fe XIII 202.044": 1.80, "Fe XIII 203.826": 1.82,
           "Fe XIV 264.787": 0.90, "Fe XIV 270.519": 0.90, "Fe XV 284.160": 0.81,
           "S XIII 256.686": 0.91, "Fe XVI 262.984": 1.04, "Ar XIV 194.396": 1.36,
           "Ca XIV 193.874": 1.31, "Ca XV 200.972": 1.43, "Ca XVI 208.604": 0.75,
           "AIA 94 FeXVIII": 0.98}
print(f"{'line':<20}{'I_obs':>10}{'R (既定)':>10}{'R (MCMC事前)':>13}{'論文 R':>8}")
for lab, o, p1, p2 in zip(labels, dn, dn_reg, dn_prior):
    pr = PAPER_R.get(lab)
    print(f"{lab:<20}{o:10.2f}{o/p1:10.2f}{o/p2:13.2f}"
          f"{(f'{pr:.2f}' if pr else ''):>8}")

# %% [markdown]
# **★ まず、R そのものが事前分布で動くことに注目。**
# 既定の解では R が全体に 1.2–4 と大きいが、
# χ² の小さい（＝よく合っている）解では 1.0–2 に収まる。
# **「どの線が合わないか」自体が、解き方に依存する。**
#
# それでも**パターンは共通**していて、しかも論文と同じ形をしている:
#
# 1. **Fe XIII の 203.826 が突出して外れる**（3–4）。
#    論文でも 1.80 / 1.82、Warren et al. (2011) でも 1.87 / 1.90 と外れている。
#    → **原子データ側の既知の問題**。モジュール 6 で見たとおり、
#      Fe XIII は**実装間の差が最大の線**（密度敏感線）でもある。
#      「コードが違うだけで 3–12% 変わる線」なのだから、
#      原子データ自体の不確かさも同程度以上あると考えるのが自然。
# 2. **Ar XIV / Ca XIV / Ca XV が揃って上がる**（1.6–2.1）。
#    論文も 1.36 / 1.31 / 1.43 と**揃って**上がっている。
#    → 高温側の輝線が系統的に「DEM で説明しきれない」。
#      組成（Ar は高 FIP、Ca は低 FIP）や電離平衡が疑われる。
# 3. Fe XIII を除く 1–2 MK の鉄の線（Fe XI, XII, XIV, XV）は
#    0.8–1.2 に収まる。ここは論文とほぼ同じ水準。
#
# **合わない線も見ておく。** 我々は Fe IX 197.862（1.47）と
# Fe X 184.536（1.45）が高いが、論文は 0.93 と 1.40 で片方しか外れていない。
# 低温側は箱の選び方（moss をどれだけ含むか）に最も敏感な部分で、
# **モジュール 5 で見た Si VII の問題と地続き**である。
#
# **★ ここが重要**: 我々の R のパターンが論文と**同じ形**をしている。
# 絶対強度が 11% 低くても、**DEM 解析の結論は論文と同じ**になる。
#
# （なお PINTofALE の MCMC_DEM で解くと、
#  Ar XIV / Ca XIV / Ca XV が **1.37 / 1.35 / 1.46** と
#  論文の 1.36 / 1.31 / 1.43 をほぼそのまま再現する。
#  正則化とは別の手法で、より論文に近い R が出る。）

# %% [markdown]
# ## 7-5. ★★ 手法・設定で答えがどれだけ動くか
#
# ここがこの章の核心。**同じ観測強度・同じ G(T)** で設定だけを変える。
#
# - `reg_tweak`: 目標 χ²（大きいほど強く平滑化）
# - `dem_norm0`: 初期の重み（＝事前分布）。**MCMC の解を入れてみる**
#
# 比較相手として、論文と同じ **PINTofALE の MCMC_DEM**
# （Kashyap & Drake 1998）の結果を同梱してある（`work/mcmc_dem_result.txt`）。

# %%
def slope(em, t0, t1):
    m = (logT >= t0) & (logT <= t1) & (em > 0)
    return np.polyfit(logT[m], np.log10(em[m]), 1)[0]


# (日本語の説明, 図の凡例（英語）, demregpy に渡す設定)
runs = [
    ("demregpy 既定",              "demregpy default",        dict()),
    ("demregpy reg_tweak=2",       "demregpy reg_tweak=2",    dict(reg_tweak=2.0)),
    ("demregpy gloci=1",           "demregpy gloci=1",        dict(gloci=1)),
    ("demregpy (MCMC を事前分布に)", "demregpy + MCMC prior",   dict(dem_norm0=norm_mcmc)),
]
results = {}
print(f"{'設定':<30}{'chi2':>7}{'EM ピーク':>12}{'alpha':>8}{'beta':>8}")
for tag, en, kw in runs:
    em_, ch, _ = solve(**kw)
    results[en] = em_
    ip = int(np.nanargmax(em_))
    print(f"{tag:<30}{ch:7.2f}{logT[ip]:8.2f} ({T[ip]/1e6:.1f}MK)"
          f"{slope(em_, 6.0, 6.6):+8.2f}{-slope(em_, 6.6, 7.0):+8.2f}")

# MCMC（PINTofALE）: DEM は [cm^-5/logK] なので EM = DEM × ΔlogT
dlt_m = float(np.median(np.diff(m0[:, 0])))
lt_m, em_m = m0[:, 0], m0[:, 2] * dlt_m
ipm = int(np.nanargmax(em_m))


def slope_m(t0, t1):
    m = (lt_m >= t0) & (lt_m <= t1) & (em_m > 0)
    return np.polyfit(lt_m[m], np.log10(em_m[m]), 1)[0]


print(f"{'MCMC_DEM (PINTofALE)':<30}{'—':>7}{lt_m[ipm]:8.2f} "
      f"({10**lt_m[ipm]/1e6:.1f}MK){slope_m(6.0, 6.6):+8.2f}{-slope_m(6.6, 7.0):+8.2f}")
print(f"{'論文 Table 1 region 7':<30}{'—':>7}{'~6.6 (4 MK)':>12}{2.9:+8.2f}{9.0:+8.2f}")

# %%
plt.figure(figsize=(8, 5.5))
for tag, em_ in results.items():
    plt.plot(logT, em_, "-", lw=1.6, label=tag)
plt.plot(lt_m, em_m, "k-", lw=2.6, marker="o", ms=4, label="MCMC_DEM (PINTofALE)")
plt.yscale("log")
plt.xlim(5.6, 7.1)
plt.ylim(1e24, 1e28)
plt.xlabel("log T [K]")
plt.ylabel(r"EM per bin [cm$^{-5}$]")
plt.title("same data, same G(T) — only the method/prior differs")
plt.legend(fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 7-6. 読み取れること
#
# 1. **EM ピークの 4 MK 付近は、どの手法・どの設定でも動かない。**
#    論文の主張の核（活動領域コアの EM は 4 MK に強くピークを持つ）は
#    **手法に依らず頑健**。
# 2. **傾き α は設定で 1.5 – 2.3 と大きく動く。**
#    正則化は既定の重みだとピークを均して α を過小評価する。
# 3. **MCMC の解を事前分布 `dem_norm0` に与えると χ² が改善し、α も MCMC に近づく。**
#    → **答えが事前分布に強く依存している**証拠。
# 4. **どの設定でも χ² が 1 に届かない。**
#    22 本・誤差 22% のデータを、滑らかな DEM で説明しきれていない。
#
# 論文自身もこう書いている:
#
# > It is also clear, however, that the detailed structure of the EM distributions
# > is much more difficult to determine with confidence.
#
# **これを実データで再現したのがこの章の成果。**

# %% [markdown]
# ## 7-7. ★ 単位の罠: per K か per logK か
#
# **コードによって DEM の単位が違う。**
#
# | | DEM の単位 | ビンあたりの EM |
# |---|---|---|
# | PINTofALE (MCMC) | cm⁻⁵ **/ logK** | DEM × ΔlogT |
# | demregpy | cm⁻⁵ **/ K** | DEM × ΔT = DEM × T ln10 ΔlogT |
#
# 揃えないと **10⁷ ずれる**（実際に踏んだ）。
#
# さらに**傾きにも効く**。論文 Eq.(3) の EM 分布 ξ(Te)dTe は
# ξ が [cm⁻⁵/K] なので
#
# $$ \xi\,dT_e = \frac{\rm DEM_{per\,logK}}{T\ln 10}\times(T\ln 10\ d\log T)
#    = {\rm DEM_{per\,logK}}\times d\log T $$
#
# **T は掛からない。** 準備段階でここを間違えて `DEM × T` を EM 分布として
# 傾きを測り、**α を 3.30 と誤って報告した（正しくは 2.30）**。
# T を余計に掛けると傾きがちょうど **+1 ずれる**。

# %%
em_correct = results["demregpy + MCMC prior"]
em_wrong = em_correct * T                     # わざと間違える
print(f"正しい   alpha = {slope(em_correct, 6.0, 6.6):+.2f}")
print(f"T を余計に掛ける = {slope(em_wrong, 6.0, 6.6):+.2f}   ← ちょうど +1 ずれる")
print("→ ピーク温度は動かないので、傾きだけ見ていると気づけない")

# %% [markdown]
# ## 7-8. 傾き α は何を語るか —— 加熱の物理
#
# ここが論文の科学的な主張。
#
# **Parker のナノフレア説**: 磁力線が対流でランダムに揺すられ、ねじれが溜まり、
# 磁気リコネクションで**間欠的に**エネルギーを放出する。
#
# ループの冷却時間 τ_cool と加熱イベントの間隔 τ_heat を比べると:
#
# | | 描像 | EM 分布 |
# |---|---|---|
# | **低頻度加熱** (τ_heat ≫ τ_cool) | ナノフレア。加熱後に十分冷える | 幅広い。**α は緩い**（≲ 2.3） |
# | **高頻度加熱** (τ_heat ≪ τ_cool) | ほぼ定常。冷える暇がない | 鋭くピーク。**α は急** |
#
# 論文が観測した α ≈ 2.9–3.4 は、
# ナノフレア・シミュレーションの上限 2.3 より**急**。
# → 「活動領域コアの高温プラズマは熱平衡に近い」
# → 「加熱は**高頻度**でなければならない」
# → Parker のナノフレア描像（低頻度）への挑戦。
#
# ### ★ 我々の結果が持つ意味（正直に）
#
# 我々の α は **1.5–2.3**。MCMC で 2.30、これは
# **ナノフレア・シミュレーションの上限とちょうど同じ**。
# 論文の 2.9 は明確に上回るが、我々の値は**境界上にある**。
#
# つまり **α の 0.6 の差が科学的な結論を左右する**。そして α は
#
# - 箱の位置（median が 0.83–0.95 で動く範囲でも変わる）
# - DEM の手法・正則化の設定（**1.5–2.3**）
# - 温度ビンの取り方
#
# に敏感。**「α が急だからナノフレアは否定される」と言うには、
# これらの系統誤差を全部押さえる必要がある。**
# これは論文への批判ではなく、**この種の測定の難しさそのもの**。
# 講習会で一番伝えたいのはここ。

# %% [markdown]
# ## 7-9. 実装上の注意（踏むと分からない）
#
# 1. **温度ビンの粗さの要求が、2 手法で正反対。**
#
#    | | MCMC_DEM | demregpy（正則化） |
#    |---|---|---|
#    | 温度ビン | **粗い方が良い**（0.1 dex） | **細かい方が良い**（0.05 dex） |
#    | 理由 | 各ビンが独立パラメータ。劣決定だと解が跳ねる | 平滑化項が劣決定を吸収する。**ビン数 > 拘束数が必須** |
#
#    demregpy に 17 ビン・23 拘束を渡すと GSVD が
#    `ValueError: operands could not be broadcast together` で破綻する。
#    **正則化は劣決定を前提にする手法**だという性格がそのまま出ている。
#
# 2. **MCMC の `dem` は χ² 最小の 1 実現**なので凸凹する。
#    論文が描いているのは MCMC アンサンブルなので、
#    `simdem` から中央値を取る（同梱ファイルはそうしてある）。
#
# 3. **誤差に較正の床を入れないと χ² が発散する**（モジュール 2）。

# %% [markdown]
# ## 7-10. 演習
#
# 1. **AIA Fe XVIII を外して**解き直す。高温側の傾き β がどうなるか。
#    （"toothpaste tube effect" を自分で見る）
# 2. 誤差を 22% → 10% にすると χ² と解はどう変わるか。
#    **誤差を小さく見積もることの危険**を体感する。
# 3. Fe XIII の 2 本を拘束から外して解き直す。
#    他の線の R は改善するか？ **1 つの外れ値が全体をどれだけ引っ張るか。**
# 4. 温度範囲を `logT <= 6.8` に狭めると何が起きるか。
# 5. **EM loci（モジュール 6）の包絡線の下に、解いた EM が収まっているか**確認する。
#    超えていたらそれだけで間違い。
#
# ## まとめ
#
# - DEM 逆問題は **ill-posed**。22 本測っても実質的な自由度は 12 程度
# - だから**必ず追加の仮定が要る**。「解いた」と言うときは仮定もセットで言う
# - **EM ピーク 4 MK は手法に依らず頑健**。**傾き α は手法で 1.5–2.3 と動く**
# - R のパターン（Fe XIII が外れ、Ar XIV/Ca XIV/Ca XV が揃って上がる）は
#   **論文と同じ形**。絶対強度が 11% 低くても結論は変わらない
# - **単位（per K / per logK）を間違えると傾きが +1 ずれる**
#
# 以上で「1 日コース」が完走。
# 発展（モジュール 8–10）では Ca XVII のブレンド分離、較正の効き、
# 他の活動領域へ進む。
