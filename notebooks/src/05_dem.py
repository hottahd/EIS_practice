# %% [markdown]
# # 第 5 章: 温度分布を出す
#
# **30 分**
#
# 第 1 章で見た「輝線ごとにまったく違う絵」の正体を、数字にします。
#
# 視線上には色々な温度のプラズマが混ざっています。
# **どの温度がどれだけあるか**を求めるのが **DEM（emission measure 分布）解析**です。

# %% [markdown]
# ## 5-1. 輝線は温度計になる
#
# 光学的に薄いコロナでは、輝線強度は視線上の足し算になります:
#
# $$ I_\lambda = \frac{1}{4\pi}\int G_\lambda(T)\, n_e n_H\, ds $$
#
# $G_\lambda(T)$ が**寄与関数**で、「その輝線がどの温度で光るか」を表します。
# 中身は **組成 × 電離平衡 × 励起**で、原子データベース（CHIANTI）から計算します。
#
# **温度の幅が狭いのは電離平衡のおかげ**です。各イオンは log T で 0.2–0.3 dex
# の範囲でしか存在できません（低温では電離しておらず、高温ではさらに電離する）。
#
# CHIANTI の計算には時間がかかるので、**事前に計算したものを同梱**しています。

# %%
import numpy as np
import matplotlib.pyplot as plt

from workshop import box_intensities


def read_gofnt(path):
    """事前計算した G(T) を読む。"""
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
print(f"{len(names)} 輝線 × {len(logT)} 温度点  "
      f"(logT {logT[0]:.1f}–{logT[-1]:.1f}, {logT[1]-logT[0]:.2f} dex 刻み)")

fig, ax = plt.subplots(figsize=(8, 4.5))
cmap = plt.get_cmap("turbo")
tpk = np.array([logT[int(np.argmax(g))] for g in G])
norm = plt.Normalize(tpk.min(), tpk.max())
for k in range(len(names)):
    ax.plot(logT, G[k] / G[k].max(), color=cmap(norm(tpk[k])), lw=1.2)
ax.set_xlim(5.4, 7.2)
ax.set_xlabel("log T [K]")
ax.set_ylabel("G(T) / max")
ax.set_title("contribution functions (color = peak temperature)")
fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax, label="log T at peak")
ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# 22 本の輝線が log T = 5.8 から 6.9 まで**少しずつずれた温度**に効いています。
# これが温度分布を測る武器です。
#
# ただし**曲線は幅広く重なっています**。だから 22 本測っても
# 22 個ぶんの独立な情報にはなりません（後述）。

# %% [markdown]
# ## 5-2. EM loci —— 解く前に答えの見当をつける
#
# 輝線ごとに「**もし視線上のプラズマが全部ちょうど温度 T にあったら、
# EM はいくら必要か**」を計算します:
#
# $$ {\rm EM}_{\rm loci}(T) = \frac{4\pi I_\lambda}{G_\lambda(T)} $$
#
# これは各温度における **EM の上限**で、真の値は必ずこの下にあります。
# **等温プラズマなら全部の曲線が 1 点で交わります。**

# %%
rows = box_intensities()        # 第 2 章と同じ箱の 22 輝線（無ければその場で作る）

iobs = np.zeros(len(names))
for ion, w, i_fit, i_paper, ratio in rows:
    k = int(np.argmin(np.abs(wvl - w)))
    if abs(wvl[k] - w) < 0.01:
        iobs[k] = i_fit
iobs[int(np.argmin(np.abs(wvl - 192.858)))] = 0.0     # Ca XVII はブレンド（付録 G）

ok = np.where(iobs > 0)[0]
fig, ax = plt.subplots(figsize=(8.5, 5.5))
tpk_ok = np.array([logT[int(np.argmax(G[k]))] for k in ok])
norm = plt.Normalize(tpk_ok.min(), tpk_ok.max())
env = np.full(len(logT), np.inf)
for k, tp in zip(ok, tpk_ok):
    g = G[k]
    m = g > g.max() * 1e-3
    loci = 4 * np.pi * iobs[k] / np.where(g > 0, g, np.nan)
    ax.plot(logT[m], loci[m], color=cmap(norm(tp)), lw=1.3, alpha=0.9)
    j = int(np.nanargmin(np.where(m, loci, np.inf)))
    ax.annotate(f"{names[k]} {wvl[k]:.1f}", (logT[j], loci[j]), fontsize=6.5,
                color=cmap(norm(tp)), xytext=(2, 2), textcoords="offset points")
    env = np.minimum(env, np.where(m, loci, np.inf))
ax.plot(logT, env, "k--", lw=1.8, label="lower envelope = upper limit")
ax.set_yscale("log")
ax.set_xlim(5.4, 7.2)
ax.set_ylim(1e25, 1e31)
ax.set_xlabel("log T [K]")
ax.set_ylabel(r"EM$_{\rm loci} = 4\pi I / G(T)$  [cm$^{-5}$]")
ax.set_title("EM loci: EM required if ALL the plasma were at temperature T")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax, label="log T at peak")
fig.tight_layout()
plt.show()

# %% [markdown]
# **1 点で交わっていません。→ 等温ではない。多温度です。**
#
# - 低温側（Si VII, Fe IX）の曲線が高いところにある = **冷たいプラズマは少ない**
# - 曲線の底が log T 6.1–6.6 に集まっている = **本体はこの温度域**
# - log T > 6.9 では拘束する輝線が無く、上限が跳ね上がる
#   → 本当は AIA 94 Å の Fe XVIII（7 MK）を足したいところ（付録 F）

# %% [markdown]
# ## 5-3. DEM を解く
#
# $I_\lambda = \frac{1}{4\pi}\int G_\lambda(T)\,\xi(T)\,dT$ を $\xi(T)$ について解きます。
# ここでは正則化インバージョン（`demregpy`）を使います。
#
# 誤差は論文と同じ **22%**（統計誤差 + 較正の不確かさ）を使います（付録 B）。

# %%
from demregpy import dn2dem

keep = (logT >= 5.5) & (logT <= 7.1)
lt, Gk = logT[keep], G[:, keep]
dlt = lt[1] - lt[0]

sel = iobs > 0
dn = iobs[sel]
edn = 0.22 * dn                       # 較正の系統誤差（付録 B）
tresp = (Gk[sel] / (4 * np.pi)).T     # I = (1/4π)∫G ξ dT なので G/(4π)
tedges = 10 ** np.append(lt - dlt / 2, lt[-1] + dlt / 2)

dem, edem, elogt, chisq, dn_reg = dn2dem(dn, edn, tresp, lt, tedges,
                                         max_iter=30, warn=False)
dem = np.atleast_1d(np.squeeze(dem))
T = 10**lt
em = dem * T * np.log(10) * dlt       # DEM [cm^-5/K] → ビンあたりの EM [cm^-5]

ipk = int(np.nanargmax(em))
print(f"使った輝線     : {int(sel.sum())} 本")
print(f"reduced chi2   : {float(np.squeeze(chisq)):.2f}")
print(f"EM のピーク    : logT {lt[ipk]:.2f}  = {T[ipk]/1e6:.1f} MK")
print(f"総 EM          : {em.sum():.2e} cm^-5")
print(f"  → n_e = 1e9 cm^-3 なら視線長 {em.sum()/1e18/1e8:.0f} Mm（妥当なオーダー）")

# %%
plt.figure(figsize=(7.5, 5))
plt.plot(lt, em, "o-", lw=2, ms=4)
plt.axvline(lt[ipk], color="r", ls="--", lw=1,
            label=f"peak: {T[ipk]/1e6:.1f} MK")
plt.yscale("log")
plt.xlim(5.6, 7.1)
plt.ylim(1e24, 1e28)
plt.xlabel("log T [K]")
plt.ylabel(r"EM per bin [cm$^{-5}$]")
plt.title("emission measure distribution (inter-moss region, NOAA 1243)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 5-4. 何が分かったか
#
# **活動領域コアの emission measure は 4 MK 付近に強くピークを持つ。**
# 論文（Warren+2012）の主要な結論と同じです。
#
# これが効いてくるのは**コロナ加熱の議論**です。
#
# - 加熱が**まれ**（ナノフレア的）なら、ループは加熱のあいだに冷える時間があり、
#   色々な温度のプラズマが視線上に溜まる → **EM 分布は幅広くなる**
# - 加熱が**頻繁**なら、ループは冷える暇がなく高温に保たれる
#   → **EM 分布は鋭くピークを持つ**
#
# 観測されたピークの鋭さは、**加熱が高頻度**であることを示唆します。
# Parker のナノフレア描像（低頻度加熱）への挑戦になっている、というのが論文の主張です。
#
# **ただし、この「鋭さ」は解き方に依存します。** ピーク温度は動きませんが、
# 傾きは手法や設定で変わります（付録 C）。
# DEM は逆問題で、**必ず何らかの仮定が入る**ことは覚えておいてください。

# %% [markdown]
# ## 5-5. 演習
#
# 1. **輝線を減らして解いてみる。** 高温側（Ca XIV–XVI）を落とすと
#    ピークはどうなるか。`sel` を書き換えて確かめる
# 2. **誤差を 22% から 10% に変える。** chi2 と解の形はどう変わるか
# 3. EM loci の**包絡線の下に**、解いた EM がちゃんと収まっているか確認する
#    （超えていたらそれだけで間違い）

# %%
# 演習 1（TODO を埋める）
#
# sel2 = sel.copy()
# for bad in ["Ca XIV", "Ca XV", "Ca XVI"]:
#     sel2[[i for i, n in enumerate(names) if n == bad]] = ____
# dem2, _, _, chi2_2, _ = dn2dem(iobs[sel2], 0.22*iobs[sel2],
#                                (Gk[sel2]/(4*np.pi)).T, lt, tedges,
#                                max_iter=30, warn=False)
# em2 = np.atleast_1d(np.squeeze(dem2)) * T * np.log(10) * dlt
# plt.plot(lt, em, "o-", label="all lines")
# plt.plot(lt, em2, "s-", label="without Ca")
# plt.yscale("log"); plt.legend(); plt.show()

# %% [markdown]
# ## 今日のまとめ
#
# | 章 | 出したもの | Solar-C で |
# |---|---|---|
# | 1 | 輝線ごとの強度マップ | 同じ輝線を 0.4″ で |
# | 2 | **輝線強度**（ガウシアンの面積） | そのまま使う |
# | 3 | **ドップラー速度**（ゼロ点は自分で決める） | そのまま使う |
# | 4 | **非熱的速度**（装置幅を引く） | そのまま使う |
# | 5 | **温度分布**（4 MK にピーク） | 温度被覆が広がってもっと良く決まる |
#
# **今日の手順は、装置が変わってもそのまま使えます。**
# 2028 年に EUVST のデータが降りてきたら、同じことをすればいい。
#
# 付録には、今日触れなかった話（欠損値、誤差、DEM の信頼度、較正、
# ブレンド分離など）をまとめてあります。自分の研究で使うときに読んでください。
