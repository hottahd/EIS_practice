# %% [markdown]
# # モジュール 6: 寄与関数 G(T) と EM loci
#
# **所要時間 40 分**
#
# **このノートで身につくこと**
#
# 1. 寄与関数 G(T) が何でできているか（組成・電離平衡・励起）を理解する
# 2. **輝線ごとに効く温度が違う**ことを数値で見る（これが DEM の武器）
# 3. ★ **1/(4π) の罠**を知る。実装によって G(T) が 12.6 倍ずれる
# 4. **EM loci** を描いて、逆問題を解く前に答えの見当をつける
#
# 前提: モジュール 2（22 輝線の強度）。
#
# ---
#
# **Colab での方針**: CHIANTI のデータベースは 257 MB あり、
# セッションが切れると消える。そこで **事前計算した G(T) をリポジトリに同梱**し、
# 既定ではそれを読む。自分で計算する経路（fiasco）も 6-4 に用意した。

# %%
!pip install -q eispac

# %%
import os
import urllib.request

import numpy as np
import matplotlib.pyplot as plt

GOFNT = "work/gofnt_chianti901.txt"          # 0.1 dex 刻み（MCMC 用）
GOFNT_FINE = "work/gofnt_chianti901_005.txt"  # 0.05 dex 刻み（demregpy 用）
print(open(GOFNT).read().split("# nT")[0])   # ヘッダ = 何を仮定して作ったか

# %% [markdown]
# ## 6-1. G(T) は 3 つの因子でできている
#
# コロナは光学的に薄いので、視線積分は**単なる足し算**:
#
# $$ I_\lambda = \frac{1}{4\pi}\int G_\lambda(T)\, n_e n_H\, ds $$
#
# $$ G_\lambda(T) = \underbrace{A(Z)}_{組成}
#    \times \underbrace{f_{\rm ion}(T)}_{電離平衡}
#    \times \frac{n_H}{n_e}
#    \times \underbrace{\frac{hc}{\lambda}\frac{n_j}{n_{\rm ion}}A_{ji}\frac{1}{n_e}}_{励起（準位占有数）} $$
#
# | 因子 | 何で決まるか | 不確かさ |
# |---|---|---|
# | 組成 $A(Z)$ | コロナ組成 (Feldman 1992) か光球組成か | **FIP 効果で 3-4 倍** |
# | 電離平衡 $f_{\rm ion}(T)$ | 電離・再結合レート | 数十 % |
# | 励起 | 衝突励起断面積、準位占有数 | 数 %〜数十 % |
#
# **G(T) が温度の狭い関数になるのは、主に電離平衡が狭いから。**
# 各イオンは log T で 0.2–0.3 dex の幅でしか存在しない
# （温度が低ければまだ電離しておらず、高ければさらに電離してしまう）。
#
# **★ これは「電離平衡が成り立っている」という仮定の上の話。**
# 急激な加熱・冷却では電離が追いつかない（非平衡電離）。論文もこの仮定に立つ。

# %% [markdown]
# ## 6-2. 事前計算した G(T) を読んで描く

# %%
def read_gofnt(path):
    """09_gofnt.pro / gofnt_fiasco.py が書く形式を読む。"""
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


logT, names, wvl, G = read_gofnt(GOFNT)
print(f"{len(names)} 輝線 x {len(logT)} 温度点  "
      f"(logT {logT[0]:.1f}-{logT[-1]:.1f}, {logT[1]-logT[0]:.2f} dex 刻み)")

print(f"\n{'line':<18}{'G のピーク':>12}{'ピーク温度 logT':>16}{'T [MK]':>9}")
for k in np.argsort([logT[np.argmax(g)] for g in G]):
    j = int(np.argmax(G[k]))
    print(f"{names[k]+' '+f'{wvl[k]:.3f}':<18}{G[k][j]:12.3e}"
          f"{logT[j]:16.2f}{10**logT[j]/1e6:9.2f}")

# %%
fig, ax = plt.subplots(figsize=(8, 5))
cmap = plt.get_cmap("turbo")
tpk = np.array([logT[int(np.argmax(g))] for g in G])
norm = plt.Normalize(tpk.min(), tpk.max())
for k in range(len(names)):
    ax.plot(logT, G[k] / G[k].max(), color=cmap(norm(tpk[k])), lw=1.2)
ax.set_xlim(5.4, 7.2)
ax.set_xlabel("log T [K]")
ax.set_ylabel("G(T) / max")
ax.set_title("contribution functions, normalized (color = peak temperature)")
fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax,
             label="log T at peak")
ax.grid(alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# **これが DEM 解析の武器そのもの。**
# 22 本の輝線が log T = 5.8 から 6.9 まで**少しずつずれた温度**に効いている。
# 各輝線の強度は「その温度あたりのプラズマ量」を測っているので、
# 22 個の測定値から温度分布を復元できる——ように見える。
#
# **★ しかし曲線の幅に注目**。どれも 0.3–0.5 dex の幅があり、大きく重なっている。
# だから「22 個の独立な情報」にはならない。これがモジュール 7 の主題。

# %% [markdown]
# ## 6-3. ★★ 1/(4π) の罠 —— 静かに 12.6 倍ずれる
#
# 準備段階で、**同じ CHIANTI 9.0.1 の ASCII ファイル**を 3 つの実装に
# 読ませて G(T) を比べた（`scripts/idl/09_gofnt.pro`,
# `scripts/gofnt_fiasco.py`, `scripts/gofnt_chiantipy.py`）。
# 比べているのは「DB のバージョン差」ではなく **実装の差だけ**。
#
# | 実装 | 単位の約束 | CHIANTI IDL との比 |
# |---|---|---|
# | CHIANTI IDL `emiss_calc` | hc/λ × N_j × A_ji（**4π で割らない**） | 1.000 |
# | fiasco `contribution_function` | 同上（4π で割らない） | 0.972 – 1.035 |
# | **ChiantiPy `ion.emiss()`** | **sr⁻¹（4π で割ってある）** | **0.077 – 0.082** |
#
# **1/(4π) = 0.0796。** 実測の median は 0.080 でぴったり一致する。
# 4π を掛け戻すと median 1.000（0.968–1.030）で他の 2 者と揃う。
#
# **なぜこれが恐ろしいか**
#
# - **エラーは一切出ない。** 静かに 12.6 倍ずれた DEM が出てくる
# - **全輝線が一律にずれる**ので、線ごとの ratio を見ている限り絶対に気づけない
# - DEM の絶対値は EM = n_e² L に直結するので、
#   **ループの長さや密度の議論が丸ごと狂う**
#
# 気づく方法はただ 1 つ、**オーダーが物理的に妥当かを見る**こと。

# %%
print("オーダーの検算（覚えておく数字）")
print(f"  Fe XII 195.119 の G ピーク = "
      f"{G[int(np.argmin(np.abs(wvl-195.119)))].max():.2e} erg cm^3 s^-1")
print("    → 10^-23 のオーダーなら正しい。10^-24 なら 4π で割られている")
print("  活動領域の EM (視線積分)   = 10^27 - 10^29 cm^-5")
print("  コロナの電子密度 n_e       = 10^9 cm^-3 （活動領域コア）")
print("  → EM = n_e^2 L より L = EM / n_e^2 = 1e28/1e18 = 1e10 cm = 100 Mm")
print("    ループの長さとして妥当。ここが 4 桁ずれたら単位を疑う")

# %% [markdown]
# **★ 実話**: 準備段階で MCMC_DEM の初期 DEM が 7.4e-7 cm⁻⁵ になった。
# `emiss_calc` を n_e で割り忘れていたため（1e9 のずれ）。
# **このオーダー感覚があったから即座に気づけた。**
#
# もう 1 つの実測結果（教材として一級品）:
# 4π を補正すると **2 つの Python 実装は 4 桁一致**（0.999–1.000）し、
# **IDL だけが最大 3.5% 違う**。差が大きいのは
# **Fe XIII 202.044 / 203.826、Fe XIV 264.787、Si X 258.375** —— すべて**密度敏感線**。
# 差がほぼ無いのは Fe IX、Ca XIV、Ca XVI、Fe XV —— **基底準位からの共鳴線**。
#
# → 差の正体は**準安定準位の占有数の扱い**（陽子励起、励起準位への電離・再結合）。
# → **Fe XIII が DEM で唯一大きく外れる**（我々 1.3–2.8、論文 1.80/1.90）ことの
#   独立な説明になっている。コードが違うだけで 3–12% 変わる線なのだから、
#   原子データ自体の不確かさも同程度以上あると考えるのが自然。

# %% [markdown]
# ## 6-4. 自分で G(T) を作る（fiasco、オプション）
#
# **fiasco** は sunpy 系の CHIANTI インターフェース。
# CHIANTI DB を持っているなら、そのまま G(T) を計算できる。
#
# | | fiasco 0.8.2 | ChiantiPy 0.16.0 |
# |---|---|---|
# | CHIANTI IDL との一致 | ✅ そのまま 3% 以内 | ⚠ 4π の補正が要る |
# | numpy 2.x | ✅ | ⚠ 単一温度で落ちる |
# | バッチ実行 | ✅ | ⚠ `chiantirc` が無いと対話を要求（Colab で必ず踏む） |
# | DB の入手 | `download_dbase()` で自動 | 手動 |
#
# → **講習会の本線は fiasco**。ChiantiPy は「IDL から来た人向けの補足」。
#
# **DB のサイズ（実測）**: v9.0.1 = **257 MB**、v10.1 = 1058 MB、
# v11.0.2（fiasco の既定）= 579 MB。
# **v9.0.1 が最も軽く、しかも論文（CHIANTI 7）に最も近い。**

# %%
USE_FIASCO = os.path.exists(os.path.expanduser("~/.fiasco/fiascorc"))
print("fiasco の設定がある" if USE_FIASCO else
      "fiasco の DB が無いので事前計算のファイルを使う（それで十分）")

if USE_FIASCO:
    import astropy.units as u
    import fiasco

    T = 10**logT * u.K
    ion = fiasco.Ion("Fe 12", T, abundance="sun_coronal_1992_feldman")
    cf = ion.contribution_function(1e9 * u.cm**-3)     # (nT, 1, n_transitions)
    k = int(np.argmin(np.abs(ion.transitions.wavelength[~ion.transitions.is_twophoton]
                            .to_value("Angstrom") - 195.119)))
    g_fiasco = cf[:, 0, k].to_value("erg cm3 / s") * 0.83   # n_H/n_e = 0.83
    g_ref = G[int(np.argmin(np.abs(wvl - 195.119)))]
    j = int(np.argmax(g_ref))
    print(f"\nFe XII 195.119 のピーク値")
    print(f"  fiasco       : {g_fiasco[int(np.argmax(g_fiasco))]:.3e}")
    print(f"  同梱ファイル : {g_ref[j]:.3e}   （CHIANTI IDL 由来）")
    print(f"  比           : {g_fiasco[int(np.argmax(g_fiasco))]/g_ref[j]:.3f}")
else:
    print("\n自分で作るなら:")
    print("  from fiasco.util import download_dbase")
    print("  download_dbase('http://download.chiantidatabase.org/"
          "CHIANTI_9.0.1_database.tar.gz', '/content/chianti')")
    print("  （257 MB。scripts/gofnt_fiasco.py が 22 輝線ぶんを作る）")

# %% [markdown]
# **22 輝線すべてで fiasco と CHIANTI IDL は median 1.000、最大 3.5% 差、
# 形成温度は完全一致**することを確認済み（`scripts/gofnt_fiasco.py`）。
# → **講習会は Python だけで寄与関数を出せる。**

# %% [markdown]
# ## 6-5. ★ EM loci —— DEM を解く前に必ず描く図
#
# 輝線 λ の観測強度 $I_\lambda$ に対して、
# 「**もし視線上のプラズマが全部ちょうど温度 T にあったら**、必要な EM はいくらか」:
#
# $$ {\rm EM}_{\rm loci,\lambda}(T) = \frac{4\pi I_\lambda}{G_\lambda(T)} $$
#
# - これは各温度における **EM の上限**。真の EM は必ずこの曲線より**下**にある
# - 全輝線の曲線を重ねると、**下側の包絡線**が DEM の目安になる
# - **等温プラズマなら全曲線が 1 点で交わる。** 交わらなければ多温度

# %%
# モジュール 2 が書いた強度を読む（無ければその場で作る）
import csv
import sys
sys.path.insert(0, "scripts")

if not os.path.exists("work/box_intensities.csv"):
    print("work/box_intensities.csv が無いので 22 輝線をフィットする（約 5 秒）")
    import eispac
    from lines_warren2012 import LINES, pick_component
    from fit_box_spectra import average_spectrum
    BOX = dict(y0=244, y1=274, x0=32, x1=40)
    rows = []
    for ion, w, tname, ip, sp in LINES:
        ww, I, s, _ = average_spectrum("data/eis/eis_20110702_030712.data.h5", w, **BOX)
        t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
        comp, _ = pick_component(t, w)
        f = eispac.fit_spectra(I, t, wave=ww, errs=s, ncpu=1, ignore_warnings=True)
        rows.append((ion, w, float(np.atleast_1d(f.fit["int"][..., comp]).ravel()[0]), ip,
                     float(np.atleast_1d(f.fit["int"][..., comp]).ravel()[0]) / ip))
    os.makedirs("work", exist_ok=True)
    with open("work/box_intensities.csv", "w", newline="") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["ion", "wvl", "I_fit", "I_paper", "ratio"])
        w_.writerows(rows)

iobs = np.zeros(len(names))
for r in csv.DictReader(open("work/box_intensities.csv")):
    k = int(np.argmin(np.abs(wvl - float(r["wvl"]))))
    if abs(wvl[k] - float(r["wvl"])) < 0.01:
        iobs[k] = float(r["I_fit"])

# ★ Ca XVII 192.858 はブレンドしたままの値（論文の 5 倍）なので EM loci から外す。
#   モジュール 8 で分離したら戻す。
iobs[np.argmin(np.abs(wvl - 192.858))] = 0.0
print(f"EM loci に使う輝線: {int((iobs > 0).sum())} 本")

# %%
logTf, namesf, wvlf, Gf = read_gofnt(GOFNT_FINE)      # 細かい格子の方が見やすい
iobsf = np.zeros(len(namesf))
for r in csv.DictReader(open("work/box_intensities.csv")):
    k = int(np.argmin(np.abs(wvlf - float(r["wvl"]))))
    if abs(wvlf[k] - float(r["wvl"])) < 0.01:
        iobsf[k] = float(r["I_fit"])
iobsf[np.argmin(np.abs(wvlf - 192.858))] = 0.0

fig, ax = plt.subplots(figsize=(8.5, 5.8))
cmap = plt.get_cmap("turbo")
ok = np.where(iobsf > 0)[0]
tpk = np.array([logTf[int(np.argmax(Gf[k]))] for k in ok])
norm = plt.Normalize(tpk.min(), tpk.max())

env = np.full(len(logTf), np.inf)
for k, tp in zip(ok, tpk):
    g = Gf[k]
    m = g > g.max() * 1e-3                     # G が十分ある温度だけ描く
    loci = 4 * np.pi * iobsf[k] / np.where(g > 0, g, np.nan)
    ax.plot(logTf[m], loci[m], color=cmap(norm(tp)), lw=1.3, alpha=0.9)
    j = int(np.nanargmin(np.where(m, loci, np.inf)))
    ax.annotate(f"{namesf[k]} {wvlf[k]:.1f}", (logTf[j], loci[j]), fontsize=6.5,
                color=cmap(norm(tp)), xytext=(2, 2), textcoords="offset points")
    env = np.minimum(env, np.where(m, loci, np.inf))

ax.plot(logTf, env, "k--", lw=1.8, label="lower envelope = upper limit on EM")
ax.set_yscale("log")
ax.set_xlim(5.4, 7.2)
ax.set_ylim(1e25, 1e31)
ax.set_xlabel("log T [K]")
ax.set_ylabel(r"EM$_{\rm loci} = 4\pi I / G(T)$   [cm$^{-5}$]")
ax.set_title("EM loci: EM required if ALL the plasma were at temperature T")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax,
             label="log T at peak of G(T)")
fig.tight_layout()
plt.show()

print(f"{'logT':>6} {'T [MK]':>8} {'包絡線 = EM の上限 [cm^-5]':>28}")
for t in (5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0):
    j = int(np.argmin(np.abs(logTf - t)))
    print(f"{t:6.1f} {10**t/1e6:8.2f} {env[j]:28.2e}")

# %% [markdown]
# **読み方**
#
# 1. **曲線は 1 点で交わらない** → **等温ではない**。多温度のプラズマがある。
#    これが「DEM を解く」動機そのもの。
# 2. **包絡線は各温度における EM の上限**。真の EM 分布は必ずこの下にある。
# 3. **log T ≈ 5.8 の深い谷**（Si VII が作っている）。
#    「もし全部 0.6 MK なら EM は 7.6e25 しか要らない」= **低温プラズマは非常に少ない**。
#    inter-moss（ループ上部）を選んだ効果がここに出ている。
# 4. 包絡線は低温側から **6.8 まで単調に上がる**（7.6e25 → 7.4e27）。
#    上限が緩いほど、その温度には EM があってよい。
#    **多数の曲線の底が log T 6.1–6.6 に集まっている**のが本体の在りか。
#    モジュール 7 で解く DEM のピーク（log T 6.6）は、
#    そこでの上限 7.3×10²⁷ の下にちゃんと収まる。**これが検算になる。**
# 5. **log T = 7.0 で包絡線が 1.3×10²⁹ に跳ね上がる**。
#    そこに効く輝線がもう無い = **上限が事実上つかない**ということ。
#    → **7 MK 以上はこのデータだけでは決まらない**。だから AIA Fe XVIII を足す
#      （入れないと DEM が高温側に漏れる。PINTofALE の文書が
#       "toothpaste tube effect" と呼ぶ現象）。
# 6. EM のオーダーは **10²⁷–10²⁸ cm⁻⁵**。活動領域として妥当（6-3 の検算どおり）。
#
# **★ この図を先に描いておくと、DEM の解が変になったときの検算に使える。**
# 逆問題の答えが包絡線を超えていたら、それだけで間違い。

# %% [markdown]
# ## 6-6. 演習
#
# 1. `iobs` を全部 1.3 倍して EM loci を描き直す（較正が 30% 違ったら、の想定）。
#    **形は変わるか、それとも上下に平行移動するだけか？**
#    → 「較正誤差では説明できないパターン」の意味が体で分かる。
# 2. Ca XVII をブレンドしたままの値（論文の 5 倍）で EM loci に入れてみる。
#    包絡線がどう壊れるか。**1 本の誤った線が高温側の結論をどう変えるか。**
# 3. G(T) のピーク温度が最も近い 2 本（Fe XI 180.401 と Fe XI 188.216）の
#    EM loci 曲線を比べる。**同じイオンなら重なるはず**。重ならなければ何が違うのか。
# 4. 組成をコロナ組成から光球組成に変えたら G(T) はどう変わるか（発展、fiasco が要る）。
#    Fe と Ca は低 FIP 元素、S と Ar は高 FIP 元素。
#
# ## まとめ
#
# - G(T) = 組成 × 電離平衡 × 励起。**温度が狭いのは電離平衡のおかげ**
# - 22 輝線が logT 5.8–6.9 に少しずつずれて効く。**でも幅広く重なっている**
# - **1/(4π) の約束は実装ごとに違う**。エラーは出ない。オーダーで検算する
# - **EM loci は逆問題を解く前に必ず描く**。答えの見当と検算に使える
#
# 次（モジュール 7）で、いよいよ **DEM の逆問題**を解く。
