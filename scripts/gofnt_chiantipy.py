"""Warren+2012 の 22 輝線の寄与関数 G(T) を **ChiantiPy** で作り、
CHIANTI IDL / fiasco の結果と突き合わせる（3 実装の比較）。

3 者に **同じ CHIANTI 9.0.1 の ASCII ファイル**を読ませるのが肝:
  - CHIANTI IDL : `!xuvtop = /opt/ssw/packages/chianti/dbase`
  - fiasco      : `~/.fiasco/fiascorc` の ascii_dbase_root
  - ChiantiPy   : 環境変数 `XUVTOP`
→ 比較しているのは「DB のバージョン差」ではなく **実装の差だけ**。

定義（3 者で厳密に揃える）:

    G(T) = Ab(Z) * ioneq(Z,ion,T) * (N_H/N_e) * eps(T) / N_e   [erg cm^3 s^-1]

    eps = hc/λ × N_j × A_ji（上準位のイオンに対する存在比 × 遷移確率）
    ChiantiPy では `ion.emiss()` → `ion.Emiss['emiss']`、形は (nWvl, nT)。

★★ 単位の落とし穴（この検証で見つけた最重要の知見）
    **ChiantiPy の `emiss` は最初から sr^-1（= 4π で割ってある）。**
    CHIANTI IDL の `emiss_calc` と fiasco の `contribution_function` は
    4π で割っていない。そのため素で比べると ChiantiPy だけ **1/(4π) = 0.0796 倍**
    になる。実測 0.077-0.082、4π を掛け戻すと median 1.000（0.968-1.030）で
    他の 2 者と一致する。

    エラーは一切出ず、**全輝線が一律にずれる**ので線比を見ている限り気づけない。
    DEM の絶対値は EM = n_e^2 L に直結するので、ここを外すとループ長や密度の
    議論が丸ごと 12.6 倍狂う。

準備:
  export XUVTOP=/opt/ssw/packages/chianti/dbase
  ~/.config/chiantirc に既定を書いておく（無いと import 時に対話で聞いてくる）:
      [chianti]
      abundfile = sun_coronal_1992_feldman
      ioneqfile = chianti
      wavelength = angstrom
      flux = energy
      gui = False

    python scripts/gofnt_chiantipy.py
"""
import os
import sys

import numpy as np

OUT = "work/gofnt_chiantipy.txt"
REF_IDL = "work/gofnt_chianti901.txt"
REF_FIASCO = "work/gofnt_fiasco.txt"

LOGNE = 9.0
NHNE = 0.83

# (表示名, ChiantiPy のイオン名, 波長[A])
LINES = [
    ("Si VII", "si_7", 275.368), ("Fe IX", "fe_9", 188.497),
    ("Fe IX", "fe_9", 197.862), ("Fe X", "fe_10", 184.536),
    ("Fe XI", "fe_11", 180.401), ("Fe XI", "fe_11", 188.216),
    ("S X", "s_10", 264.233), ("Si X", "si_10", 258.375),
    ("Fe XII", "fe_12", 192.394), ("Fe XII", "fe_12", 195.119),
    ("Fe XIII", "fe_13", 202.044), ("Fe XIII", "fe_13", 203.826),
    ("Fe XIV", "fe_14", 264.787), ("Fe XIV", "fe_14", 270.519),
    ("Fe XV", "fe_15", 284.160), ("S XIII", "s_13", 256.686),
    ("Fe XVI", "fe_16", 262.984), ("Ar XIV", "ar_14", 194.396),
    ("Ca XIV", "ca_14", 193.874), ("Ca XV", "ca_15", 200.972),
    ("Ca XVI", "ca_16", 208.604), ("Ca XVII", "ca_17", 192.858),
]


def read_gofnt(path):
    """IDL / fiasco が書き出した G(T) テキストを読む。"""
    with open(path) as f:
        lines = f.readlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# nT nline"))
    nT, nline = (int(v) for v in lines[i + 1].split())
    k = i + 2

    def skip_to(tag):
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

    skip_to("# logT")
    logT = take(nT)
    skip_to("# ion")
    wtar = []
    for _ in range(nline):
        wtar.append(float(lines[k].split()[-2 if path == REF_IDL else -1]))
        k += 1
    skip_to("# G(T)")
    G = np.array([take(nT) for _ in range(nline)])
    return logT, np.array(wtar), G


def main():
    if "XUVTOP" not in os.environ:
        sys.exit("XUVTOP を設定すること: export XUVTOP=/opt/ssw/packages/chianti/dbase")
    print("XUVTOP =", os.environ["XUVTOP"])

    logT, wtar_idl, G_idl = read_gofnt(REF_IDL)
    have_fiasco = os.path.exists(REF_FIASCO)
    if have_fiasco:
        _, wtar_fi, G_fi = read_gofnt(REF_FIASCO)

    import ChiantiPy.core as ch
    import ChiantiPy.tools.data as chdata

    T = 10 ** logT
    ne = 10 ** LOGNE
    print(f"nT={len(logT)}  logT {logT[0]:.2f}-{logT[-1]:.2f}  logNe={LOGNE}")
    print(f"abundfile = {chdata.Defaults['abundfile']}, "
          f"ioneqfile = {chdata.Defaults['ioneqfile']}")

    hdr = f"\n{'line':<10} {'wvl':>9} {'G_ChiantiPy':>13} {'G_IDL':>12} {'比':>7}"
    if have_fiasco:
        hdr += f" {'G_fiasco':>12} {'比':>7}"
    print(hdr + f" {'logTmax':>8}")

    cache, rows, r_idl, r_fi = {}, [], [], []
    for label, ionname, wvl in LINES:
        if ionname not in cache:
            ion = ch.ion(ionname, temperature=T, eDensity=ne)
            ion.emiss()
            cache[ionname] = ion
        ion = cache[ionname]
        w = np.array(ion.Emiss["wvl"])
        e = np.array(ion.Emiss["emiss"])          # (nWvl, nT)

        # ±0.03 A の中で放射率が最大の遷移（IDL / fiasco と同じ規則）
        m = np.abs(w - wvl) < 0.03
        if not m.any():
            print(f"{label:<10} {wvl:9.3f}  遷移が見つからない")
            continue
        idx = np.where(m)[0]
        j = idx[int(np.argmax([e[i].max() for i in idx]))]

        # ioneq と組成。ChiantiPy は ion 生成時に IoneqOne / Abundance を持つ
        f_ioneq = np.asarray(ion.IoneqOne, dtype=float)
        ab = float(ion.Abundance)
        # ★ ChiantiPy の emiss は sr^-1 なので 4π を掛けて他の 2 者に単位を揃える
        g = ab * f_ioneq * NHNE * e[j] / ne * (4.0 * np.pi)

        k = int(np.argmin(np.abs(wtar_idl - wvl)))
        ri = g.max() / G_idl[k].max()
        r_idl.append(ri)
        line = (f"{label:<10} {wvl:9.3f} {g.max():13.4e} {G_idl[k].max():12.4e} {ri:7.3f}")
        if have_fiasco:
            kf = int(np.argmin(np.abs(wtar_fi - wvl)))
            rf = g.max() / G_fi[kf].max()
            r_fi.append(rf)
            line += f" {G_fi[kf].max():12.4e} {rf:7.3f}"
        line += f" {logT[int(np.argmax(g))]:8.2f}"
        print(line)
        rows.append((label, wvl, g))

    def summarise(name, a):
        a = np.array(a)
        print(f"{name}: median={np.median(a):.3f} min={a.min():.3f} max={a.max():.3f} "
              f"15%以内 {int(((a >= 0.85) & (a <= 1.15)).sum())}/{len(a)}")

    print()
    summarise("G_ChiantiPy / G_IDL   ", r_idl)
    if have_fiasco:
        summarise("G_ChiantiPy / G_fiasco", r_fi)

    with open(OUT, "w") as f:
        f.write(f"# G(T) from ChiantiPy, XUVTOP={os.environ['XUVTOP']}\n")
        f.write(f"# abund={chdata.Defaults['abundfile']} "
                f"ioneq={chdata.Defaults['ioneqfile']} logNe={LOGNE} N_H/N_e={NHNE}\n")
        f.write("# nT nline\n")
        f.write(f"{len(logT)} {len(rows)}\n# logT\n")
        f.write(" ".join(f"{x:.4f}" for x in logT) + "\n# ion wvl\n")
        for label, wvl, _ in rows:
            f.write(f"{label} {wvl:.4f}\n")
        f.write("# G(T)\n")
        for _, _, g in rows:
            f.write(" ".join(f"{x:.6e}" for x in g) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
