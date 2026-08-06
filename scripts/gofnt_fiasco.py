"""Warren+2012 の 22 輝線の寄与関数 G(T) を **fiasco (Python)** で作り、
CHIANTI IDL (`scripts/idl/09_gofnt.pro`) の結果と突き合わせる。

★リストの項目 5。講習会は Python/Colab を本線にするので、
SSW を使わずに G(T) を出せることを確認しておく必要がある。

**同じ原子データを読ませる**のが肝。fiasco の ascii_dbase_root を
SSW 同梱の CHIANTI 9.0.1 (`/opt/ssw/packages/chianti/dbase`) に向けてあるので
（`~/.fiasco/fiascorc`）、この比較は「DB のバージョン差」ではなく
**2 つのコードの G(T) 計算の差だけ**を見るテストになる。

定義（IDL 側と厳密に合わせる。scripts/idl/09_gofnt.pro のコメント参照）:

    G(T) = Ab(Z) * ioneq(Z,ion,T) * (N_H/N_e) * eps(T) / N_e     [erg cm^3 s^-1]

    eps は「hc/λ × N_j × A_ji」（上準位 j のイオンに対する存在比 × 遷移確率）。
    fiasco の Ion.contribution_function() は
        G = (hc/λ) * (N_j/N_ion) * A_ji * (N_ion/N_elem) * (N_elem/N_H) * (1/N_e)
    を返す（`couple_density_to_temperature=False` なら (T, ne, λ) の 3 次元）。
    N_H/N_e は含まれないので、こちらで 0.83 を掛ける。

    python scripts/gofnt_fiasco.py
"""
import sys

import numpy as np
import astropy.units as u

OUT = "work/gofnt_fiasco.txt"
REF = "work/gofnt_chianti901.txt"

LOGNE = 9.0
NHNE = 0.83
ABUND = "sun_coronal_1992_feldman"
IONEQ = "chianti"

# (表示名, fiasco のイオン名, 波長[A])  ... IDL 側と同じ 22 本
LINES = [
    ("Si VII", "Si 7", 275.368), ("Fe IX", "Fe 9", 188.497),
    ("Fe IX", "Fe 9", 197.862), ("Fe X", "Fe 10", 184.536),
    ("Fe XI", "Fe 11", 180.401), ("Fe XI", "Fe 11", 188.216),
    ("S X", "S 10", 264.233), ("Si X", "Si 10", 258.375),
    ("Fe XII", "Fe 12", 192.394), ("Fe XII", "Fe 12", 195.119),
    ("Fe XIII", "Fe 13", 202.044), ("Fe XIII", "Fe 13", 203.826),
    ("Fe XIV", "Fe 14", 264.787), ("Fe XIV", "Fe 14", 270.519),
    ("Fe XV", "Fe 15", 284.160), ("S XIII", "S 13", 256.686),
    ("Fe XVI", "Fe 16", 262.984), ("Ar XIV", "Ar 14", 194.396),
    ("Ca XIV", "Ca 14", 193.874), ("Ca XV", "Ca 15", 200.972),
    ("Ca XVI", "Ca 16", 208.604), ("Ca XVII", "Ca 17", 192.858),
]


def read_idl(path):
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
    names, wtar = [], []
    for _ in range(nline):
        p = lines[k].split()
        names.append(p[0] + " " + p[1])
        wtar.append(float(p[2]))
        k += 1
    skip_to("# G(T)")
    G = np.array([take(nT) for _ in range(nline)])
    return logT, names, np.array(wtar), G


def main():
    import fiasco

    logT_idl, names_idl, wtar_idl, G_idl = read_idl(REF)
    T = (10 ** logT_idl) * u.K
    ne = (10 ** LOGNE) * u.cm ** -3
    print(f"IDL 側: nT={len(logT_idl)}  logT {logT_idl[0]:.2f}-{logT_idl[-1]:.2f}")

    print(f"\n{'line':<10} {'wvl':>9} {'G_fiasco':>12} {'G_IDL':>12} {'比':>7} "
          f"{'logT(fiasco)':>13} {'logT(IDL)':>10}")
    rows, ratios = [], []
    cache = {}
    for label, ionname, wvl in LINES:
        if ionname not in cache:
            cache[ionname] = fiasco.Ion(ionname, T, abundance=ABUND, ioneq_filename=IONEQ, ask_before=False)
        ion = cache[ionname]
        try:
            g = ion.contribution_function(ne)      # (T, ne, transitions)
        except Exception as e:
            print(f"{label:<10} {wvl:9.3f}  FAILED: {e}")
            continue
        lam = ion.transitions.wavelength[~ion.transitions.is_twophoton]
        g = g[:, 0, :]                              # 密度 1 点なので落とす
        # ±0.03 A の中で放射率が最大の遷移を選ぶ（IDL 側と同じ規則）
        m = np.abs(lam.to_value(u.angstrom) - wvl) < 0.03
        if not m.any():
            print(f"{label:<10} {wvl:9.3f}  遷移が見つからない")
            continue
        idx = np.where(m)[0]
        j = idx[np.argmax([g[:, i].max().value for i in idx])]
        gf = g[:, j].to_value(u.erg * u.cm ** 3 / u.s) * NHNE

        k = int(np.argmin(np.abs(wtar_idl - wvl)))
        gi = G_idl[k]
        r = gf.max() / gi.max()
        rows.append((label, wvl, gf))
        ratios.append(r)
        print(f"{label:<10} {wvl:9.3f} {gf.max():12.4e} {gi.max():12.4e} {r:7.3f} "
              f"{logT_idl[int(np.argmax(gf))]:13.2f} {logT_idl[int(np.argmax(gi))]:10.2f}")

    if ratios:
        a = np.array(ratios)
        print(f"\nG_fiasco / G_IDL :  median={np.median(a):.3f}  "
              f"min={a.min():.3f}  max={a.max():.3f}  "
              f"15%以内 {int(((a>=0.85)&(a<=1.15)).sum())}/{len(a)}")

    with open(OUT, "w") as f:
        f.write(f"# G(T) from fiasco, CHIANTI at /opt/ssw/packages/chianti/dbase\n")
        f.write(f"# abund={ABUND} ioneq={IONEQ} logNe={LOGNE} N_H/N_e={NHNE}\n")
        f.write("# nT nline\n")
        f.write(f"{len(logT_idl)} {len(rows)}\n# logT\n")
        f.write(" ".join(f"{x:.4f}" for x in logT_idl) + "\n# ion wvl\n")
        for label, wvl, _ in rows:
            f.write(f"{label} {wvl:.4f}\n")
        f.write("# G(T)\n")
        for _, _, gf in rows:
            f.write(" ".join(f"{x:.6e}" for x in gf) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
