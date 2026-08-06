"""demregpy（正則化インバージョン）で DEM を解き、PINTofALE の MCMC_DEM と比べる。

講習会モジュール 7。論文は PINTofALE の MCMC_DEM（Kashyap & Drake 1998）を使うが、
それには IDL が要る。Python だけで完結する経路として demregpy
（Hannah & Kontar 2012 の正則化インバージョン）を使い、
**同じ観測強度・同じ寄与関数**で解いて結果を比べる。

  「手法によって DEM がどう変わるか」自体を教材にする。

入力（IDL 側の成果物をそのまま使う）:
  work/gofnt_chianti901.txt      ... G(T)。scripts/idl/09_gofnt.pro
  work/idl_intensities_tied.csv  ... 22 輝線の強度。scripts/idl/04_fit_box_tied.pro
  work/idl_ca_intensities.csv    ... Ca 線（論文の拘束つき）。08_fit_ca.pro
  work/aia94_fe18_response.txt   ... AIA 94 の Fe XVIII 応答（任意）
  work/mcmc_dem_result.txt       ... 比較対象の MCMC_DEM の結果

★ 単位のつじつま（3 度も踏んだので明示しておく）
  観測強度 I は erg cm^-2 s^-1 sr^-1、G(T) は erg cm^3 s^-1 で
      I = (1/4π) ∫ G(T) n_e n_H ds
  なので demregpy に渡す応答は **G/(4π)**、DEM は [cm^-5 K^-1] で返る。
  AIA の応答 R(T) は既に 1/(4π) と画素立体角込みなので **そのまま渡す**。

    python scripts/dem_demregpy.py [--aia 6.668]
"""
import argparse
import csv
import sys

import numpy as np

GOFNT = "work/gofnt_chianti901_005.txt"   # demregpy は細かい温度グリッドが要る
INTEN = "work/idl_intensities_tied.csv"
CAINT = "work/idl_ca_intensities.csv"
AIARESP = "work/aia94_fe18_response.txt"
MCMC = "work/mcmc_dem_result.txt"
OUT = "work/demregpy_result.txt"


def read_gofnt(path):
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
        names.append(" ".join(p[:-2]))
        wtar.append(float(p[-2]))
        k += 1
    skip_to("# G(T)")
    G = np.array([take(nT) for _ in range(nline)])
    return logT, names, np.array(wtar), G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aia", type=float, default=None,
                    help="AIA 94 Fe XVIII の観測値 [DN/s]。指定すると拘束に加える")
    ap.add_argument("--tlo", type=float, default=5.5)
    ap.add_argument("--thi", type=float, default=7.1)
    ap.add_argument("--reg", type=float, default=1.0, help="reg_tweak（目標 chi2）")
    ap.add_argument("--gloci", type=int, default=0,
                    help="0=自己正規化の重み, 1=EM loci を重みに使う")
    ap.add_argument("--norm-from-mcmc", action="store_true",
                    help="MCMC の解を dem_norm0（初期の重み）に使う")
    a = ap.parse_args()

    logT, names, wtar, G = read_gofnt(GOFNT)
    keep = (logT >= a.tlo) & (logT <= a.thi)
    logT, G = logT[keep], G[:, keep]
    print(f"G(T): {len(names)} 輝線, logT {logT[0]:.2f}-{logT[-1]:.2f} "
          f"({len(logT)} 点, {logT[1]-logT[0]:.2f} dex)")

    # --- 観測強度
    iobs = np.zeros(len(names))
    esig = np.zeros(len(names))
    ipap = np.zeros(len(names))
    for r in csv.DictReader(open(INTEN)):
        w = float(r["wvl"])
        k = int(np.argmin(np.abs(wtar - w)))
        if abs(wtar[k] - w) < 0.01:
            iobs[k] = float(r["I_idl"])
            ipap[k] = float(r["I_paper"])
            esig[k] = float(r["sig_paper"]) / float(r["I_paper"]) * float(r["I_idl"])
    # Ca 線は論文の拘束つきの値で上書き。
    # ★ 08_fit_ca.pro の CSV は「波長, 強度, 中心, 幅, chi2r, 備考」の 6 列。
    #   ヘッダに ion 列が余分に書かれていた時期があるので**位置で読む**。
    with open(CAINT) as f:
        for line in f:
            p = [x.strip() for x in line.split(",")]
            try:
                w, ical = float(p[0]), float(p[1])
            except (ValueError, IndexError):
                continue
            if not (abs(w - 208.604) < 0.01 or abs(w - 192.858) < 0.01):
                continue
            k = int(np.argmin(np.abs(wtar - w)))
            iobs[k] = ical
            esig[k] = 0.22 * ical
            print(f"  {names[k]} {w:.3f} を論文の拘束つきの値 {ical:.2f} で置換")

    ok = iobs > 0
    labels = [f"{names[i]} {wtar[i]:.3f}" for i in np.where(ok)[0]]
    dn = iobs[ok]
    edn = esig[ok]
    # ★ I = (1/4π)∫G n_e n_H ds なので応答は G/(4π)
    tresp = (G[ok] / (4 * np.pi)).T          # (nT, nf)

    # --- AIA Fe XVIII を加える（応答は既に 1/(4π) と画素立体角込み）
    if a.aia is not None:
        d = np.loadtxt(AIARESP)
        R = np.interp(logT, d[:, 0], d[:, 1])
        tresp = np.column_stack([tresp, R])
        dn = np.append(dn, a.aia)
        edn = np.append(edn, 0.19 * a.aia)
        labels.append("AIA 94 FeXVIII")
        print(f"  AIA Fe XVIII を追加: {a.aia:.3f} DN/s")

    print(f"拘束に使う本数: {len(dn)}")

    # --- 温度ビンの端
    dlt = logT[1] - logT[0]
    tedges = 10 ** np.append(logT - dlt / 2, logT[-1] + dlt / 2)

    dem_norm0 = None
    if a.norm_from_mcmc:
        m0 = np.loadtxt(MCMC)
        dem_norm0 = np.interp(logT, m0[:, 0], m0[:, 2])
        dem_norm0 = dem_norm0 / dem_norm0.max()
        print("  dem_norm0 に MCMC の解を使う（相対値のみ効く）")

    from demregpy import dn2dem
    dem, edem, elogt, chisq, dn_reg = dn2dem(
        dn, edn, tresp, logT, tedges, reg_tweak=a.reg, max_iter=30,
        gloci=a.gloci, dem_norm0=dem_norm0, warn=True)

    dem = np.atleast_1d(np.squeeze(dem))
    edem = np.atleast_1d(np.squeeze(edem))
    dn_reg = np.atleast_1d(np.squeeze(dn_reg))
    print(f"\nreduced chi2 = {float(np.squeeze(chisq)):.3f}")

    print(f"\n{'line':<20} {'I_obs':>11} {'I_dem':>11} {'R':>7}")
    for lab, o, p in zip(labels, dn, dn_reg):
        print(f"{lab:<20} {o:11.4e} {p:11.4e} {o/p if p>0 else np.nan:7.3f}")

    # --- MCMC と比較
    T = 10 ** logT
    em = dem * T * np.log(10) * dlt          # DEM[cm^-5/K] -> ビンあたりの EM[cm^-5]
    ip = int(np.nanargmax(em))
    print(f"\ndemregpy: EM ピーク logT={logT[ip]:.2f} ({T[ip]/1e6:.2f} MK)")

    def slope(t0, t1, y):
        m = (logT >= t0) & (logT <= t1) & (y > 0)
        return np.polyfit(logT[m], np.log10(y[m]), 1)[0] if m.sum() > 2 else np.nan

    print(f"  alpha (6.0-6.6) = {slope(6.0,6.6,em):+.2f}   論文 2.9")
    print(f"  beta  (6.6-7.0) = {-slope(6.6,7.0,em):+.2f}   論文 9.0")

    try:
        m = np.loadtxt(MCMC)
        # ★ PINTofALE の DEM は [cm^-5/logK] なので EM(ビン) = DEM × ΔlogT。
        #   demregpy は [cm^-5/K] なので EM(ビン) = DEM × ΔT = DEM × T ln10 ΔlogT。
        #   ここを揃えないと 1e7 ずれる（実際に踏んだ）。
        dlt_m = float(np.median(np.diff(m[:, 0])))
        lt_m, em_m = m[:, 0], m[:, 2] * dlt_m
        print(f"\n{'logT':>6} {'EM(demregpy)':>14} {'EM(MCMC)':>14} {'比':>7}")
        for i, t in enumerate(logT):
            j = int(np.argmin(np.abs(lt_m - t)))
            if abs(lt_m[j] - t) > 1e-6:
                continue
            r = em[i] / em_m[j] if em_m[j] > 0 else np.nan
            print(f"{t:6.2f} {em[i]:14.4e} {em_m[j]:14.4e} {r:7.3f}")
        ipm = int(np.nanargmax(em_m))
        print(f"\nMCMC_DEM: EM ピーク logT={lt_m[ipm]:.2f} ({10**lt_m[ipm]/1e6:.2f} MK)")
    except Exception as e:
        print(f"(MCMC の結果を読めなかった: {e})")

    np.savetxt(OUT, np.column_stack([logT, dem, edem, elogt.squeeze(), em]),
               header="logT  DEM[cm^-5/K]  eDEM  elogT  EM=DEM*T*ln10*dlogT [cm^-5]")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
