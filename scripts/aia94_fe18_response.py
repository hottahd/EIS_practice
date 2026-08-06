"""AIA 94 Å の **Fe XVIII だけの**温度応答関数を作る（論文 p.6）。

論文:
    "we have computed a new response for this channel that only contains
     contributions from Fe XVIII... The response distributed with the official
     AIA software contains contributions from several of the known emission
     lines formed at lower temperatures."

つまり公式応答（`aiapy` や SSW の `aia_get_response`）はそのままでは使えない。
Fe XVIII の輝線だけを AIA 94 の**波長応答**に畳み込んで作り直す。

    R(T) = Σ_lines G_line(T) / (hc/λ) × A_eff(λ) × Ω_pix / (4π)

  G_line   [erg cm^3 s^-1]  ... CHIANTI から（scripts/idl/15_aia94_fe18_resp.pro）
  hc/λ     [erg/photon]     ... 光子数に直す
  A_eff(λ) [cm^2 DN/photon] ... aiapy の wavelength_response（DN 換算込み）
  Ω_pix    [sr]             ... 0.6″ 画素の立体角
  結果      [DN cm^5 s^-1 pix^-1]  →  DN/s = ∫ R(T) n_e n_H dh

    python scripts/aia94_fe18_response.py
"""
import sys

import numpy as np
import astropy.units as u
from astropy.constants import h, c

OUT = "work/aia94_fe18_response.txt"
SRC = "work/fe18_gofnt.txt"


def read_gofnt(path):
    with open(path) as f:
        lines = f.readlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("# nT nline"))
    nT, nline = (int(v) for v in lines[i + 1].split())
    vals = []
    k = i + 2
    while not lines[k].startswith("# logT"):
        k += 1
    k += 1
    def take(n):
        nonlocal k
        out = []
        while len(out) < n:
            out += [float(x) for x in lines[k].split()]
            k += 1
        return np.array(out[:n])
    logT = take(nT)
    while not lines[k].startswith("# wave"):
        k += 1
    k += 1
    wvl = take(nline)
    while not lines[k].startswith("# G(T)"):
        k += 1
    k += 1
    G = np.array([take(nT) for _ in range(nline)])
    return logT, wvl, G


def main():
    logT, wvl, G = read_gofnt(SRC)
    print(f"Fe XVIII: {len(wvl)} lines, nT={len(logT)}")

    try:
        from aiapy.response import Channel
    except ImportError:
        sys.exit("aiapy が要る: pip install aiapy")

    ch = Channel(94 * u.angstrom)
    wr = ch.wavelength_response()            # cm^2 DN / photon
    lam = ch.wavelength.to_value(u.angstrom)
    aeff = wr.to_value(u.cm ** 2 * u.DN / u.photon)
    print(f"AIA 94 波長応答: {lam.min():.1f}-{lam.max():.1f} A, "
          f"peak {aeff.max():.4e} cm^2 DN/photon at {lam[np.argmax(aeff)]:.2f} A")

    pix = (0.6 * u.arcsec).to_value(u.rad)   # AIA の画素サイズ
    omega = pix ** 2                          # sr / pixel

    R = np.zeros_like(logT)
    print(f"\n{'lambda':>9} {'A_eff':>12} {'Gmax':>12} {'寄与率':>8}")
    contrib = []
    for i, w in enumerate(wvl):
        a = float(np.interp(w, lam, aeff))
        ephot = (h * c / (w * u.angstrom)).to_value(u.erg)   # erg/photon
        term = G[i] / ephot * a * omega / (4 * np.pi)
        R += term
        contrib.append((term.max(), w, a, G[i].max()))
    tot = R.max()
    for cmax, w, a, gmax in sorted(contrib, reverse=True)[:8]:
        print(f"{w:9.3f} {a:12.4e} {gmax:12.4e} {cmax/tot:8.3f}")

    ip = int(np.argmax(R))
    print(f"\nR(T) peak = {R[ip]:.4e} DN cm^5 s^-1 pix^-1  at logT = {logT[ip]:.2f}")
    np.savetxt(OUT, np.column_stack([logT, R]),
               header="logT  R_FeXVIII [DN cm^5 s^-1 pix^-1]  (AIA 94, Fe XVIII lines only)\n"
                      "abund=Feldman1992 coronal, ioneq=chianti.ioneq, logNe=9.0")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
