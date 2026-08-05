"""候補の箱で 22 輝線をフィットし、論文 Table 2 と比較して要約する。

出力の見方:
  ratio  = 自分の強度 / 論文 Table 2 の強度
  各輝線で ratio がどれくらい 1 に近いか、そして
  「冷たい線と熱い線で ratio の傾向が揃っているか」を見る。

  * ratio が全体に一定倍率でずれる → 較正・面積規格化の違い（許容範囲）
  * 冷たい線と熱い線で ratio が逆向きにずれる → **箱の場所が違う**
    （暖かいループ寄り / 高温コア寄りを選んでいる）

使い方:
    python compare_table2.py <eis_..._data.h5> "y0:y1,x0:x1" ["y0:y1,x0:x1" ...]
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from fit_box_spectra import average_spectrum
from lines_warren2012 import LINES, pick_component
import eispac

# おおまかな形成温度（log T）。ratio の温度依存を見るために使う。
LOGT = {
    "Si VII": 5.8, "Fe IX": 5.9, "Fe X": 6.05, "Fe XI": 6.15, "S X": 6.15,
    "Si X": 6.15, "Fe XII": 6.2, "Fe XIII": 6.25, "Fe XIV": 6.3,
    "Fe XV": 6.35, "S XIII": 6.4, "Fe XVI": 6.45, "Ar XIV": 6.5,
    "Ca XIV": 6.55, "Ca XV": 6.65, "Ca XVI": 6.7, "Ca XVII": 6.75,
}


def fit_box(datafile, y0, y1, x0, x1, verbose=False):
    rows = []
    for ion, wvl, tname, i_paper, sig_paper in LINES:
        wave, inten, sig, npix = average_spectrum(datafile, wvl, y0, y1, x0, x1)
        t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
        comp, ids = pick_component(t, wvl)
        fit = eispac.fit_spectra(inten, t, wave=wave, errs=sig, ncpu=1,
                                 ignore_warnings=True)
        if fit is None:
            continue
        i_fit = float(np.atleast_1d(fit.fit["int"][..., comp]).ravel()[0])
        rows.append(dict(ion=ion, wvl=wvl, i_fit=i_fit, i_paper=i_paper,
                         ratio=i_fit / i_paper, logt=LOGT[ion]))
        if verbose:
            print(f"  {ion} {wvl:.3f}: {i_fit:9.2f} / {i_paper:8.2f} = {i_fit/i_paper:.2f}")
    return rows


def summarize(rows, label):
    # Ca XVII は eispac のテンプレートではブレンドを解けないので要約から外す
    use = [r for r in rows if r["ion"] != "Ca XVII"]
    ratios = np.array([r["ratio"] for r in use])
    logt = np.array([r["logt"] for r in use])
    cool = ratios[logt <= 6.35]
    hot = ratios[logt >= 6.45]
    # 温度依存の傾き（log ratio vs logT）。0 に近いほど箱の場所が論文と整合。
    slope = np.polyfit(logt, np.log10(ratios), 1)[0]
    print(f"{label:<26s} median={np.median(ratios):5.2f}  "
          f"cool(logT<=6.35)={np.median(cool):5.2f}  hot(logT>=6.45)={np.median(hot):5.2f}  "
          f"傾き={slope:+6.2f}  ばらつき={np.std(np.log10(ratios)):.2f}")
    return abs(slope)


if __name__ == "__main__":
    datafile = sys.argv[1]
    boxes = sys.argv[2:]
    print("傾きが 0 に近いほど、論文と同じ温度組成の場所を見ていることになる。\n")
    results = []
    for b in boxes:
        yy, xx = b.split(",")
        y0, y1 = (int(v) for v in yy.split(":"))
        x0, x1 = (int(v) for v in xx.split(":"))
        rows = fit_box(datafile, y0, y1, x0, x1)
        s = summarize(rows, f"y=[{y0}:{y1}] x=[{x0}:{x1}]")
        results.append((s, b, rows))
    print()
    best = min(results)
    print(f"★ 温度依存が最も小さい箱: {best[1]}")
    print("\n--- その箱の全輝線 ---")
    print(f"{'line':<16s} {'I_fit':>10s} {'I_Warren':>10s} {'ratio':>7s}")
    for r in best[2]:
        name = f"{r['ion']} {r['wvl']:.3f}"
        print(f"{name:<16s} {r['i_fit']:10.2f} {r['i_paper']:10.2f} {r['ratio']:7.2f}")
