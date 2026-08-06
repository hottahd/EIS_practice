"""eispac 用の Ca XVII 192.858 ブレンド分離テンプレート（Ko et al. 2009 相当）。

eispac 同梱の `ca_17_192_858.1c` は 192.700-193.200 を**単一ガウシアン**で塗るだけで、
Fe XI 192.813 と O V 複合線を丸ごと吸い込む。実測で論文値の 4.75-5.01 倍になる。

SSW 側（`scripts/idl/08_fit_ca.pro`）で論文どおりの拘束を実装したところ
**0.75-0.77** まで下がった。その「正解値」を再現できる Python 版を作る。

★ eispac の parinfo は `tied` に対応している
   （`eispac/core/generate_astropy_model.py` が MPFIT 形式の文字列を
    astropy のモデル関数に変換する。`p[2]` → `model.<param>`、`^` → `**`。
    `import` などはブラックリストで弾かれる）。
   **ただし HDF5 に保存すると `tied` が `<U4`（4 文字）に切り詰められる**ので、
   `0.32722*p[0]` のような式は**メモリ上でテンプレートを組み立てる**必要がある。

Ko et al. (2009) 相当の拘束（原子データは CHIANTI 9.0.1、
`scripts/idl/05b_ca17_ratios.pro` の出力）:

  - O V 多重線 6 本は「振幅スケール 1 つ」「共通シフト 1 つ」「共通線幅」だけ自由。
    成分間の相対強度はエネルギー比で固定:
      192.750:0.1749  192.797:0.3272  192.801:0.1312
      192.904:1.0000  192.911:0.1089  192.915:0.0088
    （光子数比 x 192.904/λ でエネルギー比に直したもの）
  - Fe XI 192.813 の強度は **同じ箱で測った Fe XI 188.216 から固定**:
      I(192.813)/I(188.216) = 0.21406 * (188.216/192.813) = 0.20896
  - Ca XVII 192.858 の線幅は **Ca XIV 193.874 の線幅に固定**（論文 p.6）

したがって 2 パス必要:
  1. Fe XI 188.216 と Ca XIV 193.874 を普通にフィット
  2. その結果を使って Ca XVII の窓を解く

    python scripts/ca17_template.py <data.h5> y0 y1 x0 x1
"""
import sys

import numpy as np

sys.path.insert(0, "scripts")
from fit_box_spectra import average_spectrum  # noqa: E402
from lines_warren2012 import pick_component  # noqa: E402
import eispac  # noqa: E402

# CHIANTI 9.0.1 から取得した O V 多重線の**光子数**比（192.904 = 1）
_OV_WVL_ALL = np.array([192.7500, 192.7970, 192.8010, 192.9040, 192.9110, 192.9150])
_OV_PH_ALL = np.array([0.17484, 0.32722, 0.13115, 1.00000, 0.10892, 0.00877])
_OV_EN_ALL = _OV_PH_ALL * (192.9040 / _OV_WVL_ALL)      # エネルギー比に直す

# ★ 分離できない成分は統合する。
#   EIS のサンプリングは 0.0223 A/画素、線幅は sigma ~ 0.030 A。
#   192.797 と 192.801 は 0.004 A 差、192.904/192.911/192.915 は 0.011 A 差で、
#   どちらも**原理的に分離不能**。強度重み付き波長で 1 本にまとめる。
#   （成分を増やしても情報は増えず、eispac の min_points >= n_params で弾かれる）
_GROUPS = [[0], [1, 2], [3, 4, 5]]
OV_WVL = np.array([np.average(_OV_WVL_ALL[g], weights=_OV_EN_ALL[g]) for g in _GROUPS])
OV_REL = np.array([_OV_EN_ALL[g].sum() for g in _GROUPS])

FE11_RATIO = 0.21406 * (188.216 / 192.813)    # I(192.813)/I(188.216)、エネルギー比
CA17_WVL = 192.853                            # CHIANTI 9.0.1 の値（論文表記は 192.858）
FE11_BLEND_WVL = 192.813

WMIN, WMAX = 192.60, 193.05


def fit_line(datafile, wvl, tname, y0, y1, x0, x1):
    """普通の 1 本フィット。(強度, 線幅) を返す。"""
    wave, inten, sig, _ = average_spectrum(datafile, wvl, y0, y1, x0, x1)
    t = eispac.read_template(eispac.data.get_fit_template_filepath(tname))
    comp, _ = pick_component(t, wvl)
    fit = eispac.fit_spectra(inten, t, wave=wave, errs=sig, ncpu=1, ignore_warnings=True)
    i = float(np.atleast_1d(fit.fit["int"][..., comp]).ravel()[0])
    w = float(np.atleast_1d(fit.fit["params"][..., 3 * comp + 2]).ravel()[0])
    return i, w


def build_template(fe11_amp, ca17_width, peak, background):
    """Ca XVII 窓の 8 成分テンプレートを組み立てる。

    成分の並び:
      0        : O V 192.905 群（振幅・中心・幅がこの複合線の基準）
      1-2      : 残りの O V 群（振幅は比で固定、中心は共通シフト、幅は共有）
      3        : Fe XI 192.813（振幅を fe11_amp に固定=tied、中心は共通シフト、幅は共有）
      4        : Ca XVII（振幅・中心は自由、幅は ca17_width に固定）

    自由パラメータは O V スケール / 共通シフト /
    Ca XVII 振幅 / Ca XVII 中心 / 背景 の **5 つだけ**（IDL 版と同じ）。
    """
    ref = 2          # 統合後の OV_WVL の中で 192.905 群の位置
    order = [ref] + [i for i in range(len(OV_WVL)) if i != ref]

    parinfo = []
    def add(value, lo=None, hi=None, fixed=0, tied=""):
        limited = [1 if lo is not None else 0, 1 if hi is not None else 0]
        parinfo.append(dict(value=float(value), fixed=int(fixed),
                            limited=np.array(limited, dtype=np.int16),
                            limits=np.array([lo or 0.0, hi or 0.0], dtype=float),
                            tied=tied))

    # --- 成分 0: O V 192.904
    add(peak * 0.3, lo=0.0)                                   # p[0] 振幅
    add(OV_WVL[ref], lo=OV_WVL[ref] - 0.05, hi=OV_WVL[ref] + 0.05)   # p[1] 中心
    # p[2] 幅。★ IDL 版 (08_fit_ca.pro) と同じく **Ca XIV の幅に固定**する。
    #   ここを自由にすると、Fe XI 成分の幅も一緒に太り、
    #   「強度を固定したはずの Fe XI」の積分強度が変わってしまう
    #   （振幅は固定でも I = A*sigma*sqrt(2pi) なので）。
    add(ca17_width, tied=f"{ca17_width:.6f}")
    line_ids = [f"O V {OV_WVL[ref]:.3f} (群)"]

    # --- 成分 1-5: 残りの O V
    for i in order[1:]:
        rel = OV_REL[i] / OV_REL[ref]
        dlam = OV_WVL[i] - OV_WVL[ref]
        add(peak * 0.3 * rel, tied=f"{rel:.5f}*p[0]")
        add(OV_WVL[i], tied=f"p[1]+{dlam:+.4f}")
        add(0.030, tied="p[2]")
        line_ids.append(f"O V {OV_WVL[i]:.3f} (群)")

    # --- 成分 3: Fe XI 192.813（振幅を外から固定）
    # ★ `fixed=1` では駄目。eispac は fit_spectra の中で scale_guess() を呼び、
    #   **fixed かどうかを見ずに**全ガウシアンの振幅を
    #   「その中心波長でのデータ値 - 背景」に置き換えてしまう
    #   （scale_guess.py の "compute new peaks" のループ）。
    #   一方 `tied` は mpfit が反復のたびに exec('p[i] = <式>') で強制するので
    #   scale_guess に上書きされない。**定数リテラルを tied に書く**のが正解。
    add(fe11_amp, tied=f"{fe11_amp:.6f}")
    add(FE11_BLEND_WVL, tied=f"p[1]+{FE11_BLEND_WVL - OV_WVL[ref]:+.4f}")
    add(0.030, tied="p[2]")
    line_ids.append("Fe XI 192.813")

    # --- 成分 7: Ca XVII（幅だけ固定）
    add(peak * 0.3, lo=0.0)
    add(CA17_WVL, lo=CA17_WVL - 0.05, hi=CA17_WVL + 0.05)
    add(ca17_width, tied=f"{ca17_width:.6f}")   # 幅も同じ理由で tied にする
    line_ids.append("Ca XVII 192.858")

    # --- 定数背景
    # ★ 初期値を 0 にしてはいけない。eispac は fit_spectra の中で
    #   scale_guess() を呼び、`scale = bkg_data / bkg_guess` で初期値を
    #   データに合わせて伸縮する。背景の初期値が 0 だとゼロ除算で inf になり、
    #   MPFIT が status=-16（"parameter or function value has become infinite"）
    #   で即死して、全パラメータが 0 のまま返ってくる。
    add(max(background, 1.0))

    tmpl = dict(n_gauss=len(line_ids), n_poly=1, wmin=WMIN, wmax=WMAX,
                line_ids=np.array(line_ids))

    # EISFitTemplate は filename/template/parinfo を受け取るコンストラクタを持つ。
    # `tied` は公式にドキュメント化された機能（user's guide の MPFIT の節）。
    return eispac.EISFitTemplate(filename="ca_17_192_858.8c (in memory)",
                                 template=tmpl, parinfo=parinfo)


def main():
    datafile = sys.argv[1]
    y0, y1, x0, x1 = (int(v) for v in sys.argv[2:6])
    print(f"box y=[{y0}:{y1}] x=[{x0}:{x1}]")

    # --- パス 1: Fe XI 188.216 と Ca XIV 193.874
    i_fe11, w_fe11 = fit_line(datafile, 188.216, "fe_11_188_216.2c.template.h5",
                              y0, y1, x0, x1)
    i_ca14, w_ca14 = fit_line(datafile, 193.874, "ca_14_193_874.2c.template.h5",
                              y0, y1, x0, x1)
    print(f"  Fe XI 188.216 : I={i_fe11:8.2f}  width={w_fe11:.5f}")
    print(f"  Ca XIV 193.874: I={i_ca14:8.2f}  width={w_ca14:.5f}")

    i_fe11_blend = i_fe11 * FE11_RATIO
    print(f"  → I(Fe XI 192.813) = {i_fe11_blend:8.2f}  (x {FE11_RATIO:.5f})")
    print(f"  → Ca XVII の線幅を {w_ca14:.5f} に固定")

    # --- パス 2: Ca XVII の窓
    wave, inten, sig, _ = average_spectrum(datafile, 192.858, y0, y1, x0, x1)
    m = (wave >= WMIN) & (wave <= WMAX)
    peak = float(np.nanmax(inten[m]))
    # 積分強度 → 振幅（I = A * sigma * sqrt(2 pi)）
    fe11_amp = i_fe11_blend / (w_ca14 * np.sqrt(2 * np.pi))

    base = float(np.nanmin(inten[m]))
    t = build_template(fe11_amp, w_ca14, peak, base)
    fit = eispac.fit_spectra(inten, t, wave=wave, errs=sig, ncpu=1, ignore_warnings=True)
    if fit is None:
        sys.exit("フィットが失敗した")

    ints = np.atleast_1d(fit.fit["int"]).ravel()
    print("\n  成分ごとの強度:")
    for name, v in zip(t.template["line_ids"], ints):
        print(f"    {name:<18} {v:9.2f}")
    i_ca17 = float(ints[len(t.template["line_ids"]) - 1])
    print(f"\n  Ca XVII 192.858 = {i_ca17:.2f}")
    print(f"  論文 Table 2    = 41.75   → ratio = {i_ca17/41.75:.3f}")
    print(f"  （参考）SSW 版 08_fit_ca.pro は同じ箱で 0.75-0.77）")
    print(f"  （参考）eispac 同梱の 1 成分テンプレートは 4.75-5.01）")


if __name__ == "__main__":
    main()
