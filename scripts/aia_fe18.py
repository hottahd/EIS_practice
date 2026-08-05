"""AIA 94 A から Fe XVIII 93.92 A の寄与を分離する（Warren et al. 2012, Appendix）。

94 A チャンネルには 7 MK の Fe XVIII が入っているが、それより遥かに低温の
輝線に汚染されている。論文は 171 A と 193 A の混合から「warm 成分」を
経験的に見積もって引く:

    x        = ( f*I_171 + (1-f)*I_193 ) / 116.54          f = 0.31
    I_94warm = 0.39 * ( a1 + a2*x + a3*x^2 + a4*x^3 )
    a        = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]
    I_FeXVIII = I_94 - I_94warm

- 強度は全て DN/s（露光時間で規格化済み）。
- x は 30 で頭打ちにする（「高強度側にデータが無いから」と論文にある）。
  実際この 3 次式は x = 27.4 で極大になるので、30 はちょうど折り返し付近。
- 較正の元データは 2010-03-22 12-13 UT の 1 時間平均（明るい輝点＋静穏領域）。

【重要: 論文 Eq.(A1) の指数について】
  論文に印刷されている式は

      I_94warm = 0.39 * sum_{i=1}^{4} a_i * x^i

  だが、この字面どおりに計算すると warm 成分が観測値を桁違いに超える。
  2011-07-02 の実データ (NOAA 1243) で両者を比べた結果:

      x の範囲   観測 I_94   sum a_i x^i   a_i x^(i-1)
        2 - 3       1.26        2.74          1.12
        5 - 8       3.62       21.17          3.47
       12 - 18     10.28      135.03          9.73
       25 - 32     22.92      441.32         16.50

  さらに Fe XVIII 総強度 (>2 DN/s) は
      sum a_i x^i     -> 1.37e3   (論文 6.18e4 の 1/45。負値が 77% を占める)
      a_i x^(i-1)     -> 7.68e4   (論文 6.18e4 とオーダー一致)

  よって **正しくは定数項つきの 3 次式** (a_i x^(i-1), i=1..4) であり、
  印刷された指数は組版上の誤りと判断した。x=1（中央値）で
  I_94warm = 0.389 になる点はどちらの読み方でも同じなので、
  そこだけでは判別できない。実データで初めて決まる。

【適用限界（講習会で必ず伝えること）】
  * フレア中は使えない。Fe XXIV 192.04 A が 193 A チャンネルに入るため。
  * 非常に明るい 1 MK 放射（moss）でも破綻する。
  * この経験式は AIA の感度劣化補正を **していない** データに対して導かれている。
    aiapy.calibrate.degradation() を掛けると係数と整合しなくなる。
    論文どおりの再現をするなら劣化補正は掛けない。
    （ただし 2010 年較正の式を 2011 年以降のデータに使うことの是非は議論の余地あり。
      これ自体が講習会の良い議論の種になる。）
"""
import numpy as np

F_MIX = 0.31
NORM_COMPOSITE = 116.54
NORM_94 = 0.39
POLY_COEFFS = [-7.31e-2, 9.75e-1, 9.90e-2, -2.84e-3]   # a_1 .. a_4 （x^0 .. x^3 に掛かる）
X_MAX = 30.0


def warm_94(i171, i193, f=F_MIX, clip=True):
    """171 A と 193 A から 94 A チャンネルの warm 成分 [DN/s] を推定する。"""
    x = (f * np.asarray(i171, float) + (1.0 - f) * np.asarray(i193, float)) / NORM_COMPOSITE
    if clip:
        x = np.clip(x, 0.0, X_MAX)
    out = np.zeros_like(x)
    for i, a in enumerate(POLY_COEFFS):     # a_1 x^0 + a_2 x^1 + a_3 x^2 + a_4 x^3
        out += a * x**i
    return NORM_94 * out


def fe18(i94, i171, i193, f=F_MIX):
    """AIA 94 A から Fe XVIII 成分 [DN/s] を取り出す。"""
    return np.asarray(i94, float) - warm_94(i171, i193, f=f)


def _selftest():
    """論文の記述と整合しているかの自己確認。"""
    # 中央値 (x=1) では I_94warm が 94 A の中央値 0.39 DN/s になるはず
    print(f"x=1 (median) -> I_94warm = {warm_94(116.54, 116.54):.4f} DN/s"
          f"  (論文の中央値 0.39 と一致すべき)")

    # 多項式の極大 = 頭打ち 30 の根拠
    xs = np.linspace(0, 40, 4001)
    p = np.array([sum(a * xx**i for i, a in enumerate(POLY_COEFFS)) for xx in xs])
    print(f"多項式の極大は x = {xs[np.argmax(p)]:.1f}  (論文の頭打ち 30 の根拠)")
    print(f"x=30 での I_94warm = {NORM_94*p[np.searchsorted(xs, 30)]:.2f} DN/s"
          f"  (観測の I_94 最大 ~30 DN/s と同程度で妥当)")


if __name__ == "__main__":
    _selftest()
