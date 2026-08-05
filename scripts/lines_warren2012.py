"""Warren, Winebarger & Brooks (2012) ApJ 759, 141 が使った EIS 輝線リスト。

論文 Table 2 の 22 本の EIS 輝線と、eispac 同梱のフィットテンプレートの対応。
Table 2 の観測値は region 7 (2011-07-02 03:07:12, NOAA 1243) のもので、
講習会ではこの値と自分のフィット結果を照合する。

強度の単位: erg cm^-2 s^-1 sr^-1
"""

# (ion_label, rest_wavelength[A], eispac template, I_obs(Table2), sigma_I(Table2))
#
# テンプレートが複数ガウシアン成分を持つ場合、どれが目的の線かは
# template['line_ids'] を波長で照合して自動判定する（pick_component() 参照）。
# comp を手で書くと必ず間違える: 例えば fe_13_203_826.2c の第0成分は
# Fe XII 203.720 であって目的の Fe XIII 203.826 ではない。
LINES = [
    ("Si VII", 275.368, "si_07_275_368.1c.template.h5",   66.85,   14.76),
    ("Fe IX",  188.497, "fe_09_188_497.1c.template.h5",   71.23,   15.74),
    ("Fe IX",  197.862, "fe_09_197_862.1c.template.h5",   39.80,    8.79),
    ("Fe X",   184.536, "fe_10_184_536.1c.template.h5",  258.70,   57.04),
    ("Fe XI",  180.401, "fe_11_180_401.1c.template.h5",  795.27,  175.33),
    ("Fe XI",  188.216, "fe_11_188_216.2c.template.h5",  498.81,  109.77),
    ("S X",    264.233, "s__10_264_233.1c.template.h5",   55.27,   12.22),
    ("Si X",   258.375, "si_10_258_375.1c.template.h5",  213.53,   47.05),
    ("Fe XII", 192.394, "fe_12_192_394.1c.template.h5",  357.50,   78.67),
    ("Fe XII", 195.119, "fe_12_195_119.2c.template.h5", 1147.35,  252.44),
    ("Fe XIII",202.044, "fe_13_202_044.1c.template.h5", 1076.80,  236.95),
    ("Fe XIII",203.826, "fe_13_203_826.2c.template.h5", 1839.12,  404.73),
    ("Fe XIV", 264.787, "fe_14_264_787.1c.template.h5",  653.64,  143.84),
    ("Fe XIV", 270.519, "fe_14_270_519.2c.template.h5",  336.02,   73.95),
    ("Fe XV",  284.160, "fe_15_284_160.1c.template.h5", 5931.55, 1305.03),
    ("S XIII", 256.686, "s__13_256_686.1c.template.h5",  462.30,  101.78),
    ("Fe XVI", 262.984, "fe_16_262_984.1c.template.h5",  630.81,  138.82),
    ("Ar XIV", 194.396, "ar_14_194_396.2c.template.h5",   62.34,   13.74),
    ("Ca XIV", 193.874, "ca_14_193_874.2c.template.h5",  182.64,   40.21),
    ("Ca XV",  200.972, "ca_15_200_972.2c.template.h5",  127.92,   28.21),
    ("Ca XVI", 208.604, "ca_16_208_604.2c.template.h5",   31.12,    7.86),
    # ↓ このテンプレートは 192.700-193.200 を単一ガウシアンで塗るだけで、
    #   Fe XI 192.813 と O V のブレンドを解かない。論文の値の 4-5 倍が出る。
    #   scripts/ca17_blend.py の自作テンプレートを使うこと。
    ("Ca XVII",192.858, "ca_17_192_858.1c.template.h5",   41.75,    9.35),
]

# 論文 Table 2 の AIA 94A（Fe XVIII 分離後）: I_obs = 7.20 DN/s, sigma = 1.40


def pick_component(template, target_wvl):
    """template['line_ids'] を見て、目的波長に一番近いガウシアン成分の番号を返す。"""
    ids = [str(s) for s in template.template["line_ids"]]
    best, bestd = 0, 1e9
    for i, s in enumerate(ids):
        try:
            w = float(s.split()[-1])
        except ValueError:
            continue
        if abs(w - target_wvl) < bestd:
            best, bestd = i, abs(w - target_wvl)
    return best, ids


# 論文が注意している難所（検証で実際に確認したものは [確認済] を付けた）
NOTES = """
- [確認済] Ca XVII 192.858 は Fe XI 192.813 と O V 複合線にブレンドしている。
  論文は Ko et al. (2009) の方法で分離するが、**eispac にはこのブレンドを解く
  テンプレートが同梱されていない**。同じ波長域のテンプレートは 3 つあるが:
      ca_17_192_858.1c : [Ca XVII 192.858]                 ← 単一ガウシアン。ブレンドを吸う
      fe_11_192_813.2c : [Fe XI 192.813, O V 192.906]      ← Ca XVII が無い
      o__05_192_906.2c : [Fe XI 192.813, O V 192.906]      ← 同上
  どれも 3 成分同時フィットにならない。→ 自作テンプレートが必要。講習会の発展課題。

- [確認済] 複数成分テンプレートの成分順は波長順とは限らない。目的の線が
  第 0 成分でないものが 4 つある:
      fe_13_203_826.2c : [Fe XII 203.720, Fe XIII 203.826]  → comp 1
      fe_14_270_519.2c : [Mg VI 270.394, Fe XIV 270.519]    → comp 1
      ar_14_194_396.2c : [Mn X 194.327, Ar XIV 194.396]     → comp 1
      ca_14_193_874.2c : [Fe X 193.715, Ca XIV 193.874]     → comp 1
  必ず line_ids で照合すること（pick_component）。

- 論文は Ca 線の線幅を相互に拘束している:
    * Ca XVI 208.604 の幅 = Ca XV 200.972 の幅
    * Ca XVII 192.858 の幅 = Ca XIV 193.874 の幅 ±0.05 mA 以内
  eispac のデフォルトテンプレートはこの拘束を持たない（parinfo の 'tied' が空）。
  拘束を入れるには parinfo を書き換える。これも講習会の発展課題。

- Fe XIII 202.044 / 203.826 は論文でも I_obs/I_dem ~ 1.8 と大きくずれる。
  密度診断ペアであり、原子データ・密度効果の議論に使える。
"""
