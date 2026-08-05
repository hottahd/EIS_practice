"""EIS ラスターのクイックルック画像を作る（フィット無し・数秒）。

各スペクトルウィンドウの波長方向の積分値をとるだけ。輝線強度としては
連続光や隣接線が混ざるので不正確だが、「どこに何があるか」を掴んで
解析する箱（inter-moss 領域）を選ぶには十分速くて便利。

使い方:
    python quicklook_raster.py <eis_..._data.h5> [出力png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import eispac

# 見せたいウィンドウ: (波長, ラベル, 形成温度の目安)
PANELS = [
    (275.368, "Si VII 275.4",  "0.6 MK  (moss/足元)"),
    (184.536, "Fe X 184.5",    "1.1 MK"),
    (195.119, "Fe XII 195.1",  "1.6 MK"),
    (202.044, "Fe XIII 202.0", "1.8 MK"),
    (262.984, "Fe XVI 263.0",  "2.8 MK"),
    (193.874, "Ca XIV 193.9",  "3.5 MK"),
    (200.972, "Ca XV 201.0",   "4.5 MK"),
    (192.858, "Ca XVII 192.9", "5.6 MK  (ブレンド有)"),
]


def main(datafile, outpng="quicklook.png"):
    fig, axes = plt.subplots(1, len(PANELS), figsize=(3.0 * len(PANELS), 8.5))

    for ax, (wvl, label, temp) in zip(axes, PANELS):
        cube = eispac.read_cube(datafile, wvl)
        # 波長方向に積分（NaN は無視）
        img = np.nansum(cube.data, axis=2)
        vmin, vmax = np.nanpercentile(img, [1, 99.5])
        ax.imshow(np.sqrt(np.clip(img, 0, None)), origin="lower", aspect="auto",
                  cmap="inferno", vmin=np.sqrt(max(vmin, 0)), vmax=np.sqrt(vmax))
        ax.set_title(f"{label}\n{temp}", fontsize=9)
        ax.set_xlabel("x [pix]")
        if ax is axes[0]:
            ax.set_ylabel("y [pix]")
        else:
            ax.set_yticklabels([])

    idx = cube.meta["index"]
    fig.suptitle(f"{datafile.split('/')[-1]}   {idx['date_obs']} – {idx['date_end']}   "
                 f"xcen={idx['xcen']:.1f}\" ycen={idx['ycen']:.1f}\"", fontsize=11)
    fig.tight_layout()
    fig.savefig(outpng, dpi=110)
    print("wrote", outpng, "  raster shape (y,x):", img.shape)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "quicklook.png")
