"""箱を総当たりし、**輝線ごとの ratio を全部記録する**。

「その輝線は箱を変えれば論文に合うのか」を輝線ごとに答えるため。
集計スコアだけ見ていると、上位の箱が同じ場所に固まって
「どの箱でも合わない」という誤った結論になる（実際に一度そう誤った）。
"""
import sys
import numpy as np
sys.path.insert(0, "scripts")
from compare_table2 import fit_box
from lines_warren2012 import LINES

def main():
    datafile = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "work/perline_scan.csv"
    names = [f"{i} {w:.3f}" for i, w, *_ in LINES]
    with open(out, "w") as f:
        f.write("y0,y1,x0,x1," + ",".join(n.replace(" ", "_") for n in names) + "\n")
        for dy, dx in [(20, 6), (30, 8), (45, 12)]:
            for y0 in range(120, 430 - dy, 20):
                for x0 in range(4, 52 - dx, 8):
                    rows = fit_box(datafile, y0, y0 + dy, x0, x0 + dx)
                    d = {f"{r['ion']} {r['wvl']:.3f}": r["ratio"] for r in rows}
                    f.write(f"{y0},{y0+dy},{x0},{x0+dx}," +
                            ",".join(f"{d.get(n, float('nan')):.4f}" for n in names) + "\n")
                    f.flush()
                    print(f"y=[{y0}:{y0+dy}] x=[{x0}:{x0+dx}] done", flush=True)

if __name__ == "__main__":
    main()
