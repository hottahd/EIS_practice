#!/bin/bash
# 参照論文の arXiv 版を papers/ に落とす。
#
# PDF はリポジトリに入れていない（再配布の許諾が無い。papers/README.md 参照）。
# arXiv からの個人利用のダウンロードは自由なので、各自の手元にはこれで揃う。
#
#   bash papers/fetch_papers.sh
#
# 出版社版 (ApJ) がほしい場合は papers/README.md の DOI リンクから。
# 所属機関の購読か、NASA ADS 経由で入手できる。

set -u
cd "$(dirname "$0")"
mkdir -p refs

get() {   # get <arxiv_id> <出力ファイル名>
    if [ -s "$2" ]; then
        echo "  既にある: $2"
        return
    fi
    echo "  取得中: arXiv:$1 -> $2"
    # arXiv は連続アクセスを嫌うので間隔を空ける（公式の案内どおり）
    curl -sSL --max-time 120 -o "$2" "https://arxiv.org/pdf/$1" || {
        echo "  ! 失敗: $1  （https://arxiv.org/abs/$1 から手動で取得してください）"
        rm -f "$2"
        return
    }
    sleep 3
}

echo "主教材 (Warren, Winebarger & Brooks 2012, ApJ 759, 141)"
get 1204.3220 "Warren_2012_arxiv_1204.3220.pdf"

echo "参照論文"
get 1009.5976 "refs/arxiv_1009.5976.pdf"    # Warren+2011, ApJ 734, 90
get 1106.5057 "refs/arxiv_1106.5057.pdf"    # Winebarger+2011, ApJ 740, 2
get 1107.4480 "refs/arxiv_1107.4480.pdf"    # Tripathi+2011, ApJ 740, 111

echo
echo "完了。取得したもの:"
ls -lh ./*.pdf refs/*.pdf 2>/dev/null || echo "  (何も取得できていません)"
echo
echo "※ scripts/extract_paper_boxes.py は Figure から箱の座標を実測するので、"
echo "   図の解像度が出版社版と arXiv 版で違うと結果がずれることがあります。"
