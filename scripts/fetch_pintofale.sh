#!/bin/bash
# PINTofALE (Kashyap & Drake 1998, ApJ 503, 450) を取得する。
#
# 公式サイトは「ダウンロードは廃止、著者に連絡して MEGA から」と書いているが、
# doc/installation.txt に載っている旧 URL は 2026-08-06 時点で全部生きている。
# いつ消えてもおかしくないので、必要になったら早めに取っておくこと。
#
#   usage: bash scripts/fetch_pintofale.sh [/path/to/install/dir]
#
# 既定の展開先は $SCRATCH か、無ければカレント。
# 詳細は docs/00_log.md「PINTofALE は普通に入手できた」を参照。

set -u
BASE=https://hea-www.harvard.edu/PINTofALE
DEST=${1:-${SCRATCH:-$PWD}}
DIST="$DEST/PINTofALE_dist"

mkdir -p "$DIST" || exit 1
cd "$DIST" || exit 1

# PoA_current      = pro + ardb + doc + 連続光/SPEX/APED 放射率   (236 MB)
# PoA_doc_current  = ドキュメント一式（AIA_DEM.html はここ）      (3.3 MB)
# PoA_chianti_current = CHIANTI 輝線放射率（PoA 形式）            (184 MB)
# ※ PoA_atomdb_current.tar.gz は 404。存在しない。
for f in PoA_current PoA_doc_current PoA_chianti_current; do
    echo "--- $f.tar.gz"
    wget -c "$BASE/$f.tar.gz" || echo "!! $f の取得に失敗"
done

cd "$DEST" || exit 1
# ardb/pimms/ のパーミッションが壊れていて tar が途中で止まるので、
# chmod を挟んで 2 回展開する。
for pass in 1 2; do
    for f in "$DIST"/*.tar.gz; do
        tar xzf "$f" 2>/dev/null
    done
    chmod -R u+rwX PINTofALE 2>/dev/null
done

echo
echo "install: $DEST/PINTofALE"
du -sh "$DEST/PINTofALE" 2>/dev/null
echo
echo "IDL/GDL からは:"
echo "  setenv PINTofALE $DEST/PINTofALE      # 自動検出は当てにならない"
echo "  IDL> .run \$PINTofALE/pro/scrypt/initale"
echo
echo "MCMC_DEM だけ回すなら放射率 DB (emissivity/, 568 MB) は不要。"
echo "  → pro/fitting/mcmc_dem.pro とその依存 27 ファイルだけで足りる。"
echo "  → G(T) は Python (fiasco) 側で作って emis[nT,nline] として渡す。"
