#!/bin/bash
# SSWIDL をヘッドレスで走らせるラッパ。
#   usage: bash scripts/idl/run_sswidl.sh prog.pro [logfile]
#
# prog.pro は IDL のバッチファイル。末尾に必ず `exit` を書くこと
# （書かないとプロンプトで止まる）。図を描くなら先頭で set_plot,'z'。
#
# このマシンの実情:
#   IDL  9.2.0 (NV5, Nagoya University ライセンス)  /usr/local/nv5/idl92
#   SSW  /opt/ssw -> /lustre/sc/ssw  (2021-04 版, hinode/eis + packages/poa + chianti 9.0.1)
#   SSWDB /opt/sswdb -> /lustre/sc/sswdb
# lustre が遅いので SSW の起動だけで 1-3 分かかる。長めに待つこと。
set -u
export IDL_DIR=${IDL_DIR:-/usr/local/nv5/idl92}
export SSW=${SSW:-/opt/ssw}
export SSWDB=${SSWDB:-/opt/sswdb}
export SSW_INSTR=${SSW_INSTR:-"eis"}      # 増やすほど起動が遅くなる

# PINTofALE を使うときは、**IDL 起動前に** IDL_PATH に入れておくこと。
# ssw_batch は `.run <file>` = メインプログラムとして一括コンパイルするので、
# スクリプト内で実行時に !path を書き換えても遅い。コンパイル時に関数が
# 見つからないと `varsmooth(...)` が「未定義変数の添字」と解釈され、
# キーワード指定のところで「% Syntax error」になる（原因が分かりにくい）。
if [ -n "${PINTofALE:-}" ]; then
  export IDL_PATH="+$PINTofALE/pro${IDL_PATH:+:$IDL_PATH}"
fi
PROG=$(readlink -f "$1")
LOG=${2:-$PROG.log}

RC=$(mktemp /tmp/sswidl_XXXXXX.csh)
cat > "$RC" <<CSH
setenv IDL_DIR "$IDL_DIR"
setenv SSW "$SSW"
setenv SSWDB "$SSWDB"
setenv SSW_INSTR "$SSW_INSTR"
setenv IDL_PATH "${IDL_PATH:-}"
unsetenv DISPLAY
source \$SSW/gen/setup/setup.ssw /quiet
\$SSW/gen/bin/ssw_batch "$PROG" "$LOG"
CSH
csh -f "$RC"
rc=$?
rm -f "$RC"
exit $rc
