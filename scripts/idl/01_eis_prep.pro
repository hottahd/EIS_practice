; Warren+2012 (region 7) の level-0 を eis_prep で level-1 にする。
;
;   bash scripts/idl/run_sswidl.sh scripts/idl/01_eis_prep.pro logs/01_eis_prep.log
;
; 注意: ssw_batch は `.run <file>` で実行する = **メインプログラム扱い**。
;       したがってファイル末尾に `end` が必須（無いと
;       "End of file encountered before end of program" で落ちる）。
;       逆に .run なので複数行の begin/endfor は普通に書ける。
;
; 較正の方針: Warren+2012 は 2012 年の論文なので、Del Zanna (2013) /
;   Warren et al. (2014) の打ち上げ後較正より前。つまり **打ち上げ前較正**。
;   eis_prep の既定 (/default, correct_sensitivity 無し) がそれに当たる。
;   これは eispac の `radcal=..._pre` と同じ土俵。→ 直接比較できる。

set_plot,'z'

l0dir  = getenv('EIS_L0DIR')
outdir = getenv('EIS_L1DIR')
f      = l0dir + '/eis_l0_20110702_030712.fits.gz'

print,'>>> level0 : '+f
print,'>>> exists : '+string(file_test(f))
print,'>>> outdir : '+outdir

t0 = systime(1)
eis_prep, f, /default, /save, /quiet, outdir=outdir
print,'>>> eis_prep took '+string((systime(1)-t0)/60.,format='(f6.2)')+' min'

print,'>>> produced:'
ff = file_search(outdir+'/eis_*_20110702_030712.fits*', count=n)
for i=0,n-1 do print,'>>>   '+ff[i]+'  '+string((file_info(ff[i])).size/1024L/1024L)+' MB'

exit

end
