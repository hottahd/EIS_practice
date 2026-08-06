; despike (宇宙線除去) を切って eis_prep をやり直す。
; Fe XII 195.119 の箱内平均スペクトルに 195.104 A で不自然なくぼみがあり、
; 積分強度が NRL level-1 (eispac) より 24% 小さい。
; eis_prep の despike が明るい線のコアを削っている疑いを確かめる。
set_plot,'z'
f = getenv('EIS_L0DIR')+'/eis_l0_20110702_030712.fits.gz'
print,'>>> nocr run'
eis_prep, f, /default, /save, /quiet, /nocr, outdir=getenv('EIS_L1DIR')
print,'>>> DONE'
exit
end
