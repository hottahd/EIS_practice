; 打ち上げ後の EIS 較正カーブを作る（★リストの項目 2）。
;
;   OUTDIR=work bash scripts/idl/run_sswidl.sh scripts/idl/11_calcurve.pro logs/11.log
;
; eispac は打ち上げ前較正しか持っていない（radcal/<win>_pre）。
; SSW には Warren 本人の eis_recalibrate_intensity があり、
;   - /gdz なし → Warren et al. (2014) の NRL 補正
;   - /gdz あり → Del Zanna (2013) の補正
; が使える。観測日 2011-07-02 について波長ごとの倍率を出し、
; Python 側 (eispac の read_cube(radcal=...)) に渡せる形で保存する。
;
; 動機: 箱をどこに動かしても
;   Si VII 275.368 = 0.36-0.39,  Fe XVI 262.984 = 0.62-0.70,
;   S XIII 256.686 = 0.71-0.77
; が論文より低いまま動かない。波長依存の較正が効いていないか確かめる。

set_plot,'z'
!quiet = 1

outdir = getenv('OUTDIR')
date   = '2011-07-02T03:38:08'

;; --- 論文 Table 2 の 22 波長
wv = [275.368d, 188.497d, 197.862d, 184.536d, 180.401d, 188.216d, $
      264.233d, 258.375d, 192.394d, 195.119d, 202.044d, 203.826d, $
      264.787d, 270.519d, 284.160d, 256.686d, 262.984d, 194.396d, $
      193.874d, 200.972d, 208.604d, 192.858d]
nm = ['Si VII','Fe IX','Fe IX','Fe X','Fe XI','Fe XI','S X','Si X', $
      'Fe XII','Fe XII','Fe XIII','Fe XIII','Fe XIV','Fe XIV','Fe XV', $
      'S XIII','Fe XVI','Ar XIV','Ca XIV','Ca XV','Ca XVI','Ca XVII']

print,'>>> date = '+date
print,'>>> line          wvl      NRL(W14)   GDZ(2013)'
for k=0,n_elements(wv)-1 do begin
  fn = eis_recalibrate_intensity(date, wv[k], 1.0d, /quiet)
  fg = eis_recalibrate_intensity(date, wv[k], 1.0d, /gdz, /quiet)
  print,'>>> '+string(nm[k],format='(a-9)')+string(wv[k],format='(f10.3)')+ $
        string(fn,format='(f11.4)')+string(fg,format='(f12.4)')
endfor

;; --- 連続的なカーブも吐く（SW 165-215, LW 245-291）
wsw = 170.d + 0.1d*dindgen(421)   ; SW 170-212 A（範囲外だと eis_ltds が落ちる）
wlw = 246.d + 0.1d*dindgen(451)   ; LW 246-291 A
openw,lun,outdir+'/eis_calcurve_20110702.txt',/get_lun
printf,lun,'# EIS post-launch calibration factors for '+date
printf,lun,'# I_corrected = I_preflight * factor'
printf,lun,'# columns: wavelength[A]  NRL_Warren2014  DelZanna2013'
;; eis_ea_gdz は SW/LW の端で落ちるので 1 点ずつ catch する
wall = [wsw, wlw]
for k=0,n_elements(wall)-1 do begin
  catch, err
  if err ne 0 then begin & catch,/cancel & continue & endif
  fn = eis_recalibrate_intensity(date, wall[k], 1.0d, /quiet)
  fg = eis_recalibrate_intensity(date, wall[k], 1.0d, /gdz, /quiet)
  catch,/cancel
  printf,lun, wall[k], fn, fg, format='(f9.2,2e15.6)'
endfor
free_lun,lun
print,'>>> wrote '+outdir+'/eis_calcurve_20110702.txt'

print,'>>> DONE'
exit
end
