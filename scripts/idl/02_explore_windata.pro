; level-1 の windata 構造を調べる。
; 目的は「eispac の (y,x) インデックスが SSW のどの軸・どの向きに対応するか」を
; 確定させること。ここを間違えると論文の箱と違う場所を測ることになる。
;
;   bash scripts/idl/run_sswidl.sh scripts/idl/02_explore_windata.pro logs/02.log

set_plot,'z'

l1 = getenv('EIS_L1')
print,'>>> l1 = '+l1

wd = eis_getwindata(l1, 195.119, /quiet)

print,'>>> windata tags:'
t = tag_names(wd)
for i=0,n_elements(t)-1 do print,'>>>   '+t[i]

s = size(wd.int)
print,'>>> int dims = '+strjoin(strtrim(s[1:s[0]],2),' x ')
print,'>>> nl='+strtrim(wd.nl,2)+' nx='+strtrim(wd.nx,2)+' ny='+strtrim(wd.ny,2)

print,'>>> wvl[0]='+strtrim(wd.wvl[0],2)+'  wvl[-1]='+strtrim(wd.wvl[wd.nl-1],2)
print,'>>> hdr line_id = '+wd.line_id

; 空間座標
print,'>>> solar_x: n='+strtrim(n_elements(wd.solar_x),2)+ $
      '  [0]='+string(wd.solar_x[0],format='(f8.2)')+ $
      '  [-1]='+string(wd.solar_x[n_elements(wd.solar_x)-1],format='(f8.2)')
print,'>>> solar_y: n='+strtrim(n_elements(wd.solar_y),2)+ $
      '  [0]='+string(wd.solar_y[0],format='(f8.2)')+ $
      '  [-1]='+string(wd.solar_y[n_elements(wd.solar_y)-1],format='(f8.2)')

; 波長方向に積分した強度マップの形（eispac のクイックルックと見比べる用）
img = total(wd.int,1)
print,'>>> img dims = '+strjoin(strtrim((size(img))[1:2],2),' x ')
print,'>>> img min/max = '+string(min(img),format='(e11.3)')+' / '+string(max(img),format='(e11.3)')

; テキストで書き出して Python 側と突き合わせる
openw,lun,getenv('OUTDIR')+'/idl_fe12_map.txt',/get_lun
printf,lun,(size(img))[1],(size(img))[2]
printf,lun,wd.solar_x
printf,lun,wd.solar_y
printf,lun,img
free_lun,lun
print,'>>> wrote idl_fe12_map.txt'

print,'>>> DONE'
exit
end
