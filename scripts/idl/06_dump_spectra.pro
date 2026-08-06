; 箱内平均スペクトルを 22 窓ぶんテキストに吐く。
; 目的: 「IDL と Python の差」が **データ（eis_prep vs NRL level-1）** 由来なのか
;       **フィッター** 由来なのかを切り分けること。
; 同じスペクトルを Python 側でも同じモデルでフィットして比べる。
;
;   EIS_L1=... BOX=244,274,32,40 OUTDIR=work \
;     bash scripts/idl/run_sswidl.sh scripts/idl/06_dump_spectra.pro logs/06.log

set_plot,'z'
!quiet = 1

l1     = getenv('EIS_L1')
outdir = getenv('OUTDIR')
box    = fix(strsplit(getenv('BOX'),',',/extract))
y0 = box[0] & y1 = box[1]-1 & x0 = box[2] & x1 = box[3]-1

;; テンプレート表（波長だけ使う）
nmax = 200
c_ion=strarr(nmax) & c_wvl=dblarr(nmax)
openr,tlun,outdir+'/eispac_templates.csv',/get_lun
line='' & readf,tlun,line
nline=0
while ~eof(tlun) do begin
  readf,tlun,line
  if strtrim(line,2) eq '' then continue
  f = strsplit(line,',',/extract)
  c_ion[nline]=f[0] & c_wvl[nline]=double(f[1]) & nline=nline+1
endwhile
free_lun,tlun

openw,lun,outdir+'/idl_box_spectra.txt',/get_lun
printf,lun,'# ion wvl_target nl npix ; then nl rows of: wavelength intensity error intensity_noshift'
printf,lun,'# box y['+strtrim(y0,2)+':'+strtrim(y1,2)+'] x['+strtrim(x0,2)+':'+strtrim(x1,2)+']'

for k=0,nline-1 do begin
  wd = eis_getwindata(l1, c_wvl[k], /refill, /quiet)
  nl = wd.nl & wref = wd.wvl
  acc = dblarr(nl) & accv = dblarr(nl) & nacc = 0L
  for jj=y0,y1 do for ii=x0,x1 do begin
    sp = reform(wd.int[*,ii,jj]) & er = reform(wd.err[*,ii,jj])
    ok = where(finite(sp) and sp gt -1e29, nok)
    if nok lt nl/2 then continue
    dl = 0.0d
    if tag_exist(wd,'WAVE_CORR') then dl = wd.wave_corr[ii,jj]
    acc  = acc  + interpol(sp, wref-dl, wref)
    accv = accv + interpol(er^2, wref-dl, wref)
    nacc = nacc + 1L
  endfor
  spec = acc/nacc & esp = sqrt(accv)/nacc

  ;; ★ 比較用: wave_corr のシフトを **かけずに** 素朴に平均したもの。
  ;;   EIS の分光サンプリングは 0.0223 A / 画素、線幅は sigma ~ 0.030 A なので
  ;;   sigma あたり 1.3 点しか無い。画素ごとに異なるシフトをかけて線形内挿すると
  ;;   エイリアシングで偽の谷が立つ恐れがある（Fe XII 195.119 で疑われた）。
  acc2 = dblarr(nl) & n2 = 0L
  for jj=y0,y1 do for ii=x0,x1 do begin
    sp2 = reform(wd.int[*,ii,jj])
    if total(finite(sp2)) lt nl/2 then continue
    acc2 = acc2 + sp2 & n2 = n2 + 1L
  endfor
  spec_ns = acc2/n2
  printf,lun,c_ion[k]+' '+string(c_wvl[k],format='(f9.4)')+' '+strtrim(nl,2)+' '+strtrim(nacc,2)
  for i=0,nl-1 do printf,lun,wref[i],spec[i],esp[i],spec_ns[i],format='(f12.5,3e16.7)'
  print,'>>> '+c_ion[k]+' nl='+strtrim(nl,2)+' npix='+strtrim(nacc,2)+ $
        ' peak='+string(max(spec),format='(e11.4)')
endfor
free_lun,lun

print,'>>> wrote '+outdir+'/idl_box_spectra.txt'
print,'>>> DONE'
exit
end
