; Warren+2012 Table 2 の 22 輝線を SSWIDL で出す。
;
;   EIS_L1=... BOX="ymin,ymax,xmin,xmax" OUTDIR=work \
;     bash scripts/idl/run_sswidl.sh scripts/idl/03_fit_box.pro logs/03.log
;
; 手順は論文 (§2) と同じ「箱の中で平均 → その 1 本のスペクトルをフィット」。
; フィット窓・ガウシアン成分数・初期中心波長は **eispac のテンプレートと同一**
; にしてある（work/eispac_templates.csv、Python 側で書き出したもの）。
; そうしないと「IDL と Python の差」なのか「モデルの差」なのか分からなくなる。
;
; ★ 箱は **配列インデックス**で指定する（eispac と同じ y0:y1, x0:x1）。
;   検証済み: eis_getwindata の int は [nl,nx,ny] で、行・列の並びは
;   eispac の (ny,nx) と完全に一致する（相互相関のラグ 0、相関 0.997）。
;   ただし solar_y は eispac の方が +16.16 arcsec ずれている
;   （eispac が EIS→AIA のポインティング補正を入れているため。solar_x は一致）。
;   したがって solar 座標で箱を切ると別の場所を測ってしまう。インデックスで切ること。

; ---------------------------------------------------------------- model
function w12_model, x, p
  common w12_cb, ngauss
  y = x*0.0d + p[3*ngauss]                      ; 定数背景 (eispac の n_poly=1 と同じ)
  for i=0,ngauss-1 do y = y + p[3*i]*exp(-0.5d*((x-p[3*i+1])/p[3*i+2])^2)
  return, y
end

; ---------------------------------------------------------------- main
common w12_cb, ngauss

set_plot,'z'
!quiet = 1

l1     = getenv('EIS_L1')
outdir = getenv('OUTDIR')
box    = fix(strsplit(getenv('BOX'),',',/extract))   ; y0,y1,x0,x1 (Python 流の半開区間)
y0 = box[0] & y1 = box[1]-1 & x0 = box[2] & x1 = box[3]-1

print,'>>> l1  = '+l1
print,'>>> box index y['+strtrim(y0,2)+':'+strtrim(y1,2)+'] x['+strtrim(x0,2)+':'+strtrim(x1,2)+']'

;; テンプレート表を読む。read_csv のタグ名 (FIELD1 か FIELD01 か) は
;; IDL のバージョンで揺れるので、自前でパースする。
nmax = 200
c_ion=strarr(nmax) & c_wvl=dblarr(nmax) & c_ng=intarr(nmax)
c_wmin=dblarr(nmax) & c_wmax=dblarr(nmax) & c_cen=dblarr(3,nmax)
c_wid=dblarr(nmax) & c_ip=dblarr(nmax) & c_sp=dblarr(nmax)
openr,tlun,outdir+'/eispac_templates.csv',/get_lun
line='' & readf,tlun,line          ; ヘッダ
nline=0
while ~eof(tlun) do begin
  readf,tlun,line
  if strtrim(line,2) eq '' then continue
  f = strsplit(line,',',/extract)
  c_ion[nline]  = f[0]
  c_wvl[nline]  = double(f[1])
  c_ng[nline]   = fix(f[2])
  c_wmin[nline] = double(f[4])
  c_wmax[nline] = double(f[5])
  c_cen[0,nline]= double(f[6]) & c_cen[1,nline]=double(f[7]) & c_cen[2,nline]=double(f[8])
  c_wid[nline]  = double(f[9])
  c_ip[nline]   = double(f[10])
  c_sp[nline]   = double(f[11])
  nline = nline+1
endwhile
free_lun,tlun
print,'>>> '+strtrim(nline,2)+' lines to fit'

openw,lun,outdir+'/idl_intensities_tied.csv',/get_lun
printf,lun,'ion,wvl,ngauss,icomp,npix,I_idl,eI_idl,cen_fit,wid_fit,chi2red,I_paper,sig_paper,ratio'

for k=0,nline-1 do begin

  ion    = c_ion[k]
  wvl    = c_wvl[k]
  ngauss = c_ng[k]
  wmin   = c_wmin[k]
  wmax   = c_wmax[k]
  cen0   = reform(c_cen[*,k])
  wid0   = c_wid[k]
  ipaper = c_ip[k]
  spaper = c_sp[k]

  wd = eis_getwindata(l1, wvl, /refill, /quiet)

  ;; --- 箱の画素（インデックス指定）
  gx = x0 + indgen(x1-x0+1) & nxg = n_elements(gx)
  gy = y0 + indgen(y1-y0+1) & nyg = n_elements(gy)
  if max(gx) ge wd.nx or max(gy) ge wd.ny then begin
    print,'>>> '+ion+': box outside raster - skipped'
    continue
  endif

  ;; --- 箱内平均スペクトル
  ;; ★ 欠損は **サンプル単位** で外すこと。eis_getwindata は欠損画素に
  ;;   wd.missing (= -100) を入れる。スペクトル単位でしか判定しないと、
  ;;   検出器の不良列（この観測では Fe XII 195.119 窓の 195.1039 A で
  ;;   箱内 240 画素中 176 画素が欠損）の -100 が平均に混ざり、
  ;;   線のコアに偽の谷ができて積分強度が 24% 小さくなる。
  ;; ★ wave_corr のシフト＋線形内挿はしない。EIS のサンプリングは
  ;;   0.0223 A/画素で線幅 sigma ~ 0.030 A（sigma あたり 1.3 点）しかなく、
  ;;   内挿すると壊れた点を隣に塗り広げてしまう。eispac 側も素直に平均している。
  nl   = wd.nl
  wref = wd.wvl
  acc = dblarr(nl) & accv = dblarr(nl) & cnt = lonarr(nl)
  for jj=y0,y1 do for ii=x0,x1 do begin
    sp = reform(wd.int[*,ii,jj]) & er = reform(wd.err[*,ii,jj])
    ok = where(finite(sp) and sp gt wd.missing and finite(er) and er gt wd.missing, nok)
    if nok eq 0 then continue
    acc[ok]  = acc[ok]  + sp[ok]
    accv[ok] = accv[ok] + er[ok]^2
    cnt[ok]  = cnt[ok]  + 1L
  endfor
  nacc = max(cnt)
  gd = where(cnt gt 0, ngd)
  if ngd eq 0 then begin
    print,'>>> no good samples - skipped' & continue
  endif
  spec = dblarr(nl) & esp = dblarr(nl)
  spec[gd] = acc[gd]/cnt[gd]
  esp[gd]  = sqrt(accv[gd])/cnt[gd]
  cntmap = cnt

  ;; --- フィット窓
  g = where(wref ge wmin and wref le wmax and cntmap gt 0, ng)
  x = wref[g] & y = spec[g] & e = esp[g]
  bad = where(e le 0 or ~finite(e), nbad)
  if nbad gt 0 then e[bad] = max(e)

  ;; --- 初期値と拘束
  ;; ★ v1 との違い: **同じ窓の全ガウシアンで線幅を共有する**（tied）。
  ;;   v1 では線幅を成分ごとに自由にしたため、
  ;;     - Fe XII 195.119: 自己ブレンド 195.179 を太い 1 本が飲み込み chi2r=9738
  ;;     - Ca XVI 208.604: 幅 0.067 A まで太って Fe XIII 208.679 を飲み込み ratio=3.11
  ;;   となった。eispac のテンプレートも全成分の初期幅が同一で、
  ;;   Warren+2012 §2 も Ca 線の幅を互いに拘束している。
  ;;   幅の上限は 0.05 A（EIS の装置幅 ~0.056 A FWHM = sigma 0.024 に
  ;;   熱幅を足した程度）。
  npar = 3*ngauss + 1
  p0 = dblarr(npar)
  parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d], tied:''}, npar)
  base = min(y)
  for i=0,ngauss-1 do begin
    p0[3*i]   = (max(y)-base) > 1.0
    p0[3*i+1] = cen0[i]
    p0[3*i+2] = wid0 > 0.02
    parinfo[3*i].limits   = [0.d, 1d30]
    parinfo[3*i+1].limits = [cen0[i]-0.06d, cen0[i]+0.06d]
    parinfo[3*i+2].limits = [0.020d, 0.050d]
    if i gt 0 then parinfo[3*i+2].tied = 'P[2]'      ; 幅を第 0 成分に固定
  endfor
  p0[3*ngauss] = base
  parinfo[3*ngauss].limits = [-1d30, 1d30]
  parinfo.value = p0

  p = mpfitfun('w12_model', x, y, e, parinfo=parinfo, perror=perror, $
               bestnorm=bestnorm, dof=dof, status=st, /quiet)

  ;; --- 目的の輝線に対応する成分を選ぶ（成分順は波長順とは限らない）
  dmin = 1e9 & ic = 0
  for i=0,ngauss-1 do if abs(cen0[i]-wvl) lt dmin then begin & dmin=abs(cen0[i]-wvl) & ic=i & endif

  amp = p[3*ic] & cen = p[3*ic+1] & wid = p[3*ic+2]
  ii   = amp*wid*sqrt(2*!dpi)
  ea   = (n_elements(perror) ge npar) ? perror[3*ic] : 0.d
  ew   = (n_elements(perror) ge npar) ? perror[3*ic+2] : 0.d
  eii  = ii*sqrt((ea/amp)^2 + (ew/wid)^2)
  chi2 = (dof gt 0) ? bestnorm/dof : -1.

  printf,lun,ion+','+string(wvl,format='(f8.3)')+','+strtrim(ngauss,2)+','+strtrim(ic,2)+ $
         ','+strtrim(nacc,2)+','+string(ii,format='(e12.5)')+','+string(eii,format='(e12.5)')+ $
         ','+string(cen,format='(f9.4)')+','+string(wid,format='(f8.5)')+ $
         ','+string(chi2,format='(f10.3)')+','+string(ipaper,format='(f9.2)')+ $
         ','+string(spaper,format='(f8.2)')+','+string(ii/ipaper,format='(f8.3)')

  print,'>>> '+string(ion,format='(a-9)')+string(wvl,format='(f8.3)')+ $
        '  I='+string(ii,format='(e11.4)')+'  paper='+string(ipaper,format='(f9.2)')+ $
        '  ratio='+string(ii/ipaper,format='(f7.3)')+'  chi2r='+string(chi2,format='(f8.2)')+ $
        '  npix='+strtrim(nacc,2)

endfor

free_lun,lun
print,'>>> wrote '+outdir+'/idl_intensities_tied.csv'
print,'>>> DONE'
exit
end
