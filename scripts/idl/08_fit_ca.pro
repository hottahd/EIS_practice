; Ca 線を Warren+2012 p.6 の拘束どおりにフィットし直す。
;
;   EIS_L1=... BOX=244,274,32,40 OUTDIR=work \
;     bash scripts/idl/run_sswidl.sh scripts/idl/08_fit_ca.pro logs/08.log
;
; 論文 p.6:
;   "To ensure consistency between the fits to the Ca lines we use the widths
;    measured for the Ca XIV 193.874 A and Ca XV 200.972 A lines to constrain
;    the fits to the other Ca lines. The width of the Ca XVI 208.604 A is set
;    equal to that of Ca XV 200.972 A. The width of Ca XVII 192.858 A is
;    limited to be within 0.05 mA of the width of Ca XIV 193.874 A."
;   "The Ca XVII 192.858 A line is blended with Fe XI 192.813 A and a complex
;    of O V lines. We use the method outlined by Ko et al. (2009) to disentangle
;    this blend."
;
; Ko et al. (2009) 相当の実装（原子データは CHIANTI 9.0.1 から取得済み、
; scripts/idl/05b_ca17_ratios.pro の出力）:
;   - O V 多重線は 1 つの振幅スケールと 1 つの共通シフトで動かす。
;     成分間の相対強度は固定（光子数比、192.904 を 1 として）:
;       192.750:0.17484  192.797:0.32722  192.801:0.13115
;       192.904:1.00000  192.911:0.10892  192.915:0.00877
;   - Fe XI 192.813 の強度は同じ箱で測った Fe XI 188.216 から固定する。
;     I(192.813)/I(188.216) = 0.21406 * (188.216/192.813) = 0.20896  (energy 単位)
;     ※ この比は logT 6.0-6.3 で変化しない（logNe=9 固定）
;   - Ca XVII の幅は Ca XIV 193.874 の幅に固定
;   - 背景は定数

; ---------------------------------------------------------------- 汎用モデル
function ca_model, x, p
  common ca_cb, ngauss
  y = x*0.0d + p[3*ngauss]
  for i=0,ngauss-1 do y = y + p[3*i]*exp(-0.5d*((x-p[3*i+1])/p[3*i+2])^2)
  return, y
end

; -------------------------------------------- Ca XVII ブレンド専用モデル
;  p = [ ov_scale, ov_shift, feXI_amp(固定), ca17_amp, ca17_cen, backg ]
;  ca17_wid と feXI_amp は common 経由で外から与える（mpfit の自由度から外す）
function ca17_model, x, p
  common ca17_cb, ov_wvl, ov_rel, ca17_wid, fe11_amp, ov_wid
  y = x*0.0d + p[5]
  for i=0,n_elements(ov_wvl)-1 do $
    y = y + p[0]*ov_rel[i]*exp(-0.5d*((x-(ov_wvl[i]+p[1]))/ov_wid)^2)
  y = y + p[2]*exp(-0.5d*((x-192.813d - p[1])/ov_wid)^2)     ; Fe XI 192.813
  y = y + p[3]*exp(-0.5d*((x-p[4])/ca17_wid)^2)              ; Ca XVII
  return, y
end

; ---------------------------------------------------------------- main
common ca_cb, ngauss
common ca17_cb, ov_wvl, ov_rel, ca17_wid, fe11_amp, ov_wid

set_plot,'z'
!quiet = 1

l1     = getenv('EIS_L1')
outdir = getenv('OUTDIR')
box    = fix(strsplit(getenv('BOX'),',',/extract))
y0 = box[0] & y1 = box[1]-1 & x0 = box[2] & x1 = box[3]-1
print,'>>> l1 = '+l1
print,'>>> box y['+strtrim(y0,2)+':'+strtrim(y1,2)+'] x['+strtrim(x0,2)+':'+strtrim(x1,2)+']'

;; ---- 箱内平均スペクトルを返すヘルパ（インラインで書く）
;; （関数にすると common の扱いが面倒なので、必要な窓ごとに素直に書く）

;; =============== 1) Fe XI 188.216 と Ca XIV 193.874 を測る ===============
;; 使う量: I(Fe XI 188.216) と width(Ca XIV 193.874)

nwin = 5
wtar = [188.216d, 193.874d, 200.972d, 208.604d, 192.858d]
wlo  = [188.016d, 193.600d, 200.810d, 208.404d, 192.600d]
whi  = [188.390d, 194.000d, 201.280d, 208.940d, 193.050d]
ncmp = [2,        2,        2,        2,        0]
cn   = dblarr(3,5)
cn[0,0]=188.217d & cn[1,0]=188.303d
cn[0,1]=193.717d & cn[1,1]=193.873d
cn[0,2]=200.988d & cn[1,2]=201.119d
cn[0,3]=208.598d & cn[1,3]=208.714d
itg  = [0, 1, 0, 0, 0]                 ; 目的の線がどの成分か

res_I = dblarr(5) & res_W = dblarr(5) & res_C = dblarr(5) & res_X = dblarr(5)

ov_wvl = [192.7500d, 192.7970d, 192.8010d, 192.9040d, 192.9110d, 192.9150d]
ov_rel = [0.17484d,  0.32722d,  0.13115d,  1.00000d,  0.10892d,  0.00877d]
;; 光子数比 -> エネルギー比（1/lambda を掛ける。192.904 を基準に規格化）
ov_rel = ov_rel * (192.9040d/ov_wvl)
fe11_ratio = 0.21406d * (188.216d/192.813d)     ; energy 比

openw,lun,outdir+'/idl_ca_intensities.csv',/get_lun
printf,lun,'ion,wvl,I_idl,cen_fit,wid_fit,chi2red,note'

for k=0,nwin-1 do begin

  wd = eis_getwindata(l1, wtar[k], /refill, /quiet)
  ;; 欠損はサンプル単位で外す（wd.missing = -100）。wave_corr の内挿はしない。
  ;; 理由は scripts/idl/04_fit_box_tied.pro のコメント参照。
  nl = wd.nl & wref = wd.wvl
  acc = dblarr(nl) & accv = dblarr(nl) & cnt = lonarr(nl)
  for jj=y0,y1 do for ii=x0,x1 do begin
    sp = reform(wd.int[*,ii,jj]) & er = reform(wd.err[*,ii,jj])
    ok = where(finite(sp) and sp gt wd.missing and finite(er) and er gt wd.missing, nok)
    if nok eq 0 then continue
    acc[ok] = acc[ok] + sp[ok] & accv[ok] = accv[ok] + er[ok]^2 & cnt[ok] = cnt[ok] + 1L
  endfor
  spec = dblarr(nl) & esp = dblarr(nl)
  gd = where(cnt gt 0, ngd)
  spec[gd] = acc[gd]/cnt[gd] & esp[gd] = sqrt(accv[gd])/cnt[gd]

  g = where(wref ge wlo[k] and wref le whi[k] and cnt gt 0, ng)
  x = wref[g] & y = spec[g] & e = esp[g]
  bad = where(e le 0 or ~finite(e), nbad) & if nbad gt 0 then e[bad]=max(e)
  print,'>>> win '+string(wtar[k],format='(f9.3)')+ $
        ' : nl='+strtrim(nl,2)+' [' +string(wref[0],format='(f8.3)')+','+ $
        string(wref[nl-1],format='(f8.3)')+']  npts_in_fit='+strtrim(ng,2)

  if k lt 4 then begin
    ;; ---------- 通常の多成分ガウシアン ----------
    ngauss = ncmp[k]
    npar = 3*ngauss+1
    parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d], tied:''}, npar)
    base = min(y)
    for i=0,ngauss-1 do begin
      parinfo[3*i].value    = (max(y)-base) > 1.0
      parinfo[3*i].limits   = [0.d,1d30]
      parinfo[3*i+1].value  = cn[i,k]
      parinfo[3*i+1].limits = [cn[i,k]-0.06d, cn[i,k]+0.06d]
      parinfo[3*i+2].value  = 0.030d
      parinfo[3*i+2].limits = [0.020d, 0.060d]
      if i gt 0 then parinfo[3*i+2].tied = 'P[2]'
    endfor
    parinfo[3*ngauss].value  = base
    parinfo[3*ngauss].limits = [-1d30,1d30]

    note = 'free width'
    ;; ★ Ca XVI 208.604 は幅を Ca XV 200.972 の幅に固定（論文の指定）
    if k eq 3 then begin
      for i=0,ngauss-1 do begin
        parinfo[3*i+2].value = res_W[2]
        parinfo[3*i+2].fixed = 1
        parinfo[3*i+2].tied  = ''
      endfor
      note = 'width fixed to CaXV = '+string(res_W[2],format='(f7.5)')
    endif

    p = mpfitfun('ca_model', x, y, e, parinfo=parinfo, perror=perror, $
                 bestnorm=bn, dof=dof, /quiet)
    ic = itg[k]
    res_I[k] = p[3*ic]*p[3*ic+2]*sqrt(2*!dpi)
    res_C[k] = p[3*ic+1]
    res_W[k] = p[3*ic+2]
    res_X[k] = (dof gt 0) ? bn/dof : -1.

  endif else begin
    ;; ---------- Ca XVII: Ko et al. (2009) 相当 ----------
    ca17_wid = res_W[1]                       ; Ca XIV の幅に固定
    ov_wid   = res_W[1]                       ; O V / Fe XI も同じ装置幅
    fe11_amp = res_I[0]*fe11_ratio/(ov_wid*sqrt(2*!dpi))   ; Fe XI 188.216 から固定

    npar = 6
    parinfo = replicate({value:0.d, fixed:0, limited:[1,1], limits:[0.d,0.d], tied:''}, npar)
    ;; ★ 注意: replicate の既定が limited=[1,1], limits=[0,0] なので、
    ;;   固定パラメータでも limited を明示的に落とさないと
    ;;   「初期値が limits の外」で MPFIT がその場で失敗し、
    ;;   **初期値をそのまま返す**（しかもエラーを出さない）。
    ;;   最初これに気づかず Ca XVII の値を誤って読んだ。
    parinfo[0].value = max(y)*0.3d & parinfo[0].limits = [0.d,1d30]     ; O V スケール
    parinfo[1].value = 0.d         & parinfo[1].limits = [-0.05d,0.05d] ; 共通シフト
    parinfo[2].value = fe11_amp    & parinfo[2].fixed  = 1              ; Fe XI 固定
    parinfo[2].limited = [0,0]     & parinfo[2].limits = [0.d,0.d]
    parinfo[3].value = max(y)*0.3d & parinfo[3].limits = [0.d,1d30]     ; Ca XVII 振幅
    parinfo[4].value = 192.853d    & parinfo[4].limits = [192.82d,192.89d]
    parinfo[5].value = min(y)      & parinfo[5].limited = [0,0]

    p = mpfitfun('ca17_model', x, y, e, parinfo=parinfo, perror=perror, $
                 bestnorm=bn, dof=dof, status=st, errmsg=emsg, /quiet)
    print,'>>>   mpfit status='+strtrim(st,2)+' dof='+strtrim(dof,2)+ $
          ' npar_free='+strtrim(n_elements(where(parinfo.fixed eq 0)),2)+ $
          ((emsg ne '') ? ('  errmsg='+emsg) : '')
    print,'>>>   p = '+strjoin(string(p,format='(e12.4)'),' ')
    res_I[k] = p[3]*ca17_wid*sqrt(2*!dpi)
    res_C[k] = p[4]
    res_W[k] = ca17_wid
    res_X[k] = (dof gt 0) ? bn/dof : -1.
    note = 'Ko+2009: OV tied, FeXI fixed from 188.216, width=CaXIV'
    print,'>>>   O V scale='+string(p[0],format='(e11.4)')+ $
          '  shift='+string(p[1],format='(f8.4)')+ $
          '  I(FeXI 192.813)='+string(fe11_amp*ov_wid*sqrt(2*!dpi),format='(f9.2)')
  endelse

  printf,lun,string(wtar[k],format='(f9.3)')+','+string(res_I[k],format='(e12.5)')+ $
         ','+string(res_C[k],format='(f9.4)')+','+string(res_W[k],format='(f8.5)')+ $
         ','+string(res_X[k],format='(f10.3)')+','+note
  print,'>>> '+string(wtar[k],format='(f9.3)')+'  I='+string(res_I[k],format='(e11.4)')+ $
        '  cen='+string(res_C[k],format='(f9.4)')+'  wid='+string(res_W[k],format='(f7.5)')+ $
        '  chi2r='+string(res_X[k],format='(f9.2)')+'  ('+note+')'

endfor
free_lun,lun

print,'>>> --- 論文 Table 2 との比 ---'
pap = [498.81d, 182.64d, 127.92d, 31.12d, 41.75d]
for k=0,4 do print,'>>>   '+string(wtar[k],format='(f9.3)')+ $
      '  ratio='+string(res_I[k]/pap[k],format='(f8.3)')

print,'>>> DONE'
exit
end
