; 論文と同じ手法 (PINTofALE の MCMC_DEM, Kashyap & Drake 1998) で DEM を出す。
;
;   PINTofALE=/scr/a000/c0234hotta/PINTofALE OUTDIR=work NSIM=2000 \
;     bash scripts/idl/run_sswidl.sh scripts/idl/13_mcmc_dem.pro logs/13.log
;
; 使うもの:
;   work/gofnt_chianti901.txt      ... scripts/idl/09_gofnt.pro が作った G(T)
;   work/idl_intensities_tied.csv  ... scripts/idl/04_fit_box_tied.pro の輝線強度
;   work/idl_ca_intensities.csv    ... scripts/idl/08_fit_ca.pro の Ca 線（論文の拘束つき）
;
; ★ 単位のつじつま（ここを外すと DEM の絶対値が 4pi や 0.83 だけずれる）
;   - PINTofALE の lineflx は入力放射率を [1e-23 erg cm^3/s] と仮定する
;     → G [erg cm^3/s] に 1e23 を掛ける
;   - lineflx は **4pi で割らない**（lineflx.pro の 256 行目に明記）。
;     観測強度は erg cm^-2 s^-1 sr^-1 なので、**観測側に 4pi を掛けて**渡す。
;   - lineflx は既定で N_H/N_e = 0.83 を掛ける。G(T) 側でも掛けてあるので
;     二重適用を避けるため **nhne=1.0** を渡す。
;   - 組成も G(T) に畳み込み済みなので **Z=1**（水素に見せかける）＋ /noabund。
;   - /noph: 光子数に変換させない（G はエネルギー単位）。
;
; 出力の DEM は [cm^-5 / logK]。論文 Figure 6-8 の xi(Te)dTe と同じ土俵。

set_plot,'z'
!quiet = 1

poa    = getenv('PINTofALE')
outdir = getenv('OUTDIR')
nsim   = long(getenv('NSIM')) > 500L
print,'>>> PINTofALE '+poa
;; パスは run_sswidl.sh が IDL 起動前に IDL_PATH へ入れている。
;; ここで !path を書き換えても、.run は既に一括コンパイル済みなので手遅れ。
print,'>>> mcmc_dem -> '+file_which(!path,'mcmc_dem.pro')
print,'>>> nsim = '+strtrim(nsim,2)

;; ===================== G(T) を読む =====================
openr,u,outdir+'/gofnt_chianti901.txt',/get_lun
line=''
repeat readf,u,line until strmid(line,0,6) eq '# nT n'
nT=0L & nline=0L & readf,u,nT,nline
readf,u,line                                   ; '# logT'
logT = dblarr(nT) & readf,u,logT
readf,u,line                                   ; '# ion ...'
ionnm = strarr(nline) & wv = dblarr(nline)
for k=0,nline-1 do begin
  readf,u,line
  f = strsplit(line,' ',/extract)
  ionnm[k] = f[0]+' '+f[1]
  wv[k]    = double(f[2])
endfor
readf,u,line                                   ; '# G(T) ...'
gofT = dblarr(nT, nline)
tmp  = dblarr(nT)
for k=0,nline-1 do begin & readf,u,tmp & gofT[*,k] = tmp & endfor
free_lun,u
print,'>>> G(T): nT='+strtrim(nT,2)+' nline='+strtrim(nline,2)

;; ★ 輝線が拘束しない温度域まで解を許すと、そこに EM が漏れて
;;   高温側の傾き beta が出鱈目になる（PINTofALE の AIA スレッドが
;;   "toothpaste tube effect" と呼ぶ現象）。使う輝線の感度範囲に絞る。
tlo = double(getenv('TLO')) & thi = double(getenv('THI'))
if tlo le 0 then tlo = 5.5d
if thi le 0 then thi = 7.5d
keep = where(logT ge tlo and logT le thi, nT)
logT = logT[keep] & gofT = gofT[keep,*]
print,'>>> logT を '+string(tlo,format='(f4.2)')+'-'+string(thi,format='(f4.2)')+ $
      ' に制限 (nT='+strtrim(nT,2)+')'

;; ===================== 観測強度を読む =====================
;; 04_fit_box_tied.pro の CSV（ion,wvl,...,I_idl,...,I_paper,sig_paper,ratio）
iobs = dblarr(nline) & esig = dblarr(nline) & ipap = dblarr(nline)
openr,u,outdir+'/idl_intensities_tied.csv',/get_lun
readf,u,line
while ~eof(u) do begin
  readf,u,line
  f = strsplit(line,',',/extract)
  w = double(f[1])
  d = min(abs(wv-w), k)
  if d lt 0.01 then begin
    iobs[k] = double(f[5])
    ipap[k] = double(f[10])
    esig[k] = double(f[11])/double(f[10]) * double(f[5])   ; 論文と同じ相対誤差
  endif
endwhile
free_lun,u

;; Ca 線は論文の拘束つきの値で上書きする（08_fit_ca.pro）
openr,u,outdir+'/idl_ca_intensities.csv',/get_lun
readf,u,line
while ~eof(u) do begin
  readf,u,line
  f = strsplit(line,',',/extract)
  w = double(f[0])
  d = min(abs(wv-w), k)
  ;; 208.604 と 192.858 だけ差し替え（他は同じ値）
  if d lt 0.01 and (abs(w-208.604d) lt 0.01 or abs(w-192.858d) lt 0.01) then begin
    iobs[k] = double(f[1])
    esig[k] = 0.22d * iobs[k]
    print,'>>> replaced '+ionnm[k]+' with Ca-constrained value '+string(iobs[k],format='(f10.3)')
  endif
endwhile
free_lun,u

;; ★ AIA 94 A の Fe XVIII を 23 本目の拘束として加える（論文 Table 2 の最終行）。
;;   我々の最高温の EIS 輝線は Ca XVII (logT 6.75) までで、
;;   それより上が拘束されず EM が漏れて beta が出鱈目になる。
;;   Fe XVIII は logT 6.90 にピークがあるのでちょうど効く。
;;   単位: AIA の応答 R(T) は [DN cm^5 s^-1 pix^-1] で 1/(4pi) と画素立体角
;;   込みなので、観測値 (DN/s) に 4pi を掛けてはいけない。
aia_obs = double(getenv('AIA_FE18'))       ; [DN/s]
use_aia = 0
if aia_obs gt 0 and file_test(outdir+'/aia94_fe18_response.txt') then begin
  nr = file_lines(outdir+'/aia94_fe18_response.txt')
  openr,u,outdir+'/aia94_fe18_response.txt',/get_lun
  aline='' & atl=dblarr(nr) & arr=dblarr(nr) & na=0L
  while ~eof(u) do begin
    readf,u,aline
    if strmid(strtrim(aline,2),0,1) eq '#' then continue
    f = strsplit(aline,' ',/extract)
    atl[na]=double(f[0]) & arr[na]=double(f[1]) & na=na+1
  endwhile
  free_lun,u
  aresp = interpol(arr[0:na-1], atl[0:na-1], logT) > 0.d
  use_aia = 1
  print,'>>> AIA Fe XVIII 拘束を追加: I='+string(aia_obs,format='(f8.3)')+' DN/s'
  mx = max(aresp, im)
  print,'>>>   R(T) peak = '+string(mx,format='(e11.4)')+' at logT='+ $
        string(logT[im],format='(f5.2)')
endif

ok = where(iobs gt 0, nok)
print,'>>> usable lines: '+strtrim(nok,2)+' / '+strtrim(nline,2)
for i=0,nok-1 do print,'>>>   '+string(ionnm[ok[i]],format='(a-14)')+ $
    ' I='+string(iobs[ok[i]],format='(e11.4)')+' +-'+string(esig[ok[i]],format='(e11.4)')+ $
    '  (paper '+string(ipap[ok[i]],format='(f9.2)')+')'

emis = gofT[*,ok] * 1d23                ; [1e-23 erg cm^3/s]
wvl  = wv[ok]
flx  = iobs[ok] * 4.d0*!dpi             ; erg cm^-2 s^-1 sr^-1 -> erg cm^-2 s^-1
efl  = esig[ok] * 4.d0*!dpi
zz   = intarr(nok) + 1                  ; 組成は emis に畳み込み済み

if use_aia then begin
  emis = [[emis], [aresp*1d23]]
  wvl  = [wvl, 93.932d]
  flx  = [flx, aia_obs]                  ; ★ 4pi を掛けない
  efl  = [efl, 0.19d*aia_obs]            ; 論文 Table 2 の 1.40/7.20
  zz   = [zz, 1]
  nok  = nok + 1
  ionnm = [ionnm, 'AIA 94 FeXVIII']
  ok    = [ok, n_elements(ionnm)-1]
  ipap  = [ipap, 7.20d]
endif

;; ===================== 初期値・探索範囲・平滑化 =====================
;; EM loci の最小値から出発する（PINTofALE の AIA スレッドと同じやり方）
emc = dblarr(nT, nok)
;; ★ emis は [1e-23 erg cm^3/s] なので、EM loci を作るときは 1e-23 を戻すこと。
;;   忘れると初期 DEM が 1e23 倍ずれる。
for i=0,nok-1 do emc[*,i] = flx[i] / ((emis[*,i]*1d-23) > 1d-40)
ymin = min(emc[where(finite(emc) and emc gt 0)])
diffem = dblarr(nT) + ymin
demrng = dblarr(nT,2) & demrng[*,0] = diffem/1d5 & demrng[*,1] = diffem*1d10
print,'>>> initial DEM = '+string(ymin,format='(e11.4)')+' [cm^-5/logK]'

sampenv = dblarr(nT)
for i=0,nT-1 do sampenv[i] = stddev(emis[i,*])
stepenv = dblarr(nT)
for i=0,nok-1 do stepenv = stepenv + (emis[*,i] gt 1d-5*max(emis[*,i]))
sampenv = sampenv * stepenv
smooscl = findscale(sampenv)
;; nutilde (実効自由度) が 1 を超えるまで平滑化スケールを縮める。
;; PINTofALE の AIA スレッドが「nu>1 を必ず確認せよ」と強調している点。
for itry=0,6 do begin
  jnk = varsmooth(sampenv, smooscl, nueff=nueff, nutilde=nutilde)
  print,'>>> smooscl try '+strtrim(itry,2)+': nueff='+string(nueff,format='(f8.2)')+ $
        '  nutilde='+string(nutilde,format='(f8.2)')
  if nutilde gt 1.0 then break
  smooscl = smooscl/2.0
endfor

ulim = intarr(nok)

;; ===================== MCMC =====================
t0 = systime(1)
dem = mcmc_dem(wvl, flx, emis, Z=zz, logt=logT, diffem=diffem, fsigma=efl, $
               ulim=ulim, nsim=nsim, nburn=long(nsim/5), demrng=demrng, $
               sampenv=sampenv, smooscl=smooscl, /noph, noabund=1, nhne=1.0, $
               savfil=outdir+'/mcmc_dem.save', $
               storpar=storpar, storidx=storidx, simprb=simprb, $
               simdem=simdem, demerr=demerr, simflx=simflx, simprd=simprd, $
               verbose=0)   ; ★ verbose>0 だと window,0 を呼ぶので Z バッファでは落ちる
print,'>>> mcmc_dem took '+string((systime(1)-t0)/60.,format='(f7.2)')+' min'

;; ===================== 結果 =====================
fpred = pred_flx(emis, logT, wvl, dem, Z=zz, /noph, noabund=1, nhne=1.0)
print,'>>> --- 観測 vs モデル（論文 Table 2 の I_obs / I_dem / R に相当）---'
print,'>>> line            I_obs      I_dem        R      I_obs(paper)  I_dem(paper)'
;; AIA の行は 4pi で割ってはいけない（応答に 1/(4pi) が入っているため）
for i=0,nok-1 do begin
  sc = 4*!dpi
  if use_aia and i eq nok-1 then sc = 1.d
  print,'>>> '+string(ionnm[ok[i]],format='(a-14)')+ $
    string(flx[i]/sc,format='(e12.4)')+string(fpred[i]/sc,format='(e12.4)')+ $
    string(flx[i]/fpred[i],format='(f9.3)')+string(ipap[ok[i]],format='(f14.2)')
endfor

;; ★ mcmc_dem が返す dem は「chi2 最小の 1 実現」なので凸凹する。
;;   論文が示しているのは MCMC アンサンブル（whisker plot）なので、
;;   simdem から温度ごとの中央値と 50% 区間を出す。
sz = size(simdem)
med = dblarr(nT) & q25 = dblarr(nT) & q75 = dblarr(nT)
if sz[0] eq 2 and sz[1] eq nT then begin
  nsam = sz[2]
  print,'>>> simdem: nT='+strtrim(sz[1],2)+' x nsim='+strtrim(nsam,2)
  for i=0,nT-1 do begin
    v = reform(simdem[i,*])
    v = v[sort(v)]
    med[i] = v[nsam/2]
    q25[i] = v[long(0.25*nsam)]
    q75[i] = v[long(0.75*nsam)]
  endfor
endif else begin
  print,'>>> simdem の形が想定外: '+strjoin(strtrim(sz,2),' ')
  med = dem
endelse

openw,u,outdir+'/mcmc_dem_result.txt',/get_lun
printf,u,'# logT  DEM_best  DEM_median  DEM_q25  DEM_q75   [cm^-5/logK]'
for i=0,nT-1 do printf,u,logT[i],dem[i],med[i],q25[i],q75[i],format='(f7.3,4e15.6)'
free_lun,u
print,'>>> wrote '+outdir+'/mcmc_dem_result.txt'

print,'>>> DONE'
exit
end
