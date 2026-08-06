; Warren+2012 の 22 輝線の寄与関数 G(T) を CHIANTI から作る。
;
;   OUTDIR=work bash scripts/idl/run_sswidl.sh scripts/idl/09_gofnt.pro logs/09.log
;   （SSW_INSTR に chianti が必要）
;
; 論文 p.6 の指定に合わせる:
;   - 原子データ: CHIANTI（論文は v7、このマシンは **v9.0.1**。ここだけ違う）
;   - 組成: Feldman et al. (1992) コロナ組成 -> sun_coronal_1992_feldman.abund
;   - 電離平衡: Dere et al. (2009) -> chianti.ioneq（v6 以降変わっていない）
;
; 定義:
;   G(T) = Ab(Z) * ioneq(Z,ion,T) * (N_H/N_e) * eps(T) / N_e
;   ★ emiss_calc が返すのは hc/lambda * N_j * A_ji **だけ**（N_e で割っていない）。
;     emiss_calc.pro のヘッダに明記されている。N_j は「そのイオンに対する
;     上準位の存在比」。したがって寄与関数にするには **N_e で割る**必要がある。
;     最初これを忘れて G が 1e9 倍大きくなり、DEM の初期値が 1e-6 という
;     あり得ない値になった（正しくは ~1e27 cm^-5）。
;   これで I = (1/4pi) * int G(T) * N_e N_H ds/dlogT dlogT となる。
;   mcmc_dem には Z=1, /noph で「組成込みの放射率」として渡す。

set_plot,'z'
!quiet = 1

outdir = getenv('OUTDIR')
;; ★ 温度ビンの幅は DLOGT で変えられる。輝線が 22 本しか無いので
;;   0.05 dex (61 ビン) だと劣決定になり DEM が隣同士で 1-2 桁跳ねる。
;;   PINTofALE の AIA スレッドも既定は 0.1 dex。
dlogt = double(getenv('DLOGT'))
if dlogt le 0 then dlogt = 0.10d
nT0   = fix((8.0d - 5.0d)/dlogt) + 1
logT  = 5.0d + dlogt*dindgen(nT0)          ; 5.0 - 8.0
nT     = n_elements(logT)
print,'>>> dlogT = '+string(dlogt,format='(f5.3)')
logNe  = 9.0d                              ; 活動領域コアの典型値
nhne   = 0.83d                             ; N_H/N_e

abfile = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_1992_feldman.abund')
ioneqf = !ioneq_file

print,'>>> CHIANTI  '+!xuvtop
print,'>>> abund    '+abfile
print,'>>> ioneq    '+ioneqf
print,'>>> log Ne   '+string(logNe,format='(f4.1)')
print,'>>> logT     '+string(logT[0],format='(f4.2)')+' .. '+ $
                      string(logT[nT-1],format='(f4.2)')+'  ('+strtrim(nT,2)+' points)'

read_abund, abfile, abund, abund_ref
read_ioneq, ioneqf, ioneq_logt, ioneq, ioneq_ref

;; --- 22 輝線: ion名, Z, ion段, 波長
nline = 22
ionnm = ['Si VII','Fe IX','Fe IX','Fe X','Fe XI','Fe XI','S X','Si X', $
         'Fe XII','Fe XII','Fe XIII','Fe XIII','Fe XIV','Fe XIV','Fe XV', $
         'S XIII','Fe XVI','Ar XIV','Ca XIV','Ca XV','Ca XVI','Ca XVII']
zz    = [14, 26, 26, 26, 26, 26, 16, 14, 26, 26, 26, 26, 26, 26, 26, $
         16, 26, 18, 20, 20, 20, 20]
ii    = [ 7,  9,  9, 10, 11, 11, 10, 10, 12, 12, 13, 13, 14, 14, 15, $
         13, 16, 14, 14, 15, 16, 17]
wv    = [275.368d, 188.497d, 197.862d, 184.536d, 180.401d, 188.216d, $
         264.233d, 258.375d, 192.394d, 195.119d, 202.044d, 203.826d, $
         264.787d, 270.519d, 284.160d, 256.686d, 262.984d, 194.396d, $
         193.874d, 200.972d, 208.604d, 192.858d]

gofT = dblarr(nT, nline)
wfit = dblarr(nline)

zprev = -1 & iprev = -1
for k=0,nline-1 do begin

  if (zz[k] ne zprev) or (ii[k] ne iprev) then begin
    em = emiss_calc(zz[k], ii[k], temp=logT, dens=logNe, /quiet)
    zprev = zz[k] & iprev = ii[k]
  endif

  ;; ★ 波長が最も近い遷移を選ぶだけでは駄目。CHIANTI には同じ波長付近に
  ;;   非常に弱い遷移が何本も入っていて、そちらを拾うと G(T) が 6 桁小さくなる
  ;;   （最初 Fe IX 188.497 / 197.862 と Fe XVI 262.984 でこれを踏んだ）。
  ;;   ±0.03 A の中で **放射率が最大** のものを選ぶ。
  cand = where(abs(em.lambda - wv[k]) lt 0.03d, ncand)
  if ncand eq 0 then begin
    d = min(abs(em.lambda - wv[k]), j)
    print,'>>> WARNING '+ionnm[k]+' no line within 0.03 A; nearest dl='+string(d,format='(f6.3)')
  endif else begin
    best = -1.d & j = cand[0]
    for c=0,ncand-1 do begin
      pk = max(em[cand[c]].em)
      if pk gt best then begin & best = pk & j = cand[c] & endif
    endfor
    if ncand gt 1 then print,'>>>   ('+ionnm[k]+': '+strtrim(ncand,2)+ $
        ' transitions within 0.03 A, took the strongest)'
  endelse
  wfit[k] = em[j].lambda

  ;; イオン分布を logT グリッドへ
  f = interpol(reform(ioneq[*, zz[k]-1, ii[k]-1]), ioneq_logt, logT) > 0.d
  gofT[*,k] = abund[zz[k]-1] * f * nhne * reform(em[j].em) / (10.d0^logNe)

  mx = max(gofT[*,k], im)
  print,'>>> '+string(ionnm[k],format='(a-9)')+string(wv[k],format='(f9.3)')+ $
        ' -> '+string(wfit[k],format='(f9.3)')+ $
        '  Gmax='+string(mx,format='(e11.4)')+ $
        '  at logT='+string(logT[im],format='(f5.2)')
endfor

;; --- 書き出し（Python 側でも使えるプレーンテキスト）
gname = getenv('GOFNT_OUT')
if gname eq '' then gname = 'gofnt_chianti901.txt'
openw,lun,outdir+'/'+gname,/get_lun
printf,lun,'# G(T) for Warren+2012 lines. CHIANTI '+!xuvtop
printf,lun,'# abund = sun_coronal_1992_feldman   ioneq = chianti.ioneq   logNe = '+ $
       string(logNe,format='(f4.1)')
printf,lun,'# G = Ab(Z) * ioneq * (N_H/N_e=0.83) * emiss_calc  [erg cm^3 s^-1]
printf,lun,'# nT nline'
printf,lun,nT, nline
printf,lun,'# logT'
printf,lun,logT
printf,lun,'# ion wvl_target wvl_chianti'
for k=0,nline-1 do printf,lun,ionnm[k]+' '+string(wv[k],format='(f9.4)')+' '+ $
       string(wfit[k],format='(f9.4)')
printf,lun,'# G(T) : nline blocks of nT values'
for k=0,nline-1 do printf,lun,gofT[*,k]
free_lun,lun

print,'>>> wrote '+outdir+'/'+gname
print,'>>> DONE'
exit
end
