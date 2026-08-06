; AIA 94 Å の **Fe XVIII だけの**温度応答を作る（論文 p.6 の要求）。
;
;   OUTDIR=work bash scripts/idl/run_sswidl.sh scripts/idl/15_aia94_fe18_resp.pro logs/15.log
;
; 論文 p.6:
;   "To utilize the subtracted AIA Fe XVIII intensities, we have computed a
;    new response for this channel that only contains contributions from
;    Fe XVIII, which would contribute to the inversion. The response
;    distributed with the official AIA software contains contributions from
;    several of the known emission lines formed at lower temperatures."
;
; つまり公式の AIA 94 応答は使えない。Fe XVIII の輝線だけを
; AIA 94 の有効面積に畳み込んで作り直す。
;
; R(T) = sum_lines [ Ab(Fe) * ioneq(Fe XVIII,T) * (N_H/N_e) * eps_line(T)/N_e ]
;        * A_eff(lambda_line) / (hc/lambda_line) * (電子/DN) * (画素立体角) / (4 pi)
;
; 単位: [DN cm^5 s^-1 pix^-1]  →  DN/s = int R(T) n_e n_H dh
;
; AIA の有効面積は aiapy (Python) から出すのが確実なので、
; ここでは **Fe XVIII の G(T) と波長リストだけ**を書き出し、
; 畳み込みは scripts/aia94_fe18_response.py で行う。

set_plot,'z'
!quiet = 1

outdir = getenv('OUTDIR')
dlogt  = double(getenv('DLOGT'))
if dlogt le 0 then dlogt = 0.10d
nT   = fix((8.0d - 5.0d)/dlogt) + 1
logT = 5.0d + dlogt*dindgen(nT)
logNe = 9.0d
nhne  = 0.83d

abfile = concat_dir(concat_dir(!xuvtop,'abundance'),'sun_coronal_1992_feldman.abund')
read_abund, abfile, abund, abund_ref
read_ioneq, !ioneq_file, ioneq_logt, ioneq, ioneq_ref

;; AIA 94 のバンドはおよそ 90-98 A。Fe XVIII (Z=26, ion=18) の輝線を全部拾う。
em = emiss_calc(26, 18, temp=logT, dens=logNe, /quiet)
g  = where(em.lambda gt 88.0 and em.lambda lt 100.0, ng)
print,'>>> Fe XVIII lines in 88-100 A : '+strtrim(ng,2)

f = interpol(reform(ioneq[*, 25, 17]), ioneq_logt, logT) > 0.d   ; Fe(Z=26) XVIII(ion=18)
mx = max(f, im)
print,'>>> Fe XVIII ioneq peak at logT = '+string(logT[im],format='(f5.2)')+ $
      '  (frac='+string(mx,format='(f7.4)')+')'

openw,u,outdir+'/fe18_gofnt.txt',/get_lun
printf,u,'# Fe XVIII G(T) for the AIA 94 A band. CHIANTI '+!xuvtop
printf,u,'# abund = sun_coronal_1992_feldman  ioneq = chianti.ioneq  logNe = '+ $
      string(logNe,format='(f4.1)')
printf,u,'# G = Ab(Fe) * ioneq * (N_H/N_e=0.83) * emiss_calc / N_e   [erg cm^3 s^-1]
printf,u,'# nT nline'
printf,u,nT, ng
printf,u,'# logT'
printf,u,logT
printf,u,'# wavelengths [A]'
printf,u,em[g].lambda
printf,u,'# G(T): nline blocks of nT'
for i=0,ng-1 do begin
  gg = abund[25] * f * nhne * reform(em[g[i]].em) / (10.d0^logNe)
  printf,u,gg
  pk = max(gg, ip)
  if pk gt 0 then print,'>>>   '+string(em[g[i]].lambda,format='(f9.4)')+ $
      '  Gmax='+string(pk,format='(e11.4)')+'  at logT='+string(logT[ip],format='(f5.2)')
endfor
free_lun,u

print,'>>> wrote '+outdir+'/fe18_gofnt.txt'
print,'>>> DONE'
exit
end
