; Si VII 275.368 の窓に何が居るか CHIANTI 9.0.1 で確認する。
; どの箱でも論文の 0.36-0.45 倍にしかならない唯一の輝線なので、
; モデル化していないブレンドが無いかを潰す。
set_plot,'z'
!quiet=1
logT = 4.5d + 0.05d*dindgen(51)
dens = 1d9
;; Si VII (Z=14, ion=7) の 275 付近
em = emiss_calc(14, 7, temp=logT, dens=alog10(dens), /quiet)
g = where(em.lambda gt 274.8 and em.lambda lt 276.0, n)
print,'>>> Si VII lines 274.8-276.0 : '+strtrim(n,2)
it = (where(abs(logT-5.80) lt 0.001))[0]
if n gt 0 then begin
  ref = max(em[g].em[it])
  for i=0,n-1 do print,'>>>   SiVII '+string(em[g[i]].lambda,format='(f9.4)')+ $
      '  rel='+string(em[g[i]].em[it]/ref,format='(f9.5)')
endif
;; 全イオンで 275.1-275.7 を洗う（時間はかかるが決定的）
ch_synthetic, 275.10, 275.70, output=out, err_msg=err, density=dens, $
              logt_isothermal=[5.5,5.8,6.2,6.5], logem_isothermal=[27.,27.,27.,27.], $
              ioneq_name=!ioneq_file, /all
if err eq '' then begin
  l = out.lines
  s = reverse(sort(l.int))
  nn = n_elements(l) < 15
  print,'>>> --- 275.10-275.70 の全輝線（強い順）---'
  for i=0,nn-1 do print,'>>>   '+string(l[s[i]].wvl,format='(f9.4)')+'  '+ $
      string(l[s[i]].snote,format='(a-20)')+'  I_rel='+ $
      string(l[s[i]].int/l[s[0]].int,format='(f9.5)')
endif else print,'>>> ch_synthetic error: '+err
print,'>>> DONE'
exit
end
