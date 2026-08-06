; Ca XVII 192.858 のブレンドを解くために、192.70-193.05 A に何が居るかを
; CHIANTI 9.0.1 から出す。Ko et al. (2009) の方法を自分で組むための下ごしらえ。
;
;   bash scripts/idl/run_sswidl.sh scripts/idl/05_ca17_lines.pro logs/05.log
;   （SSW_INSTR に chianti が要る。run_sswidl.sh の既定に入っている環境なら不要）
;
; Warren+2012 p.6 の記述:
;   "The Ca XVII 192.858 A line is blended with Fe XI 192.813 A and a complex
;    of O V lines. We use the method outlined by Ko et al. (2009) to disentangle
;    this blend. ... The width of Ca XVII 192.858 A is limited to be within
;    0.05 mA of the width of Ca XIV 193.874 A."

set_plot,'z'

wmin = 192.60d & wmax = 193.10d
logt = [5.2, 5.4, 6.1, 6.7]          ; O V / O V / Fe XI / Ca XVII の形成温度あたり
dens = 1d9                            ; 活動領域コアの典型密度

print,'>>> CHIANTI !xuvtop = '+!xuvtop
print,'>>> window '+string(wmin,format='(f7.2)')+' - '+string(wmax,format='(f7.2)')+' A'
print,'>>> log Ne = '+string(alog10(dens),format='(f4.1)')

for it=0,n_elements(logt)-1 do begin

  ch_synthetic, wmin, wmax, output=out, err_msg=err, density=dens, $
                logt_isothermal=logt[it], logem_isothermal=27.0, $
                ioneq_name=!ioneq_file, /photons, /all

  if err ne '' then begin
    print,'>>> ERROR at logT='+string(logt[it],format='(f4.1)')+': '+err
    continue
  endif

  l = out.lines
  s = reverse(sort(l.int))
  n = n_elements(l) < 12
  print,'>>> --- logT = '+string(logt[it],format='(f4.2)')+' -----------------------------'
  for i=0,n-1 do begin
    j = s[i]
    print,'>>>   '+string(l[j].wvl,format='(f9.4)')+'  '+ $
          string(l[j].snote,format='(a-22)')+ $
          '  I_rel='+string(l[j].int/l[s[0]].int,format='(f8.5)')
  endfor

endfor

print,'>>> DONE'
exit
end
