; Ca XVII 192.858 のブレンドを解くのに必要な原子データ比を CHIANTI 9.0.1 から出す。
; ch_synthetic を全イオンで回すと遅いので、必要なイオンだけ emiss_calc する。
;
;   bash scripts/idl/run_sswidl.sh scripts/idl/05b_ca17_ratios.pro logs/05b.log
;
; 欲しいもの:
;   (a) O V の 192.7-192.95 の多重線の相対強度（分岐比。密度・温度にほぼ依存しない）
;   (b) Fe XI 192.813 / Fe XI 188.216 の比（192.813 を 188.216 から見積もるため）
;   (c) Ca XVII 192.858 の位置確認
; Ko et al. (2009) はこの (a)(b) を使って 192.8 のブレンドを解く。

set_plot,'z'
!quiet = 1

logT = 4.0 + 0.05*findgen(81)          ; 4.0 - 8.0
temp = 10.d0^logT
dens = 1.d9

print,'>>> CHIANTI '+!xuvtop
print,'>>> log Ne = 9.0'

;; ---------------- O V (Z=8, ion=5) ----------------
em = emiss_calc(8, 5, temp=logT, dens=alog10(dens), /no_de, /quiet)
w  = em.lambda
g  = where(w gt 192.6 and w lt 193.1, ng)
print,'>>> O V lines in 192.6-193.1 : '+strtrim(ng,2)
it = (where(abs(logT-5.4) lt 0.001))[0]        ; O V の形成温度あたり
if ng gt 0 then begin
  ref = 0.d
  for i=0,ng-1 do if em[g[i]].em[it] gt ref then ref = em[g[i]].em[it]
  for i=0,ng-1 do print,'>>>   OV  '+string(w[g[i]],format='(f9.4)')+ $
      '  em(logT=5.4)='+string(em[g[i]].em[it],format='(e11.4)')+ $
      '  rel='+string(em[g[i]].em[it]/ref,format='(f8.5)')
endif

;; ---------------- Fe XI (Z=26, ion=11) ----------------
em2 = emiss_calc(26, 11, temp=logT, dens=alog10(dens), /no_de, /quiet)
w2  = em2.lambda
it2 = (where(abs(logT-6.10) lt 0.001))[0]
;; 192.813 と 188.216 を拾う
d1 = min(abs(w2-192.813), j1)
d2 = min(abs(w2-188.216), j2)
print,'>>> Fe XI 192.813 -> '+string(w2[j1],format='(f9.4)')+ $
      '  em='+string(em2[j1].em[it2],format='(e11.4)')
print,'>>> Fe XI 188.216 -> '+string(w2[j2],format='(f9.4)')+ $
      '  em='+string(em2[j2].em[it2],format='(e11.4)')
print,'>>> ratio 192.813/188.216 (logT=6.10, logNe=9) = '+ $
      string(em2[j1].em[it2]/em2[j2].em[it2],format='(f9.5)')
;; 温度依存も見ておく（比が温度でどれだけ動くか＝この方法の弱点）
for lt0=6.0,6.3,0.1 do begin
  k = (where(abs(logT-lt0) lt 0.001))[0]
  print,'>>>   at logT='+string(lt0,format='(f4.2)')+ $
        '  ratio='+string(em2[j1].em[k]/em2[j2].em[k],format='(f9.5)')
endfor

;; ---------------- Ca XVII (Z=20, ion=17) ----------------
em3 = emiss_calc(20, 17, temp=logT, dens=alog10(dens), /no_de, /quiet)
w3  = em3.lambda
g3  = where(w3 gt 192.5 and w3 lt 193.2, n3)
print,'>>> Ca XVII lines in 192.5-193.2 : '+strtrim(n3,2)
for i=0,n3-1 do print,'>>>   CaXVII '+string(w3[g3[i]],format='(f9.4)')

print,'>>> DONE'
exit
end
