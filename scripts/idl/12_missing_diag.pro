; Fe XII 195.119 の 195.1039 A のサンプルだけが SSW level-1 で 1/5 に落ちている。
; 箱内の画素ごとの値を見て、欠損値の扱いが原因かを確かめる。
set_plot,'z'
!quiet=1
l1 = getenv('EIS_L1')
wd = eis_getwindata(l1, 195.119, /refill, /quiet)
print,'>>> missing value tag = ',wd.missing
d = min(abs(wd.wvl-195.1039), j)
print,'>>> sample index ',j,'  wvl=',wd.wvl[j]
v = reform(wd.int[j,32:39,244:273])
print,'>>> n=',n_elements(v),' min=',min(v),' max=',max(v),' median=',median(v)
h = where(v lt 0, nneg)
print,'>>> negative samples: ',nneg
h2 = where(v lt 1000, nlow)
print,'>>> samples < 1000: ',nlow
print,'>>> sorted values (first 20):'
s = v[sort(v)]
print,s[0:19]
;; 隣のサンプルと比べる
v2 = reform(wd.int[j-1,32:39,244:273])
v3 = reform(wd.int[j+1,32:39,244:273])
print,'>>> neighbour medians: ',median(v2), median(v), median(v3)
;; 誤差配列も見る
e = reform(wd.err[j,32:39,244:273])
print,'>>> err  min/med/max = ',min(e),median(e),max(e)
print,'>>> DONE'
exit
end
