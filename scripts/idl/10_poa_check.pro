; 自前で落とした PINTofALE (2015 版) が IDL 9.2 でコンパイルできるか確かめる。
;
;   PINTofALE=/scr/a000/c0234hotta/PINTofALE \
;     bash scripts/idl/run_sswidl.sh scripts/idl/10_poa_check.pro logs/10.log
;
; SSW 同梱の /opt/ssw/packages/poa は fitting/ 以下のパーミッションが
; drwxr--r-- で入れないうえ 2004 年版なので使わない。

set_plot,'z'

poa = getenv('PINTofALE')
print,'>>> PINTofALE = '+poa

;; PoA のパスを !path の先頭に足す（SSW 版より前に来るように）
!path = expand_path('+'+poa+'/pro') + path_sep(/search_path) + !path

r = ['mcmc_dem','mcmc_abund','likeli','findscale','varsmooth','mixie','pred_flx', $
     'lineflx','mk_dem','getabund','rdabund','rd_line','rd_list','rd_ioneq', $
     'fold_ioneq','cat_ln','ionabs','ismtau','bamabs','getpoadef','hastogram', $
     'rebinw','rebinx','roofn','syze','findex','is_keyword_set','wvlt_scale', $
     'mcmc_dem_whiskerplot','initale_pro']
nbad = 0
for i=0,n_elements(r)-1 do begin
  f = file_which(!path, r[i]+'.pro')
  if f eq '' then begin
    print,'>>> MISSING  '+r[i] & nbad = nbad+1 & continue
  endif
  ;; 実際にコンパイルしてみる（構文が IDL 9.2 で通るか）
  catch, err
  if err ne 0 then begin
    catch,/cancel
    print,'>>> FAIL     '+r[i]+'  ('+!error_state.msg+')'
    nbad = nbad+1
    continue
  endif
  resolve_routine, r[i], /either, /compile_full_file, /no_recompile
  catch,/cancel
  print,'>>> ok       '+r[i]+'  <- '+f
endfor

print,'>>> failures = '+strtrim(nbad,2)
print,'>>> DONE'
exit
end
