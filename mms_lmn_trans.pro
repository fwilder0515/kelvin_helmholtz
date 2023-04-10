
pro mms_lmn_trans, str_in, str_out, sc=sc, t1=t1, t2=t2, savelmn=savelmn, use=use

;keyword savelmn stores congrid version of Nvec and GSE->LMN matrix
;keyword use will use previously calculated Nvec stored as nvec_gse_'sc'

if not keyword_set(sc) then sc='mms1'

if keyword_set(use) then goto,skip


;get Fairfield normal vector
re=6378.16		;Earth equatorial radius [km]
get_data, sc+'_mec_r_gse',data=p	;s/c GSE position in Re

nn=n_elements(p.x)
nvec=fltarr(nn,3)
for i=0,nn-1 do begin
  data=p.y(i,*)
  mpnormal,data,vec
  nvec(i,*)=vec
endfor
store_data,'nvec_gse_'+sc,data={x:p.x, y:nvec}, $
           dlimit={spec:0,ystyle:1,ytitle:'n fair'}

print,'*** Got Fairfield normal in gse ***'


skip:
print,' input quantity: ',str_in
print,'output quantity: ',str_out


;Align Ngse with Data interval
get_data,'nvec_gse_'+sc,data=n
get_data,str_in,data=d

if not keyword_set(t1) then tt1=d.x(0) else tt1=time_double(t1)
if not keyword_set(t2) then tt2=d.x(n_elements(d.x)-1) else tt2=time_double(t2)


ind1=where((d.x ge tt1) and (d.x le tt2),c1)
ind2=where((n.x ge tt1) and (n.x le tt2),c2)
t_new=d.x(ind1)
d_new=d.y(ind1,*)
n_new=congrid(n.y(ind2,*),c1,3)

if keyword_set(savelmn) then begin
  store_data,'nvec_gse_'+sc+'_congrid',data={x:t_new,y:n_new}
  options,'nvec_gse_'+sc+'_congrid','ytitle','n fair congrid'
endif

;help,n_new,/str
;help,d_new,/str
;print,c1,c2

vec_out=fltarr(c1,3)
lmn_matrix=fltarr(c1,3,3)
i=long(0)
for i=0l, long(c1-1) do begin
  lmnmat, n_new(i,*), gselmn
  lmn_matrix(i,*,*)=gselmn
  vec_out(i,*)= transpose(gselmn)#reform(d_new(i,*))
endfor

;phi= atan(vec_out(*,1),vec_out(*,0))/!dtor
;lamda= atan(vec_out(*,2),sqrt(vec_out(*,0)^2+vec_out(*,1)^2))/!dtor
;vec_mag= sqrt(vec_out(*,0)^2+vec_out(*,1)^2+vec_out(*,2)^2)

store_data,str_out,data={x:t_new,y:vec_out}
if keyword_set(savelmn) then store_data,sc+'_gselmn',data={x:t_new,y:lmn_matrix}

;if keyword_set(str_mag) then store_data,str_mag,data={ytitle:str_mag,x:t_new,y:vec_mag}
;if keyword_set(str_phi) then store_data,str_phi,data={ytitle:str_phi,x:t_new,y:phi},dlim=dlim
;if keyword_set(str_lamda) then store_data,str_lamda,data={ytitle:str_lamda,xtitle:'Time',x:d.x,y:lamda},dlim=dlim

end
