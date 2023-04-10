; Proposal walen test
; V assumed to be l component only

function symmetric_walen, v, density, b, Tperp, Tpara, btime, ptime, t1, t2, t3, backward=backward

; We interpolate based on velocity

; interpolate bl to plasma time
ev2k = 1/(8.621738e-5)
mu0 = 4*!pi*1e-7
kb = 1.38064852e-23
mi = 1.6726219e-27

mult1 = -1
mult2 = 1

if keyword_set(backward) then begin
  mult1 = 1
  mult2 = -1
endif


bli = interpol(b[*,0], btime, ptime, /nan)
bmi = interpol(b[*,1], btime, ptime, /nan)
bni = interpol(b[*,2], btime, ptime, /nan)
bt2 = bli^2 + bmi^2 + bni^2

loc1 = where(ptime ge time_double(t1) and ptime le time_double(t2), count1)
loc2 = where(ptime ge time_double(t2) and ptime le time_double(t3), count2)

;-----------------------------------------------------
; First walen test
;-----------------------------------------------------

; Initial conditions
vl0 = 1e3*v[loc1[0]]
bl0 = 1e-9*bli[loc1[0]]
Tperp0 = ev2k*Tperp[loc1[0]]
Tpara0 = ev2k*Tpara[loc1[0]]
n0 = 1e6*density[loc1[0]]
rho0 = n0*mi
alpha0 = n0*kb*(Tpara0-Tperp0)*mu0/(bt2[loc1[0]])

; Now we do the time series
vl = 1e3*v[loc1]
bl = 1e-9*bli[loc1]
Tperpw = ev2k*Tperp[loc1]
Tparaw = ev2k*Tpara[loc1]
n = 1e6*density[loc1]
alpha = n*kb*(Tparaw-Tperpw)*mu0/(bt2[loc1])

; Now evaluate the walen relation
numerator = bl*(1-alpha) - bl0*(1-alpha0)
denominator = sqrt(rho0*mu0*(1-alpha0))
dVal1 = numerator/denominator
Vwl1 = 1e-3*(vl0 + mult1*dval1)

;-----------------------------------------------------
; Second walen test
;-----------------------------------------------------

; Initial conditions
vl0 = 1e3*v[loc2[count2-1]]
bl0 = 1e-9*bli[loc2[count2-1]]
Tperp0 = ev2k*Tperp[loc2[count2-1]]
Tpara0 = ev2k*Tpara[loc2[count2-1]]
n0 = 1e6*density[loc2[count2-1]]
rho0 = n0*mi
alpha0 = n0*kb*(Tpara0-Tperp0)*mu0/(bt2[loc2[count2-1]])

; Now we do the time series
vl = 1e3*v[loc2]
bl = 1e-9*bli[loc2]
Tperpw = ev2k*Tperpw[loc2]
Tparaw = ev2k*Tparaw[loc2]
n = 1e6*density[loc2]
alpha = n*kb*(Tpara-Tperp)*mu0/(bt2[loc2])

; Now evaluate the walen relation
numerator = bl*(1-alpha) - bl0*(1-alpha0)
denominator = sqrt(rho0*mu0*(1-alpha0))
dVal2 = numerator/denominator
Vwl2 = 1e-3*(vl0 + mult2*dval2)

outstr = {t1: ptime[loc1], $
            t2: ptime[loc2], $
            vwl1: vwl1, $
            vwl2: vwl2}

return, outstr


end