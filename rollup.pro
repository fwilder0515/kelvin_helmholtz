; Do scatter plots to determine how rolled up a KH interval is.

function rollup, velnames, densnames, tempperp, temppara, lmnname, trange, flipm=flipm

get_data, velnames, data = vstr

get_data, densnames, data = denstr

get_data, tempperp, data = tperpstr

get_data, temppara, data = tparastr

get_data, lmnname, data = lmnstr

tloc = where(vstr.x gt time_double(trange[0]) and vstr.x lt time_double(trange[1]), counttime)

if counttime eq 0 then begin
  print, 'Time is out of range!!!'
  return, -1
endif

tempavg = (tparastr.y(tloc) + 2*tperpstr.y(tloc))/3.

vx = vstr.y(tloc,0)
vy = vstr.y(tloc,1)
vz = vstr.y(tloc,2)

mx = interpol(lmnstr.y(*, 0, 1), lmnstr.x, vstr.x)
my = interpol(lmnstr.y(*, 1, 1), lmnstr.x, vstr.x)
mz = interpol(lmnstr.y(*, 2, 1), lmnstr.x, vstr.x)

vm = mx*vx + my*vy + mz*vz

if keyword_set(flipm) then vm *= -1

n = denstr.y(tloc)

loc_low = where(tempavg lt 250)
loc_mid = where(tempavg ge 250 and tempavg le 900)
loc_high = where(tempavg gt 900)

outstruct = {high: loc_high, $
             mid: loc_mid, $
             low: loc_low, $
             vm: vm, $
             n: n, $
             ti: tempavg}

return, outstruct

end