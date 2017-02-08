FUNCTION sptype_dm,rz,nosubtype=nosubtype
;; returns the spectral type of M0-M9 based on the r-z color
;; if NOSUBTYPE specified, then return ONLY the type

; all possible SpTypes - all types have 10 subtypes except for the
;                        undefined K6, K8, and K9
spName = ['O0', 'O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', $
          'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', $
          'A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', $
          'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', $
          'G0', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', $
          'K0', 'K1', 'K2', 'K3', 'K4', 'K5', 'K7' , $
          'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6','M7', 'M8', 'M9','L']

n = n_elements(rz)
spType = strarr(n)
spNum = fltarr(n)

;; B8-K5 -- Covey et al. (2007)
ind1 = where(rz lt 0.6,c1)
if c1 gt 0 then $
  spNum[ind1] = 33.5674 + 47.8054*rz[ind1] + 20.0984*rz[ind1]^2 - 33.6663*rz[ind1]^3 - 58.2037*rz[ind1]^4

; K5 or later -- Bochanski et al. (2010)
ind2 = where(rz ge 0.6 and rz le 4.5,c2)
if c2 gt 0 then $
  spNum[ind2] = 50.3763 + 9.03538*rz[ind2] - 2.97758*rz[ind2]^2 + 0.516391*rz[ind2]^3 - 0.0282781*rz[ind2]^4 

; L0 or later -- just call them L
ind3 = where(rz gt  4.5,c3)
if c3 gt 0 then spType[ind3] = 'L'

type = spName[floor(spNum)]
subtype = spNum - floor(spNum)

; if subtype > 0.95, increase type
ind = where(subtype gt 0.95)
if ind[0] ne -1 then begin
    for j=0,n_elements(ind)-1 do begin
        type[ind[j]] = spName[floor(spNum[ind[j]])+1]
        subtype[ind[j]] = 0
    endfor
endif

IF NOT keyword_set(nosubtype) THEN $
   spType = type + strmid(string(subtype,format='(f3.1)'),1,2) $
ELSE spType = type

return, spType

END
 
