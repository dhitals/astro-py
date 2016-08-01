function dist_ms,psfg,psfr,psfi,psfz
  
;; This function calculates the absolute magnitude and the distance to 
;; main sequence dwarfs based on their colors: g-i (F0-K5), r-z (K5-M8), and i-z (>M8)

;SpTy=[0.60,0.49,0.44,0.39,0.36,0.20,0.15,0.12,0.11,0.08] ;mass of M0-M9

absmag = fltarr(n_elements(psfr)) & dist = absmag
gi = psfg - psfi
rz = psfr - psfz
iz = psfi - psfz

;don't calculate the distance for points outside the range
ind0 = where(rz lt -0.01 OR  iz gt 3.20,ct0)
ind1 = where(rz ge -0.01 and rz lt 0.50,ct1)
ind2 = where(rz ge  0.50 and rz le 4.53,ct2)
ind3 = where(iz ge  1.70 and iz le 3.20,ct3)

;; Covey et al. 2007: F0 - <K5
if ct1 ge 1 then begin
   absmag[ind1] = 2.845 + 1.656*gi[ind1] + 3.863*gi[ind1]^2 - 1.795*gi[ind1]^3
   dist[ind1] = 10^((psfg[ind1] - absmag[ind1])/5. + 1)
endif

;; Bochanski et al. (2010) - about K5 - M9
if ct2 ge 1 then begin
   absmag[ind2] = 5.190 + 2.474*rz[ind2] + 0.434*rz[ind2]^2 - 0.086*rz[ind2]^3
   dist[ind2] = 10^((psfr[ind2] - absmag[ind2])/5. + 1)
endif

;; Schmidt et al. (2010) - >M9
if ct3 ge 1 then begin
   absmag[ind3] = -23.27 + 38.4*iz[ind3] - 11.11*iz[ind3]^2 + 1.064*iz[ind3]^3
   dist[ind3] = 10^(-1.*(psfi[ind3] - absmag[ind3])/5. + 1)
endif

;; don't calculate the distance for points outside the range
if ct0 ge 1 then dist[ind0] = -1

;; return scalar is scalar was supplied
dist = (n_elements(psfr) eq 1) ?  dist[0]:dist

return,dist
end
