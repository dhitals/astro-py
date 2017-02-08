FUNCTION dist_sdm,psfg,psfr,psfi,psfz

;; statistical parallax relations from Bochanski et al. (2012)

absmag = fltarr(n_elements(psfr)) & dist = absmag

gr = psfg - psfr
gi = psfg - psfi
iz = psfi - psfz
rz = psfr - psfz

;; Bochanski et al. (2012)
;; 6th-order polynomial for median r-z vs. g-r
;; valid for 0.8 < r-z < 3.6 (K5-L1)
c = [-14.8952, +47.8830, -56.3038, +34.1232, -11.2926, 1.94407, -0.136479]
gr_dM = c[0] + c[1]*rz + c[2]*rz^2 + c[3]*rz^3+ c[4]*rz^4 + c[5]*rz^5 + c[6]*rz^6
delta_gr = gr - gr_dM

;; polynomial for absmag vs. r-z, excess g-r
;; valid for 1.0 < r-z < 2.0 (M0-M4) ; 0.0 < delta_gr < 0.5
;; sigma_Mr = 0.41
a = [7.9547,1.8102,-0.17347,7.7038,-1.4170]

ind0 = where(rz ge 0.8 OR rz le 3.60,ct0)
if ct0 ge 1 then begin
   absmag[ind0] = a[0] + a[1]*rz[ind0]  + a[2]*rz[ind0]^2 + $
                  a[3] * delta_gr[ind0] + a[4]*rz[ind0]*delta_gr[ind0]
   dist[ind0] = 10^((psfr[ind0] - absmag[ind0])/5. + 1)
endif
   
;; return scalar is scalar was supplied
dist = (n_elements(psfr) eq 1) ?  dist[0]:dist

return,dist
END
