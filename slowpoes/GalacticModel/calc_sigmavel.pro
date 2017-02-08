FUNCTION CALC_SIGMAVEL,Z

; U_thin,V_thin,W_thin,U_thick,V_thick,W_thick
; Values obtained by fitting sigma = coeff * Z^power 
; data from Bochanski et al. (2006)
; see ~/sdss/uw/velocity_ellipsoid.pro[.ps] FOR fitting algorithm[fit]
coeff = [7.085,3.199,3.702,10.383,1.105,5.403]
power = [0.276,0.354,0.307, 0.285,0.625,0.309]

; calculate sigma_vel from the empirical power-law fit
sigmaa = coeff * abs(Z)^power

RETURN,sigmaa
END


