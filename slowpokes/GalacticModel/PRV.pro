PRO PRV, ra, dec, theta, dist, distErr, pmRA, pmDEC, pmRAerr, pmDECerr, dRV=dRV

;; This program reads in (or an input can be used) a list of binaries
;; and calculates the number of stars expected in the LOS and volume
;; occupied by the binary system by creating "fake" galaxies and
;; randomly distributing the stars in that galaxy. It uses the
;; rejection method (Press et al. 1992) to distribute the stars using
;; the stellar density profile as calculated by Bochanski et al. (2009)
;; using SDSS low-mass stars.
;;
;; This version is customized for to measure the intrinsic scatter in
;; delta_RV for CPM pairs generated using the Galactic model from
;; Dhital et al. (2010)

;;
;; Written by : Saurav Dhital
;; Written on : January 25, 2011

FORWARD_FUNCTION gen_nstars,gen_vel,gen_2Dgaussian
FORWARD_FUNCTION calc_sigmaVel,calc_rho,count_nstars
FORWARD_FUNCTION conv_to_galactic,conv_to_equatorial,calc_UVW
astrolib

COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed
Rsun = 8500.d & Tsun = 0.d & Zsun = 15.d ; R, theta, & Z coords for the sun in parsecs
rho0 = 0.0064d                  ; rho0 from Juric et al. (2008)
au2pc = 1.d/206264.806          ; conversion from AU to parsecs
cell_size = 0.5d                ; size of one cell in degrees
max_dist = 2500.d               ; maximum allowed distance for simulated star
seed = ranseed()
                                ; no. Monte Carlo steps

; **************************************************************************
; **************************************************************************
tstart = systime(0)
print,'Start Time: ',tstart
t_start = systime(1)

ra0    = ra                 ; system properties are subscripted with 0
dec0   = dec 
theta0 = theta
dist0  = 0.5*(dist[0]+dist[1])
sig_ddist0 = sqrt(disterr[0]^2+disterr[1]^2)
temp = conv_to_galactic(ra0,dec0,dist0,R0,T0,Z0)

pm0  = 0.5*[pmra[0]+pmra[1], pmdec[0]+pmdec[1]]
sigpm0  = [sqrt(pmRAerr[0]^2  + pmRAerr[0]^2),$
           sqrt(pmDECerr[1]^2 + pmDECerr[1]^2)]

; count the number of stars in each cell of length cell_size
n= n_elements(ra0)
nstars = round(count_nstars(ra0,dec0))

temp = gen_nstars(ra0,dec0,nstars,ra,dec,dist)

theta = angdist(ra0,dec0,ra,dec)
ddist = abs(dist-dist0)           

;; ************** COUNT FOR MATCHES  *******************
;; counts all stars within given theta and all d
ind1 = where(theta GE 0. AND theta LE theta0,c1) 

;; counts stars within given theta and 1-sigma of delta_distance
ind2 = where(theta GE 0. AND theta LE theta0 AND $
             ddist LE sig_ddist0 AND ddist LE 100,c2) 

;; since generating velocity in computationally intensive,
;; we'll only generate velocities for stars that are
;; already matched in position
vel = gen_vel(R0,T0,Z0,ra0,dec0,dist0,c2) ; returns [[pmra], [pmdec],[rv]]

dpmRA  = abs(vel[*,0]-pm0[0])/sigpm0[0]
dpmDEC = abs(vel[*,1]-pm0[1])/sigpm0[1]
ind3 = where(dpmRA^2+dpmDEC^2 le 2,c3)

;; now, that we have matched 3-d position & PM, we can see what
;; the distribution of dRV looks like. But (for each LOS) we need
;; only 2 stars

;; easy fix if N > 2
rv = (c3 gt 2) ? vel[ind3[0:1],2]:vel[ind3,2]
;    IF c3 gt 2 THEN rv = vel[ind3[0:1],2] ELSE rv = vel[ind3,2]

;; regenerate stars while N < 2
WHILE c3 lt 2 DO BEGIN
   temp = gen_nstars(ra0,dec0,nstars,ra,dec,dist)
   theta = angdist(ra0,dec0,ra,dec)
   ddist = abs(dist-dist0)           
   
   ind2 = where(theta GE 0. AND theta LE theta0 AND $
                ddist LE sig_ddist0 AND ddist LE 100,c2) 
   vel = gen_vel(R0,T0,Z0,ra0,dec0,dist0,c2) ; returns [[pmra], [pmdec],[rv]]
   
   dpmRA  = abs(vel[*,0]-pm0[0])/sigpm0[0]
   dpmDEC = abs(vel[*,1]-pm0[1])/sigpm0[1]
   ind3 = where(dpmRA^2+dpmDEC^2 le 2,c3)
   
   vel = gen_vel(R0,T0,Z0,ra0,dec0,dist0,c2) ; returns [[pmra], [pmdec],[rv]]
   dpmRA  = abs(vel[*,0]-pm0[0])/sigpm0[0]
   dpmDEC = abs(vel[*,1]-pm0[1])/sigpm0[1]
   ind4 = where(dpmRA^2+dpmDEC^2 le 2,c4)
   
   ;; if CPM pair, get the RV
   IF c4 gt 0 THEN BEGIN
      IF (c3 eq 0 and c4 eq 1) THEN rv = vel[ind4,2]
      IF (c3 eq 0 and c4 ge 2) THEN rv = vel[ind4[0:1],2]
      
      IF (c3 eq 1 and c4 ge 1) THEN rv = [rv,vel[ind4[0],2]]
      
      c3 =n_elements(rv)
   ENDIF
ENDWHILE

dRV = rv[0]-rv[1]

END


; **********************************************************
; ***************      FUNCTIONS             ***************
; **********************************************************

FUNCTION GEN_VEL,R0,T0,Z0,ra0,dec0,dist0,num
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

sigmaa = calc_sigmaVel(Z0)      ; calculate the UVW velocity dispersions
                                ; returns [U_thin,V_thin,W_thin,U_thick,V_thick,W_thick]
temp = calc_rho(R0,Z0,frac=frac); calc the frac of thin/thick disk stars 
                                ; returns frac = [f_thin, f_thick, f_halo]
vel = calc_UVW(R0,T0,Z0)-calc_UVW(Rsun,Tsun,Zsun) ; convert to cartesian velocities
                                ; returns [U,V,W]

    ; draw from both the thin and thick disks for UVW velocities
U = gen_2Dgaussian(vel[0],sigmaa[0],sigmaa[3],frac[0],1-frac[0],num)
V = gen_2Dgaussian(vel[1],sigmaa[1],sigmaa[4],frac[0],1-frac[0],num)
W = gen_2Dgaussian(vel[2],sigmaa[2],sigmaa[5],frac[0],1-frac[0],num)

; change UVW to pmra and pmdec
ra = fltarr(num)+ra0
dec = fltarr(num)+dec0
dist = fltarr(num)+dist0
gal_uvw_pm,pmra,pmdec,rv,U=U,V=V,W=W,/LSR,ra=ra,dec=dec,dist=dist

RETURN,[[pmra],[pmdec],[rv]]
END

FUNCTION GEN_2DGAUSSIAN,mu,sig1,sig2,f1,f2,num

n_acc = 0l                ; number of stars accepted
WHILE n_acc LT num DO BEGIN
    n_lft = num - n_acc ; no. of needed stars
    
    X = randomu(seed,n_lft)*10*sig2 - 5*sig2 ; generate random numbers between (-5,5) sigma of a normalized gaussian

;    X = num_gen(-5*sig2,5*sig2,0.1)
    z1 = (X-mu)/sig1 & z2 = (X-mu)/sig2
    G1 = 1/(sqrt(2*!dpi)*sig1)*exp(-z1^2/2) ; thin disk
    G2 = 1/(sqrt(2*!dpi)*sig2)*exp(-z2^2/2) ; thin disk
 
    Px = f1*G1+f2*G2
    Fx = 1.2*max(Px)

    rand = randomu(seed,n_lft)
    ind = where(rand LT Px/Fx, count) ; uniformm deviate comp function

    IF count NE 0 THEN BEGIN
        IF n_acc EQ 0 THEN Xarr = x[ind] ELSE Xarr = [Xarr,x[ind]]
        n_acc += count
    ENDIF
ENDWHILE

RETURN,Xarr
END



FUNCTION COUNT_NSTARS,ra,dec
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

ddist = 5                       ; steps in distance in pc
n = max_dist/ddist + 1
dist = findgen(n)*ddist         ; 0 < d < 2500 in 5 pc steps
rho  = fLTarr(n) & nstars = fLTarr(n) ; create an array to store rho FOR each d
 
; define fractional positions so that rho can be averaGEd
x = [-0.5, 0.0, 0.5,-0.25, 0.25,-0.5,0.0,0.5,-0.25,0.25,-0.5,0.0,0.5]*cell_size
y = [-0.5,-0.5,-0.5,-0.25,-0.25, 0.0,0.0,0.0, 0.25,0.25, 0.5,0.5,0.5]*cell_size

; chanGE to galactic coordinates (R,Z) from input convert (ra, dec, dist)
FOR k = 0,n_elements(dist)-1 DO BEGIN
    temp = conv_to_galactic(ra+x,dec+y,dist[k],R,T,Z)
    
    rho[k] = mean(calc_rho(R,Z)) ; calculate the stellar density 
    vol = !dpi * (1800.*dist[k]*au2pc)^2 * 5.d ; volume at that d = pi * r^2 * h
    nstars[k] = rho[k] * vol
ENDFOR
nstars_tot = total(nstars)

RETURN,nstars_tot
END


FUNCTION GEN_NSTARS,ra0,dec0,num,ra,dec,dist
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

; ra0,dec0,num - input parameters
; ra,dec,dist - output arrays for the generated stars

seed = ranseed()
n_acc = 0l                ; number of stars accepted

WHILE n_acc LT num DO BEGIN
    n_lft = num - n_acc ; no. of needed stars
    
    ra1   = ra0 + randomu(seed,n_lft)*cell_size
    dec1  = dec0+ randomu(seed,n_lft)*cell_size
    dist1 = randomu(seed,n_lft)*max_dist

    temp = conv_to_galactic(ra1,dec1,dist1,R,T,Z)
    
    rho = calc_rho(R,Z)

    rand = randomu(seed,n_lft)

    ; accept if random number is less than rho(R,Z)/rho0
    ind = where(rand LT rho/rho0, count)

    IF count NE 0 THEN BEGIN
        IF n_acc EQ 0 THEN BEGIN
            ra   = ra1[ind]
            dec  = dec1[ind]
            dist = dist1[ind]
        ENDIF ELSE BEGIN
            ra   = [ra,ra1[ind]]
            dec  = [dec,dec1[ind]]
            dist = [dist,dist1[ind]]
        ENDELSE
        n_acc += count
    ENDIF
ENDWHILE

RETURN,0
END


FUNCTION CALC_RHO,R,Z,frac=frac
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

H_thin = 260.d  & H_thick = 900.d  ; scale height in pc
L_thin = 2500.d & L_thick = 3500.d ; scale length in pc

f_thick = 0.09d & f_halo = 0.0025  ; stellar density in the solar neighborhood
f_thin = 1 - f_thick - f_halo

r_halo = 2.77                   ; halo density gradient
q = 0.64                        ; halo flattening parameter

rho_thin  = rho0 * exp(-abs(Z-Zsun)/H_thin)  * exp(-(R-Rsun)/L_thin)
rho_thick = rho0 * exp(-abs(Z-Zsun)/H_thick) * exp(-(R-Rsun)/L_thick)
rho_halo  = rho0 * (Rsun/sqrt(R^2+(Z/q)^2))^r_halo

rho  = f_thin*rho_thin + f_thick*rho_thick + f_halo*rho_halo
frac = [f_thin*rho_thin, f_thick*rho_thick, f_halo*rho_halo]/rho

RETURN,rho
END


FUNCTION CALC_UVW,R,theta,Z

Rdot = 0.d 
Tdot = (226.-0.013*Z-1.56e-5*Z^2)/R ; will later  be converted to Tdot(Z) 
Zdot = 0.d                          ; typical values for the MW in km/s

theta = theta*!DTOR

Xdot =   Rdot*COS(theta) - Tdot*SIN(theta)
Ydot = -(Rdot*SIN(theta) + Tdot*COS(theta))

RETURN,[Xdot,Ydot,Zdot]

END

FUNCTION CONV_TO_EQUATORIAL,R,theta,Z,ra,dec,d
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

dcosb = SQRT(R^2+Rsun^2-2*R*Rsun*COS(theta*!DTOR))

l = ASIN(R*SIN(theta*!DTOR)/dcosb) * !RADEG
IF l LT 0 THEN BEGIN                          ; 3rd or 4th quadrant
   IF cos(l) GT 0 THEN l = 180-l ELSE l=360+l ; 3rd else 4th
ENDIF

IF l GE 0 THEN BEGIN                          ; 1st or 2nd quadrant
   IF cos(l) GT 0 THEN l=l ELSE l=180-l       ; 1st else 2nd
ENDIF

b = ABS(ATAN((Z-Zsun)/dcosb) + ATAN(Zsun/Rsun)) * !RADEG
IF z LT 0 THEN b = -b
d = dcosb/COS(b*!DTOR) 
euler,l,b,ra,dec,2

RETURN,0
END


FUNCTION CONV_TO_GALACTIC,ra,dec,d,r,t,z
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

euler,ra,dec,l,b,1

r = SQRT( (d*COS(b*!DTOR))^2 + Rsun *(Rsun- 2*d*COS(b*!DTOR)*COS(l*!DTOR)))
t = ASIN(d*SIN(l*!DTOR)*COS(b*!DTOR)/R) * !RADEG
z = Zsun + d * SIN(b*!DTOR - ATAN(Zsun/Rsun))

RETURN, 0
END


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



