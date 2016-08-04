PRO get_avgPM, ra, dec,dist

;; This program returns the average 3D velocity (Pm, RV) for a given
;; ra, dec, and distance.
;; Adapted from the SLoWPoKES Galactic model

; Written by : Saurav Dhital
; Written on : April 17, 2009

FORWARD_FUNCTION gen_vel,gen_2Dgaussian,calc_sigmaVel
FORWARD_FUNCTION conv_to_galactic,conv_to_equatorial,calc_UVW
astrolib

COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed
Rsun = 8500.d & Tsun = 0.d & Zsun = 15.d ; R, theta, & Z coords for the sun in parsecs
rho0 = 0.0064d                  ; rho0 from Juric et al. (2008)
au2pc = 1.d/206264.806          ; conversion from AU to parsecs
cell_size = 0.5d                ; size of one cell in degrees
max_dist = 2500.d               ; maximum allowed distance for simulated star
seed = ranseed()

; **************************************************************************
; **************************************************************************
tstart = systime(0)
print,'Start Time: ',tstart
t_start = systime(1)
PM = 'no' & RV = 'no'
; storage arrays
count = lonarr(5)

IF N_PARAMS() EQ 0 THEN BEGIN
    COMMENT = "Parameters for KELT-2"
    ra = 15.*ten(6,10,39.35)
    dec = ten(30,57,25.7)
    dist = 128.9  & sig_dist = 6.9
    theta = 2.29 
    ;ddist = abs(dist1-dist2)
ENDIF

num = 1e6 ;; draw 1000 stars

;; change equatorial coords to galactic
temp = conv_to_galactic(ra,dec,dist,R,T,Z)

vel = gen_pm(R,T,Z,ra,dec,dist,num)

; pmRA
plothist,vel[*,0],x1,y1,bin=1,xtitle='pmRA (mas/yr)'
result = fit_gaussian(x1,y1,vals=[0,-10,20,1],FIT=fit,fixed=[1,0,0,0],perror=perror,/Quiet) 
oplot,x1,fit,color=cgcolor('red')
print,'pmRA (mas/yr): ',result[1],result[2],format='(a20,f8.2,f8.2)'

wait,10
; pmDEC
plothist,vel[*,1],x1,y1,bin=1,xtitle='RV (km/s)'
result = fit_gaussian(x1,y1,vals=[0,-10,20,1],FIT=fit,fixed=[1,0,0,0],perror=perror,/Quiet) 
oplot,x1,fit,color=cgcolor('red')
print,'pmDEC (mas/yr): ',result[1],result[2],format='(a20,f8.2,f8.2)'

wait,10
; RV
plothist,vel[*,2],x1,y1,bin=1,xtitle='RV (km/s)'
result = fit_gaussian(x1,y1,vals=[0,-10,20,1],FIT=fit,fixed=[1,0,0,0],perror=perror,/Quiet) 
oplot,x1,fit,color=cgcolor('red')
print,'RA (km/s)): ',result[1],result[2],format='(a20,f8.2,f8.2)'

stop
END
