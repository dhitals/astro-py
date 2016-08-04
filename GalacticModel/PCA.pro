PRO PCA, obj, nstepsMC=nstepsMC
;; PRO PCA, ra, dec, theta, dist1, dist2, distErr1, distEr2, $
;;          pmRA, pmDEC, pmRAerr, pmDECerr, outfile=outfile
;
; This program reads in (or an input can be used) a list of binaries
; and calculates the number of stars expected in the LOS and volume
; occupied by the binary system by creating "fake" galaxies and
; randomly distributing the stars in that galaxy. It uses the
; rejection method (Press et al. 1992) to distribute the stars using
; the stellar density profile as calculated by Bochanski et al. (2009)
; using SDSS low-mass stars.
;
; This version is customized for calculating the P_CA (or P_f) of wide
; binaries with Adam Burgasser.
;
; Written by : Saurav Dhital
; Written on : April 17, 2009

astrolib

COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed
Rsun = 8500.d & Tsun = 0.d & Zsun = 15.d ; R, theta, & Z coords for the sun in parsecs
rho0 = 0.0064d                  ; rho0 from Juric et al. (2008)
au2pc = 1.d/206264.806          ; conversion from AU to parsecs
cell_size = 0.5d                ; size of one cell in degrees
max_dist = 2500.d               ; maximum allowed distance for simulated star
seed = ranseed()
                                ; no. Monte Carlo steps
IF KEYWORD_SET(nstepsMC) THEN nstepsMC=nstepsMC else nstepsMC = 1e6
print,' No. of MC steps:',nstepsMC

; **************************************************************************
; **************************************************************************
tstart = systime(0)
print,'Start Time: ',tstart
t_start = systime(1)
PM = 'yes' & RV = 'no'
; storage arrays
count = lonarr(5)

CASE obj of 
   1: BEGIN
      object = 'HIP67594' & outfile=object+'.out'
      COMMENT = 'Parameters for HIP67594'

      ra1 = 15.*ten(13,51,2.889) & dec1 = ten(23,45,45.29)    
      dist1  = 43.63   & sig_dist1 = 4     
      pmra1  = 48.51   & pmRAerr1  = 1.74
      pmDec1 = -76.71  & pmDecerr1 = 1.12

      ra2 = 15.*ten(13,49,26.4) & dec2 = ten(23,40,45.90)      
      dist2 = 36.0   & sig_dist2 = 24
      pmra2  =  105  & pmRAErr2  = 21
      pmDec2 = -70   & pmDecErr2 = 21

      theta = angdist(ra1,dec1,ra2,dec2)
   END
ENDCASE

ra0    = ra2                 ; system properties are subscripted with 0
dec0   = dec2 
theta0 = theta
;dist0  = 0.5*(dist1+dist2)
dist0 = dist2
temp = conv_to_galactic(ra0,dec0,dist0,R0,T0,Z0)

IF (PM EQ 'yes') THEN BEGIN
    pm0  = 0.5*[pmra1+pmra2, pmdec1+pmdec2]
    sigpm0  = [sqrt(pmRAerr1^2  + pmRAerr2^2),$
               sqrt(pmDECerr1^2 + pmDECerr2^2)]
ENDIF

sig_ddist0 = sqrt(sig_dist1^2+sig_dist2^2)
sig_dpm0 = sqrt(sigPM0[0]^2+sigPM0[1]^2)

; count the number of stars in each cell of length cell_size
n= n_elements(ra0)
nstars = round(count_nstars(ra0,dec0))
print,nstars
volume_src = 1/3.*!dpi * (0.5*theta0*au2pc)^2 * ((dist0+sig_ddist0)^3-((dist0-sig_ddist0)^3))
print,volume_src

FOR niter = 0l,nstepsMC-1 DO BEGIN
    temp = gen_nstars(ra0,dec0,nstars,ra,dec,dist)

    theta = angdist(ra0,dec0,ra,dec)
    ddist = abs(dist-dist0)           

    ; window,1,xsize=600,ysize=600
    ; plot,ra,dec,psym=3
    ; x = where(ddist le sig_ddist0 and ddist le 100) 
    ; oplot,ra[x],dec[x],psym=8,color=cgcolor('red')
    ; tvcircle,theta0/3600.,ra0,dec0,/data

    ; ************** COUNT FOR MATCHES  *******************
    ind1 = where(theta GE 0. AND theta LE theta0,c1) ; counts all stars within given theta and all d

    ind2 = where(theta GE 0. AND theta LE theta0 AND $
                 ddist LE sig_ddist0 AND ddist LE 100,c2) ; counts stars within given theta and 1-sigma of delta_distance

    IF (PM EQ 'yes' AND c2 GT 0) THEN BEGIN
        vel = gen_pm(R0,T0,Z0,ra0,dec0,dist0,c2) ; returns [[pmra], [pmdec],[rv]]
        
        dpmRA  = abs(vel[*,0]-pm0[0])/sigpm0[0]
        dpmDEC = abs(vel[*,1]-pm0[1])/sigpm0[1]
        ind3 = where(sqrt(dpmRA^2+dpmDEC^2) le 2,c3)

        IF (RV EQ 'yes') THEN BEGIN
            dRV = abs(vel[*,2]-RV0)/sigRV0
            ind4 = where(dpmRA^2+dpmDEC^2 le 2 and dRV le 1,c4)
;; FIX            ind5 = where((dpmRA^2+dpmDEC^2+dRV^2) le 3,c5)
         ENDIF ELSE BEGIN
            c4 = 0
            c5 = 0
         ENDELSE
    ENDIF ELSE BEGIN
       c3 = 0
       c4 = 0
       c5 = 0
    ENDELSE
    
    count += [c1,c2,c3,c4,c5] * 1LL

    IF (niter MOD (0.1*nstepsMC)) EQ 0 THEN $
      print,systime(),(100.*niter/nstepsMC),'% done.';    N = ',i+1 

ENDFOR                          ; END of ONE MC STEP
print,count,format='(5(i6))'    
P_f = count/nstepsMC

IF keyword_set(outfile) THEN lun = 1 else lun = -1

IF lun eq 1 THEN openw,lun,outfile
    
printf,lun,' '

printf,lun,COMMENT
printf,lun,'No. of MC realizations : ',nstepsMC,format='(a32,g10.4)'
printf,lun,' '

printf,lun,'RA : ',ra0,format='(a32,f11.6)'
printf,lun,'DEC : ',dec0,format='(a32,f11.6)'
printf,lun,'Angular Separation (arcsec): ',theta0,format='(a32,f8.2)'
printf,lun,'Distance (pc): ',dist1,'+/-',sig_dist1,';',dist2,'+/-',sig_dist2,format='(a32,f8.2,a4,f5.2,a1,f8.2,a4,f5.2)'
printf,lun,'Distance matched within (pc): ',dist0,'+/-',sig_ddist0,format='(a32,f8.2,a4,f5.2)'
printf,lun,'Volume of Search (pc^3): ',volume_src,format='(a32,f10.4)'
printf,lun,'No. of stars in 30x30 sq arcmin and 0-2500 pc volume: ',nstars,format='(a55,i5)'

;printf,lun,'Sigma DeltaDistance (pc): ',sig_ddist0,format='(a32,f6.2)'
printf,lun,''
IF (PM eq 'yes') THEN printf,lun,'Proper Motion in RA (mas/yr): ',pmra1,'+/-',pmRAerr1[0],';',pmRA2,'+/-',pmRAerr2,format='(a32,f8.2,a4,f5.2,a1,f8.2,a4,f5.2)'
IF (PM eq 'yes') THEN printf,lun,'Proper Motion in DEC (mas/yr): ',pmdec1,'+/-',pmDECerr1[0],';',pmDEC2,'+/-',pmDECerr2,format='(a32,f8.2,a4,f5.2,a1,f8.2,a4,f5.2)'
IF (RV eq 'yes') THEN printf,lun,'Radial Velocity (km/s): ',RV0,'+/-',sigRV0,format='(a32,f8.2,a4,f5.2)'
;IF (PM eq 'yes') THEN printf,lun,'pmRA matched within (mas/yr):',pm0[0],'+/-',
printf,lun,' '

printf,lun,'P_f(s) : ',P_f[0],format='(a32,f12.6)'
printf,lun,'P_f(s,d) : ',P_f[1],format='(a32,f12.6)'
IF (PM eq 'yes') THEN printf,lun,'P_f(s,d,mu) : ',P_f[2],format='(a32,f12.6)'
IF (RV eq 'yes') THEN printf,lun,'P_f(s,d,mu,RV) : ',P_f[3],format='(a32,f12.6)'
IF (RV eq 'yes') THEN printf,lun,'P_f(s,d,mu,RV) : ',P_f[4],format='(a32,f12.6)'
printf,lun,' '

printf,lun,'Start Time : ',tstart,format='(a32,a25)'
printf,lun,'TOTAL TIME TAKEN : ',(systime(1)-t_start)/3660.,'hours',format='(a32,f7.3,a10)'
printf,lun,'TIME PER ITERATION : ',(systime(1)-t_start)/(n*nstepsMC),'seconds',format='(a32,f6.3,a10)'
printf,lun,'END TIME : ',systime(0),format='(a32,a25)'
printf,lun,' '

IF lun eq 1 THEN close,lun

END
