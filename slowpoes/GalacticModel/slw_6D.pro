PRO SLW_6D,arg1,arg2,infile=infile,outfile=outfile,chunk=chunk,nstepsMC=nstepsMC

; This program reads in (or an input can be used) a list of binaries
; and calculates the number of stars expected in the LOS and volume
; occupied by the binary system by creating "fake" galaxies and
; randomly distributing the stars in that galaxy. It uses the
; rejection method (Press et al. 1992) to distribute the stars using
; the stellar density profile as calculated by Bochanski et al. (2009)
; using SDSS low-mass stars.

; Written by : Saurav Dhital
; Written on : April 17, 2009

astrolib

COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed
Rsun = 8500.d & Tsun = 0.d & Zsun = 15.d ; R, theta, & Z coords FOR the sun in parsecs
rho0 = 0.0064d                  ; rho0 from Juric et al. (2008)
au2pc = 1.d/206264.806          ; conversion from AU to parsecs
cell_size = 0.5d                ; size of one cell in degrees
max_dist = 2500.d               ; maximum allowed distance for simulated star
seed = ranseed()

IF KEYWORD_SET(nstepsMC) THEN nstepsMC=nstepsMC else nstepsMC = 5 ; no. Monte Carlo steps
print,' No. of MC steps:',nstepsMC

; **************************************************************************
; **************************************************************************
tstart = systime(0)
print,'Start Time: ',tstart
t_start = systime(1)
PM = 'yes' & RV = 'no'
; storage arrays

IF NOT keyword_set(INFILE) THEN infile = 'kolby/bry.v2.fits'
bry = mrdfits(infile,1)
n = n_elements(bry)
bry_theta = angdist(bry.ra[0],bry.dec[0],bry.ra[1],bry.dec[1])

; **************************************************************************
chunk_size = 1000L

IF keyword_set(chunk) THEN begin
    x0 = (chunk-1)*chunk_size & x1 = x0+chunk_size-1
    if chunk eq 103 then x1 = n-1
    print,chunk,x0,x1
    bry = bry[x0:x1] 
end else chunk = '00'

IF n_params() EQ 1 THEN bry = bry[round(randomu(seed,arg1)*n)] ; randomize if not doing the complete simulation
IF n_params() EQ 2 THEN bry = bry[arg1:arg2] ; do the model on a subset
n = n_elements(bry)
; **************************************************************************

print,''
print,'No. of candidate pairs: ',n
print,'No. of MC steps       : ',nstepsMC
print,''

; storage arrays
nstars = lonarr(n)              ; stores no. of stars in each 30'x30'LOS
count_star = lon64arr(n,5)      ; stores no. of companions for each LOS

FOR i = 0l,n-1 DO BEGIN         ; loop for each LOS (binary)
    
    ra0    = bry[i].ra[0]       ; system properties are subscripted with 0
    dec0   = bry[i].dec[0] 
    theta0 = bry_theta[i]
    dist0 = 0.5*(bry[i].dist[0]+bry[i].dist[1])
    sig_ddist0 = 0.1383*sqrt(bry[i].dist[0]^2+bry[i].dist[1]^2)
    
    vel0 = 0.5*[bry[i].pmra[0]  + bry[i].pmra[1],$
                bry[i].pmdec[0] + bry[i].pmdec[1],$
                bry[i].rv[0] + bry[i].rv[1]]
    sig_vel0 = sqrt([bry[i].pmraerr[0]^2  + bry[i].pmraerr[1]^2,$
                     bry[i].pmdecerr[0]^2 + bry[i].pmdecerr[1]^2,$
                     bry[i].rverr[0]^2 + bry[i].rverr[1]^2])

    temp   = conv_to_galactic(ra0,dec0,dist0,R0,T0,Z0)
    
    ;; **********************************************************
    ;; ********************      CALC PROB    *******************
    ;; **********************************************************
                                ; storage arrays
    count_MC = lon64arr(5)        ; store data for each niter

    ;; count the number of stars in each cell of length cell_size
    nstars[i] = round(count_nstars(ra0,dec0))

    FOR niter = 0l,nstepsMC-1 DO BEGIN
        temp = gen_nstars(ra0,dec0,nstars[i],ra,dec,dist)
            
        theta = angdist(ra0,dec0,ra,dec)
        ddist = abs(dist0-dist) 
        
        ;; ************** COUNT FOR MATCHES  *******************
        ind1 = where(theta GE 0 AND theta LE theta0,c1) ; counts all stars within given theta and all d

        ind2 = where(theta GE 0 AND theta LE theta0 AND $
                     ddist LE sig_ddist0 AND ddist LE 100,c2) ; counts stars within given theta and d

        ;; if kinematics are available
        IF (c2 GT 0 AND (PM EQ 'yes' or RV EQ 'yes')) THEN BEGIN
           vel = gen_pm(R0,T0,Z0,ra0,dec0,dist0,c2) ; returns [[pmra], [pmdec],[rv]]
           
           ;; replicate vel0 to match the dimensions of generated velocity array
           ;; allows for vector arithmetic
           vel0_arr = transpose(cmreplicate(vel0,c2))
           sig_vel0_arr = transpose(cmreplicate(sig_vel0,c2))

           ;; difference in binary and simulated velocity in units of sigma
           dVel = abs(vel - vel0_arr) / sig_vel0_arr

           IF (PM EQ 'yes') THEN $ ;; PM MATCH
              ind3 = where(sqrt(dVel[0]^2 + dVel[1]^2) le 2,c3) ELSE c3 = 0
           IF (RV EQ 'yes') THEN $ ;; RV match
              ind4 = where(dVel[2] le 1,c4)  ELSE c4 = 0
           IF (PM EQ 'yes' and RV EQ 'yes') THEN $ ;; PM+RV match
              ind5 = where(dVel[0]^2+dVel[1]^2 le 2 and dVel[2] le 1,c5)

        ENDIF ELSE BEGIN
           c3 = 0
           c4 = 0
           c5 = 0
        ENDELSE

        ;; ******************** STORE DATA FOR EACH NITER  ********
        count_MC += [c1,c2,c3,c4,c5]

     ENDFOR                     ;; END of ONE MC STEP
;    print,count_MC,format='(5(i6))'
        ; *********************** STORE DATA FOR EACH STAR ***********
    count_star[i,*] = count_MC
    
    ;; print update every 5%
    IF (i MOD round(0.1*n)) EQ 0 THEN BEGIN
        print,systime(),(100.*i/n),'% done. N = ',i+1,FORMAT='(a,f8.1,a,i-5)'
    ENDIF

ENDFOR                          ; END of ONE STAR

prob  = count_star/nstepsMC

chunk=0
IF NOT keyword_set(outfile) THEN $
   outfile='kolby/slw6D.out'

COMMENT='        RA           DEC         DIST     THETA         P1            P2            P3            P4            P5      Nstars'
FORMAT = '(2(f13.6),f12.1,f8.2,5(f14.5),i8)'
   ; ### RA   --> Right Ascension of simulated ellipsoid
   ; ### DEC  --> Declination of simulated ellipsoid
   ; ### DIST --> Distance of simulated ellipsoid
   ; ### P1 --> P(chance alignment | theta)
   ; ### P2 --> P(chance alignment | theta, distance)
   ; ### P3 --> P(chance alignment | theta, distance, mu)
   ; ### P4 --> P(chance alignment | theta, distance, RV)
   ; ### P5 --> P(chance alignment | theta, distance, mu, RV)
   ; ### Nstars --> No. of stars simulated

lun = 1
openw,lun,outfile
printf,lun,COMMENT
print,' '
print,'Printing to... ',outfile
for i=0,n-1 do begin

   printf,lun,bry[i].ra[0],bry[i].dec[0],0.5*(bry[i].dist[0]+bry[i].dist[1]),bry_theta[i],$
          prob[i,0],prob[i,1],prob[i,2],prob[i,3],prob[i,4],nstars[i],$
          format=format
endfor
close,lun

print,'            *************               '
print,'            *************               '
;forprint,bry.ra[0],bry.dec[0],bry.avg_dist,bry.theta,prob[*,0],prob[*,1],nstars,format=format

print,'TOTAL TIME TAKEN   : ',(systime(1)-t_start)/3600.,' hours'
print,'TIME TAKEN PER LOS : ',(systime(1)-t_start)/(60*n),' minutes'
print,'END TIME           : ',systime(0)

END
