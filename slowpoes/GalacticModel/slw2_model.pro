PRO SLW_MODEL,arg1,arg2,infile=infile,outfile=outfile,chunk=chunk

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

nstepsMC = 2e2                  ; no. Monte Carlo steps

; **************************************************************************
; **************************************************************************
print,'Start Time: ',systime(0)
t_start = systime(1)

IF keyword_set(INFILE) THEN BEGIN
   bry = mrdfits('../bry2.'+INFILE+'.fits',1)
ENDIF ELSE BEGIN
   restore,'bry2.sav' ;; for ACCRE VM since it can't handle long integer call in mrdfits
ENDELSE

n = n_elements(bry)

; **************************************************************************
chunk_size = 5000L

IF keyword_set(chunk) THEN begin
    x0 = (chunk-1)*chunk_size & x1 = x0+chunk_size-1
    if chunk eq 103 then x1 = n-1
    print,chunk,x0,x1
    bry = bry[x0:x1]
end else chunk = '00'

txt = keyword_set(INFILE) ? INFILE:'bry2.chunk'
outfile='slw2Model.'+txt+timestamp(10)+'.dat'

IF n_params() EQ 1 THEN bry = bry[round(randomu(seed,arg1)*n)] ; randomize if not doing the complete simulation
IF n_params() EQ 2 THEN bry = bry[arg1:arg2] ; do the model on a subset
n = n_elements(bry)
; **************************************************************************

print,''
print,'No. of candidate pairs: ',n
print,'No. of MC steps       : ',nstepsMC
print,''

; storage arrays
nstars = lonarr(n)            ; stores no. of stars in each 30'x30'LOS
count_star = lonarr(n,2)       ; stores no. of companions for each LOS

FOR i = 0l,n-1 DO BEGIN      ; loop for each LOS (binary)
    
    ra0    = bry[i].ra[0]       ; system properties are subscripted with 0
    dec0   = bry[i].dec[0] 
    theta0 = bry[i].theta

    dist0 = bry[i].avg_dist
    disterr0 = sqrt(2)*13.83*dist0
    
    temp   = conv_to_galactic(ra0,dec0,dist0,R0,T0,Z0)
    
    ;; **********************************************************
    ;; ********************      CALC PROB    *******************
    ;; **********************************************************
                                ; storage arrays
    count_MC = lonarr(2)        ; store data for each niter

    ; count the number of stars in each cell of length cell_size
    nstars[i] = round(count_nstars(ra0,dec0))

    FOR niter = 0l,nstepsMC-1 DO BEGIN
        temp = gen_nstars(ra0,dec0,nstars[i],ra,dec,dist)
            
        theta = angdist(ra0,dec0,ra,dec)
        ddist = abs(dist0-dist) 
        
         ;;; ************** COUNT FOR MATCHES  *******************
        theta_lim = [0,theta0]
        dist_lim  = dist0 + disterr0 * [-1,1]
        IF dist_lim[0] LT 0 THEN dist_lim[0] = 0.

        ind1 = where(theta GE theta_lim[0] AND theta LE theta_lim[1],count1) ; counts all stars within given theta and all d
        ind2 = where(theta GE theta_lim[0] AND theta LE theta_lim[1] AND $
                     ddist GE dist_lim[0]  AND ddist LE dist_lim[1] AND ddist LE 200,count2) ; counts stars within given theta and d

         ; ******************** STORE DATA FOR EACH NITER  ********
        count_MC += [count1,count2]

    ENDFOR                      ; END of ONE MC STEP
    
        ; *********************** STORE DATA FOR EACH STAR ***********
    count_star[i,*] = count_MC
    
    ;; print update every 5%
    IF (i MOD round(0.1*n)) EQ 0 THEN BEGIN
        print,systime(),(100.*i/n),'% done. N = ',i+1,FORMAT='(a,f8.1,a,i-5)'
        ;spawn,'mail -s "slowpokes model" sauravdhital@gmail.com' < str
    ENDIF

ENDFOR                          ; END of ONE STAR

prob  = count_star/nstepsMC

COMMENT='        RA           DEC         DIST     THETA         P1            P2      Nstars'
FORMAT = '(2(f13.6),f12.1,f8.2,2(f14.5),i8)'
   ; ### RA   --> Right Ascension of simulated ellipsoid
   ; ### DEC  --> Declination of simulated ellipsoid
   ; ### DIST --> Distance of simulated ellipsoid
   ; ### P1 --> P(chance alignment | theta)
   ; ### P2 --> P(chance alignment | theta, distance)
   ; ### Nstars --> No. of stars simulated

lun = 1
openw,lun,outfile
printf,lun,COMMENT
print,' '
print,'Printing to... ',outfile
for i=0,n-1 do begin
   printf,lun,bry[i].ra[0],bry[i].dec[0],bry[i].avg_dist,bry[i].theta,prob[i,0],prob[i,1],nstars[i],$
          format=format
endfor
close,lun

print,'            *************               '
print,'            *************               '
;forprint,bry.ra[0],bry.dec[0],bry.avg_dist,bry.theta,prob[*,0],prob[*,1],nstars,format=format

print,'TOTAL TIME TAKEN   : ',(systime(1)-t_start)/3600.,' hours'
print,'TIME TAKEN PER LOS : ',(systime(1)-t_start)/(60*n),' minutes'
print,'END TIME           : ',systime(0)

;; write P_f to bry FITS file
bry.p_f = transpose(prob)
slw = bry[where(bry.P_f[1] le 0.05)]

mwrfits,slw,'../slw2.'+infile+'.fits',/CREATE
mwrfits,bry,'../bry2.'+infile+'.fits',/CREATE
END
