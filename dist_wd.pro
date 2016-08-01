FUNCTION dist_wd,ra,dec,ugriz,Au,print=print,Teff=xhTeff
  
;; INPUT:           ra, dec, (ugriz mags as 5-D array), Au -- Vectors NOT allowed
;;      OPTIONAL:   If print is set, the calculated parameters will be printed to screen
;; ROUTINES CALLED: EULER, POLY_FIT
;; FILES NEEDED:    bergeronDA8.dat
;;
;; Calculates the distance (and, optionally, mass, log g, age, & Teff) of a 
;; white dwarf (WD) based on its ugriz magnitudes by comparing to the 
;; Bergeron models using a chi-sq minimization technique. Hydrogen-dominated (DA) 
;; atmospheres are assumed.
;; 
;; Input needs to be one star at a time
;; 
;; Reddening (E) and extinction (A) along the line-of-sight calculated using an iterative process.
;;
;; Modification History:
;;
;; 03/08/2011: Changed code structure for efficient and elegance by SD
;; 06/27/2008: Written in IDL by Saurav Dhital
;; Written in FORTRAN by Hugh Harris (photdistDA.f)

if (N_Params() NE 4 OR $
   n_elements(ra) NE 1 OR n_elements(dec) NE 1 OR n_elements(Au) NE 1 OR $
   n_elements(ugriz) ne 5) then begin
   print,'Syntax - dist_wd, ra, dec, ugriz, Au, /print'
   print,'                  ra, dec, Au are scalars'
   print,'                  ugriz is a 5-D array'
   return,0
endif

;; correction factors to interchage between mags/colors
E_col = [0.030,0.015,0.015,0.030] ;; reddening: E(u-g), E(g-r), E(r-i), E(i-z)
f_mag = [0.736,0.534,0.405]       ;; extinction conversion factor for mag --- A_[gri] = f_mag * Au
f_col = [0.264,0.202,0.129,0.119] ;; reddening  correction factor for color --- E_[u-g, g-r, r-i, i-z] = f_col * Au

;; read the Bergeron Hydrogen models --- prefix 'h' -> Hydrogen model
;; If you want to solve for a DB, just read in bergeronDB8.dat although I am not sure 
;; if the solved parameters will be accurate---there is no reason they would not be.
;;
;; columns: Teff,logg,mass,Mbol,BC,UBVRIJHKugrizy,b-y,u-b,v-y,V-I,G-R,U-V,U-G,B-V,age
PATH = '~/estrella/local/datalib/wd/'
readcol,PATH+'bergeronDA8.dat',hTeff,hlogg,hmass,hMbol,hBC,$
        huu,hbb,hvv,hrr,hii,hjj,hhh,hkk,hu,hg,hr,hi,hz,hy,hby,hub,hvy,hvi,hgr,huv,hug,hbv,hage,$
        format='f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f',/silent
Hmag = transpose([[hu],[hg],[hr],[hi],[hz]]) ;; => [5,n] array
Hcol = transpose([[hu-hg],[hg-hr],[hr-hi],[hi-hz]]) ;; => [4,n] array

ind = where(hTeff ge 3500.,size)   ;; we'll only deal with WD Teff >= 3500 K
hTeff = hTeff[ind] & hMbol = hMbol[ind]
Hmag = Hmag[*,ind] & Hcol = Hcol[*,ind]

;; change the coords to Galactic coords
euler,ra,dec,ll,bb,1

;; ********** Start the iterative loop here ************
;; start with a guess for reddening. Get approx temp. Iterate...
fracabs = 0.5
niter = 0
chi1 = 999999. & chimin = chi1-1.

;; Input mags and colors --- Prefix 'I' -> Input
Imag = ugriz
Icol = transpose([ [Imag[0,*]-Imag[1,*]],[Imag[1,*]-Imag[2,*]],$
                   [Imag[2,*]-Imag[3,*]],[Imag[3,*]-Imag[4,*]] ])

;; We'll only use the gri mags
Imag = Imag[1:3,*] & Hmag = Hmag[1:3,*]

;; Now, we shall solve for fracabs and iterate at least 20 times
WHILE (chimin lt chi1 and niter le 20) DO BEGIN
   ;; calc the extinction corrected mags/colors --- Prefix 'x' ->  extinction-corrected
   Xmag = Imag - fracabs*f_mag*Au
   Xcol = Icol - fracabs*f_col*Au

   ;; find extinction law as in k = A_u/E(u-g)  --- Au = data-model
   ;; Note Hcol is a 2D array while Xcol and E_col are 1D
   k = (cmreplicate(Xcol,size)-Hcol)/cmreplicate(E_col,size)

   chi = sqrt(total(k^2,1))       ;; chi = rms (k) i.e. diff between data & model
   index = where(chi eq min(chi))

   ;; FIND A MINIMUM CHI TO GET TEMP

   ;; send in the min and the two adjacent values to polyfit - make sure  min is not at the endpoints
   if index eq 0 then index += 1 
   if index eq size-1 then index -= 1
   ind = [index-1,index,index+1]

   ;; fit a second-order polynomial for chi(hTeff) --- returns coefficients aa = [a0,a1,a2]
   aa = poly_fit(hTeff[ind],chi[ind],2)
   xhTeff = -0.5 * aa[1]/aa[2]  ;; prefix 'xh' -> interpolated --- cross between model & data
   chimin = aa[0] + aa[1]*xhTeff + aa[2]*xhTeff^2

   ;; INTERPOLATE TO GET ABS MAG
   if xhTeff lt hTeff[0] then xhTeff = hTeff[0] ; temp cannot be less than min temp from Bergeron model
   if xhTeff gt hTeff[size-1] then xhTeff = hTeff[size-1] ; temp cannot be more than max temp from Bergeron model

   ind = where(xhTeff le hTeff)
   ind = (ind[0] eq size-1) ? ind[0]-1:ind[0]

   frac = (xhTeff-hTeff[ind])/(hTeff[ind+1]-hTeff[ind])
   xhMbol = hMbol[ind] + frac*(hMbol[ind+1]-hMbol[ind])   
   XHmag = Hmag[*,ind] + frac*(Hmag[*,ind+1]-Hmag[*,ind])

   dmag = MEAN(Xmag-XHmag) ;; avg difference in extinction corrected mag (data-interpolated)
   dist = 10.^((dmag+5)/5.)       ;; distance to WD
   z    = dist*sin(bb*!DTOR)      ;; Galactic height of WD

   ;; calculate fracabs for the next iteration
   if chimin lt chi1 then begin
      if dist lt 100. then fracabs = 0. ; if dist<100 pc, reddening = 0
      if z gt 250. then fracabs = 1.    ; if z>250 pc, reddening = 1
      if (dist ge 100. and z le 250.) then fracabs = (dist-100.)/(250.*(1-1/sin(bb*!DTOR))) ; for intermediate distances, get a new reddening
      chi1 = chimin
   endif

   niter += 1
ENDWHILE

;; print the calculated parameters -- age, log g can be directly infered from 
;; Bergeron grids once you have Teff, Mbol
if keyword_set(print) then $
   print,ra,dec,chimin,niter,dist,z,fracabs,xhMbol,xhTeff

return,dist
end
