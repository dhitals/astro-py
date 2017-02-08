pro gal_uvwerr, U, V, W, sig_u, sig_v, sig_w, $
             pmra=pmra,pmdec=pmdec,vrad=vrad,$
             sig_pmra=sig_pmra, sig_pmdec=sig_pmdec, sig_vrad=sig_vrad,$
             ra=ra, dec=dec, distance=distance, sig_distance=sig_distance, $
             plx=plx, LSR=lsr, sig_plx = sig_plx
;+
; NAME:
;     GAL_UVW_ERROR
; PURPOSE:
;     Calculate the UVW space velocities and the associated errors
;     using the error propagation derived in Johnson and Soderblom (1987, AJ, 93,864)
; EXPLANATION:
;     Calculates the galactic velocities U, V, and W and their
;     standard deviations for a star given its (1) coordinates, (2)
;     distance (or parallax), (3) proper motion, (4) radial velocity,
;     (5) errors on distance (or parallax), proper motion, and radial velocity
;
; CALLING SEQUENCE:
;     GAL_UVW_ERROR, U, V, W, SIG_U, SIG_V, SIG_W, [DISTANCE=, RA=, DEC=,
;     PLX=, SIG_PMRA=, SIG_PMDEC=, SIG_VRAD=, PMRA=, PMDEC=, 
;     SIG_DISTANCE=, SIG_PLX= ] 
; OUTPUT PARAMETERS:
;     u (km/s), v (km/s), w (km/s)
;     sig_u (km/s)
;     sig_v (km/s)
;     sig_w (km/s)
;
; REQUIRED INPUT KEYWORDS:
;      User must supply a position, galactic space velocity, and distance
;      (or parallax).    Either scalars or vectors can be supplied.
;     (1) Position:
;      RA - Right Ascension in *Degrees*
;      DEC - Declination in *Degrees*
;     (2) Distance or Parallax
;      PLX - parallax with same distance units as proper motion measurements
;            typically milliarcseconds (mas)
;     (3) Proper Motion
;      PMRA - proper motion in the direction of right ascension (mas/yr)
;      PMDEC - proper motion in the direction of declination (mas/yr)
;     (4) Errors
;      SIG_DISTANCE - error on distance (pc)
;                 or
;      SIG_PLX - error on parallax with sam units as proper motion
;                  measurements, typically milliarcseconds (mas)
;      SIG_PMRA - error on proper motion in the direction of right ascension (mas/yr)
;      SIG_PMDEC - error on proper motion in the direction of declination (mas/yr)
;      SIG_VRAD - error on radial velocity (km/s)
;
; METHOD:
;      Follows the error propagation outline of Johnson & Soderblom 
;      (1987, AJ, 93,864) 
;      except that the J2000 transformation matrix to Galactic
;      coordinates is taken from the introduction to the Hipparcos catalog.   
; REVISION HISTORY:
;      Modified, S. Dhital                        September  2008
;                           Incorporated calculation of UVW as well
;      Written, A. Skemer                         May        2004
;      Using code from gal_uvw:
;      Written, W. Landsman                       December   2000
;-

if N_Params() EQ 0 then begin
    print,'Syntax - GAL_UVW_ERROR, SIG_U, SIG_V, SIG_W, '
    print, '[DISTANCE=, RA=, DEC=, PLX=, SIG_PMRA=, SIG_PMDEC=,'
    print, 'SIG_VRAD=, PMRA=, PMDEC=, SIG_DISTANCE=, SIG_PLX=] 
     print,' SIG_U, SIG_V, SIG_W - output Galactic space velocitiy error' 
    return
endif

Nra = N_elements(ra)
if (nra EQ 0) or (N_elements(dec) EQ 0) then message, $
  'ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)'
if (N_elements(sig_pmra) LT Nra) or (N_elements(sig_pmdec) LT Nra) or $
  (N_elements(sig_vrad) LT Nra) then message, $ 
  'ERROR - Uncertainties on proper motion and radial velocity must be supplied for each star'
if N_elements(distance) GT 0 then begin 
    bad = where(distance LE 0, Nbad)
    if Nbad GT 0 then message,'ERROR - All distances must be > 0'
    plx = 1e3/distance          ;Parallax in milli-arcseconds
    sig_plx=1e3*sig_distance/(distance^2) ;Parallax error in mas
endif else begin
    if N_elements(plx) EQ 0 then message, $
      'ERROR - Either a parallax or distance must be specified'
    if N_elements(sig_plx) EQ 0 then message, $
      'ERROR - Uncertainties in parallax or distance must be specified'
    bad = where(plx LE 0.0, Nbad)
    if Nbad GT 0 then message,'ERROR - Parallaxes must be > 0'
endelse

cosd = cos(dec/!RADEG)
sind = sin(dec/!RADEG)
cosa = cos(ra/!RADEG)
sina = sin(ra/!RADEG)

u = fltarr(Nra) & v = u       & w = u
sig_u = u     & sig_v = u & sig_w = u

k = 4.74047                     ;Equivalent of 1 A.U/yr in km/s   

t = [ [ 0.0548755604,  0.8734370902,  0.4838350155],$
      [ 0.4941094279, -0.4448296300,  0.7469822445],$
      [-0.8676661490, -0.1980763734,  0.4559837762] ]

for i = 0,Nra -1 do begin
    a = [ [cosa[i]*cosd[i] ,-sina[i]        ,-cosa[i]*sind[i] ],$
          [sina[i]*cosd[i] , cosa[i]        ,-sina[i]*sind[i] ],$
          [sind[i]         , 0              , cosd[i]         ] ]
    b = t##a

    vec = [[vrad[i]], [k*pmra[i]/plx[i]], [k*pmdec[i]/plx[i]] ]

    uvw = b##vec

    if keyword_set(lsr) then uvw += [[-9],[12],[7]]
    
    error_uvw=sqrt((b^2)##[[sig_vrad [i]^2], $
                       [(k/plx[i])^2*(sig_pmra[i]^2+(pmra[i]*sig_plx[i]/plx[i])^2)],    $
                       [(k/plx[i])^2*(sig_pmdec[i]^2+(pmdec[i]*sig_plx[i]/plx[i])^2)] ] $
                         +2*pmra[i]*pmdec[i]*k^2*sig_plx[i]^2/plx[i]^4* [ [b[1,0]*b[2,0]],[b[1,1]*b[2,1]],[b[1,2]*b[2,2]] ] )
    
    u[i] = uvw[0]  & sig_u[i] = error_uvw[0]
    v[i] = uvw[1]  & sig_v[i] = error_uvw[1]
    w[i] = uvw[2]  & sig_w[i] = error_uvw[2]
endfor

sz = size(ra)
if sz(0) EQ 0 then begin
    u = u[0] & sig_u = sig_u[0]  
    v = v[0] & sig_v = sig_v[0]
    w = w[0] & sig_w = sig_w[0]
endif

return
end
