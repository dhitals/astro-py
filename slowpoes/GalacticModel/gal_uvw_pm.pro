pro gal_uvw_pm, pmra, pmdec, vrad, U=U, V=V, W=W,$                
                ra=ra, dec=dec, distance=distance, plx=plx, LSR=lsr

if N_Params() EQ 0 then begin
    print,'Syntax - GAL_UVW_PM, PMRA, PMDEC, VRAD, '
    print,'[U=, V=, W=, PLX=, DISTANCE=, /LSR]'
    return
endif

Nra = N_elements(ra)
if (nra EQ 0) or (N_elements(dec) EQ 0) then message, $
  'ERROR - The RA, Dec (J2000) position keywords must be supplied (degrees)'
if (N_elements(U) LT Nra OR N_elements(V) LT Nra OR N_elements(W) LT Nra) then message, $
  'ERROR - UVW space velocities (km/s) must be supplied for each star'
if N_elements(distance) GT 0 then begin 
    bad = where(distance LE 0, Nbad)
    if Nbad GT 0 then message,'ERROR - All distances must be > 0'
    plx = 1e3/distance          ;Parallax in milli-arcseconds
endif else begin
    if N_elements(plx) EQ 0 then message, $
      'ERROR - Either a parallax or distance must be specified'
    bad = where(plx LE 0.0, Nbad)
    if Nbad GT 0 then message,'ERROR - Parallaxes must be > 0'
endelse

; convert to radians
cosd = cos(dec/!RADEG)
sind = sin(dec/!RADEG)
cosa = cos(ra/!RADEG)
sina = sin(ra/!RADEG)

vrad  = fltarr(Nra)
pmra  = fltarr(Nra)
pmdec = fltarr(Nra)

k = 4.74047                     ;Equivalent of 1 A.U/yr in km/s   
t = [ [ 0.0548755604,  0.8734370902,  0.4838350155], $
      [ 0.4941094279, -0.4448296300,  0.7469822445], $
      [-0.8676661490, -0.1980763734,  0.4559837762] ]

for i = 0l,Nra -1 do begin
    a = [ [cosa[i]*cosd[i], -sina[i], -cosa[i]*sind[i] ], $
          [sina[i]*cosd[i],  cosa[i], -sina[i]*sind[i] ], $
          [sind[i]        ,  0      ,  cosd[i]         ] ]
    b = t##a

    uvw = [ [U[i]],[V[i]],[W[i]] ]
    if keyword_set(lsr) then uvw -= [[-9],[12],[7]]

    vec = transpose(b)##uvw 

    vrad[i] = vec[0]         
    pmra[i] = vec[1]*plx[i]/k
    pmdec[i]= vec[2]*plx[i]/k
endfor

sz = size(ra)
if sz(0) EQ 0 then begin
    vrad = vrad[0] & pmra = pmra[0] & pmdec= pmdec[0] 
endif

return
end
