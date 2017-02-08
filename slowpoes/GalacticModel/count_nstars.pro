FUNCTION COUNT_NSTARS,ra,dec
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

ddist = 5                       ; steps in distance in pc
n = max_dist/ddist + 1
dist = findgen(n)*ddist         ; 0 < d < 2500 in 5 pc steps
rho  = fltarr(n) & nstars = fltarr(n) ; create an array to store rho FOR each d
 
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
