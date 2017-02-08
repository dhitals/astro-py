
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

