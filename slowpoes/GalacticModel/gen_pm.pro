FUNCTION GEN_PM,R0,T0,Z0,ra0,dec0,dist0,num
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
gal_uvw_pm,pmra,pmdec,rv,U=U,V=V,W=W,/LSR,$
  ra=replicate(ra0,num),dec=replicate(dec0,num),dist=replicate(dist0,num)

RETURN,[[pmra],[pmdec],[rv]]
END
