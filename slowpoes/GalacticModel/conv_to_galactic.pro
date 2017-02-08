FUNCTION CONV_TO_GALACTIC,ra,dec,d,r,t,z
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

euler,ra,dec,l,b,1

r = SQRT( (d*COS(b*!DTOR))^2 + Rsun *(Rsun- 2*d*COS(b*!DTOR)*COS(l*!DTOR)))
t = ASIN(d*SIN(l*!DTOR)*COS(b*!DTOR)/R) * !RADEG
z = Zsun + d * SIN(b*!DTOR - ATAN(Zsun/Rsun))

RETURN, 0
END
