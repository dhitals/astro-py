FUNCTION CONV_TO_EQUATORIAL,R,theta,Z,ra,dec,d
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

dcosb = SQRT(R^2+Rsun^2-2*R*Rsun*COS(theta*!DTOR))

l = ASIN(R*SIN(theta*!DTOR)/dcosb) * !RADEG
IF l LT 0 THEN BEGIN                          ; 3rd or 4th quadrant
   IF cos(l) GT 0 THEN l = 180-l ELSE l=360+l ; 3rd else 4th
ENDIF

IF l GE 0 THEN BEGIN                          ; 1st or 2nd quadrant
   IF cos(l) GT 0 THEN l=l ELSE l=180-l       ; 1st else 2nd
ENDIF

b = ABS(ATAN((Z-Zsun)/dcosb) + ATAN(Zsun/Rsun)) * !RADEG
IF z LT 0 THEN b = -b
d = dcosb/COS(b*!DTOR) 
euler,l,b,ra,dec,2

RETURN,0
END
