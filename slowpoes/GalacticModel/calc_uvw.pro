FUNCTION CALC_UVW,R,theta,Z

Rdot = 0.d 
Tdot = (226.-0.013*Z-1.56e-5*Z^2)/R ; will later  be converted to Tdot(Z) 
Zdot = 0.d                          ; typical values for the MW in km/s

theta = theta*!DTOR

Xdot =   Rdot*COS(theta) - Tdot*SIN(theta)
Ydot = -(Rdot*SIN(theta) + Tdot*COS(theta))

RETURN,[Xdot,Ydot,Zdot]

END
