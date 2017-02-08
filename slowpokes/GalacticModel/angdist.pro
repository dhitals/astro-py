function angdist,ra_1,dec_1,ra_2,dec_2

;This function (ra, dec) of two objects in the decimal system
;(degrees) and returns the angular distance between them in
;arcseconds.

ra1 = ra_1*!DTOR  & dec1 = dec_1*!DTOR
ra2 = ra_2*!DTOR  & dec2 = dec_2*!DTOR

dist = sqrt( (abs(ra1-ra2)^2*cos(dec1)*cos(dec2)) + (abs(dec1-dec2))^2 )

return,3600.d*dist*!RADEG
end
