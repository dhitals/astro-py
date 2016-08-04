FUNCTION GEN_NSTARS,ra0,dec0,num,ra,dec,dist
COMMON CONSTS,Rsun,Tsun,Zsun,rho0,au2pc,cell_size,max_dist,seed

; ra0,dec0,num - input parameters
; ra,dec,dist - output arrays for the generated stars

seed = ranseed()
n_acc = 0l                ; number of stars accepted

WHILE n_acc LT num DO BEGIN
    n_lft = num - n_acc ; no. of needed stars
    
    ra1   = ra0 + (randomu(seed,n_lft) - 0.5) * cell_size
    dec1  = dec0+ (randomu(seed,n_lft) - 0.5) * cell_size
    dist1 = randomu(seed,n_lft) * max_dist

    temp = conv_to_galactic(ra1,dec1,dist1,R,T,Z)

    rho = calc_rho(R,Z)

    rand = randomu(seed,n_lft)

    ; accept if random number is less than rho(R,Z)/rho0
    ind = where(rand LT rho/rho0, count)

    IF count NE 0 THEN BEGIN
        IF n_acc EQ 0 THEN BEGIN
            ra   = ra1[ind]
            dec  = dec1[ind]
            dist = dist1[ind]
        ENDIF ELSE BEGIN
            ra   = [ra,ra1[ind]]
            dec  = [dec,dec1[ind]]
            dist = [dist,dist1[ind]]
        ENDELSE
        n_acc += count
    ENDIF
ENDWHILE

RETURN,0
END
