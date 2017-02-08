FUNCTION GEN_2DGAUSSIAN,mu,sig1,sig2,f1,f2,num

n_acc = 0l                ; number of stars accepted
WHILE n_acc LT num DO BEGIN
    n_lft = num - n_acc ; no. of needed stars
    
    X = randomu(seed,n_lft)*10*sig2 - 5*sig2 ; generate random numbers between (-5,5) sigma of a normalized gaussian

;    X = num_gen(-5*sig2,5*sig2,0.1)
    z1 = (X-mu)/sig1 & z2 = (X-mu)/sig2
    G1 = 1/(sqrt(2*!dpi)*sig1)*exp(-z1^2/2) ; thin disk
    G2 = 1/(sqrt(2*!dpi)*sig2)*exp(-z2^2/2) ; thin disk
 
    Px = f1*G1+f2*G2
    Fx = 1.2*max(Px)

    rand = randomu(seed,n_lft)
    ind = where(rand LT Px/Fx, count) ; uniformm deviate comp function

    IF count NE 0 THEN BEGIN
        IF n_acc EQ 0 THEN Xarr = x[ind] ELSE Xarr = [Xarr,x[ind]]
        n_acc += count
    ENDIF
ENDWHILE

RETURN,Xarr
END
