;; clean the riz photometry for both components
;PRO cleanPhot_lms,bry

bry = mrdfits('../data/LMS.all.fit.gz',1)

size = n_elements(bry) & sz = 1e4
nLoops = ceil(bry/sz)

f = intarr(size)

FOR i=0,nLoops-1 do begin

   x0 = i*sz & x1 = (i lt nLoops-1) ? (i+1)*sz - 1:(size mod sz)
   slw = bry[x0:x1]

   f1 = intarr(sz) &   f2 = intarr(sz) &   f3 = intarr(sz) &   f4 = intarr(sz)
   
   f1[where(slw.psfmag_g le 22.2 and slw.psfmagerr_g le 0.05 and $
            ((slw.flags_g and 2LL^5)  + (slw.flags_g and 2LL^18) + (slw.flags_g and 2LL^19) + $
             (slw.flags_g and 2LL^40) + (slw.flags_g and 2LL^44) + (slw.flags_g and 2LL^47) EQ 0))] = 8
   f2[where(slw.psfmag_r le 22.2 and slw.psfmagerr_r le 0.05 and $
            ((slw.flags_r and 2LL^5)  + (slw.flags_r and 2LL^18) + (slw.flags_r and 2LL^19) + $
             (slw.flags_r and 2LL^40) + (slw.flags_r and 2LL^44) + (slw.flags_r and 2LL^47) EQ 0))] = 4
   f3[where(slw.psfmag_i le 21.3 and slw.psfmagerr_i le 0.05 and $
            ((slw.flags_i and 2LL^5)  + (slw.flags_i and 2LL^18) + (slw.flags_i and 2LL^19) + $
             (slw.flags_i and 2LL^40) + (slw.flags_i and 2LL^44) + (slw.flags_i and 2LL^47) EQ 0))] = 2
   f4[where(slw.psfmag_z le 20.5 and slw.psfmagerr_z le 0.05 and $
            ((slw.flags_z and 2LL^5)  + (slw.flags_z and 2LL^18) + (slw.flags_z and 2LL^19) + $
             (slw.flags_z and 2LL^40) + (slw.flags_z and 2LL^44) + (slw.flags_z and 2LL^47) EQ 0))] = 1
   
   f[x0:x1] = f1+f2+f3+f4
endfor

stop
add_tag, slw, 'flag_check', 0B, slw1, structyp=structyp

mwrfits,slw1,'../data/LMS.clean.fits',/CREATE

END
