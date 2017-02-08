pro slw_match,fname=fname
FORWARD_FUNCTION REFORM_ARRAY,FORM_2DARRAY,FORM_3DARRAY

; This routines matches two fits files (a fits file of primary stars
; and a fits file of secondary stars) that have been generated using
; SDSS CasJobs Query within a certain radius. This will output a text
; file containing pairs that have r > 1", d_dist > max_ddist % of
; dist, and d_mu <= max_dmu.

;Fname = ['WD','sdM','VLM']
print,' '
print,'Processing.....L'+fname+'......'+systime()

pri = mrdfits('../L'+fname+'.fit.gz',1)
pri = pri[bsort(pri.objID)]

sec = mrdfits('../N'+fname+'.fit.gz',1)
sec = sec[bsort(sec.targetID)]

priSize = n_elements(pri)

;; get ugriz extinction from extinction_r (Schlegel,Finkbeiner, & Davis 1998)
f = [1.87,1.38,1.00,0.76,0.54]
;; correct mags for extinction
umag1 = pri.psfmag_u - f[0] * pri.extinction_r & umag2 = sec.psfmag_u - f[0] * sec.extinction_r 
gmag1 = pri.psfmag_g - f[1] * pri.extinction_r & gmag2 = sec.psfmag_g - f[1] * sec.extinction_r 
rmag1 = pri.psfmag_r - f[2] * pri.extinction_r & rmag2 = sec.psfmag_r - f[2] * sec.extinction_r 
imag1 = pri.psfmag_i - f[3] * pri.extinction_r & imag2 = sec.psfmag_i - f[3] * sec.extinction_r
zmag1 = pri.psfmag_z - f[4] * pri.extinction_r & zmag2 = sec.psfmag_z - f[4] * sec.extinction_r   
    
theta = angdist(pri.ra,pri.dec,sec.ra,sec.dec) ; calculate the angular distance in arcseconds
   
maxPhotErr = (fname eq 'sdM') ? 0.05:0.10
minTheta = 1.
maxTheta = 20.
maxDdist = 100.

CASE fname OF
   'wd': BEGIN
      ;; calculate the distance assuming it is a dM
      dist1 = dist_ms(gmag1,rmag1,imag1,zmag1)
      
      ;; now calculate the distance to the candidate-WD, which were
      ;; selected from the (u-g, g-r) locus from Girven et al. (2011)
      dist2 = fltarr(priSize)
      for i = 0l,priSize-1 do begin
         mag = [umag2[i],gmag2[i],rmag2[i],imag2[i],zmag2[i]]
         dist2[i] = dist_wd(sec[i].ra,sec[i].dec,mag,sec[i].extinction_r)
      endfor

      ddist = abs(dist1-dist2)
      sig_ddist = 0.1383*sqrt(dist1^2+dist2^2) ; add the errors of the two distances in quadrature
      
      ind2 = where(theta ge minTheta and theta le maxTheta and $
                   pri.psfmagErr_g le maxPhotErr and pri.psfmagErr_r le maxPhotErr and $
                   pri.psfmagErr_i le maxPhotErr and pri.psfmagErr_z le maxPhotErr and $
                   sec.psfmagErr_g le maxPhotErr and sec.psfmagErr_r le maxPhotErr and $
                   sec.psfmagErr_i le maxPhotErr and sec.psfmagErr_z le maxPhotErr and $
                   dist1 gt 0 and dist2 gt 0 and $
                   ddist le sig_ddist and ddist le maxDdist,nCandidateBry)      
   END   
;;
   'sdm': BEGIN
      ;; calculate the distance assuming it is a sdM
      dist1 = dist_sdm(gmag1,rmag1,imag1,zmag1)
      dist2 = dist_sdm(gmag2,rmag2,imag2,zmag2)

      ddist = abs(dist1-dist2)
      sig_ddist = 0.1383*sqrt(dist1^2+dist2^2) ; add the errors of the two distances in quadrature

      ind2 = where(theta ge minTheta and theta le maxTheta and $
                   (gmag1-rmag1) ge 1.5 and (gmag2-rmag2) ge 1.5 and $
                   (rmag1-zmag1) ge 0.8 and (rmag1-zmag1) le 3.6 and $
                   (rmag2-zmag2) ge 0.8 and (rmag2-zmag2) le 3.6 and $
                   pri.psfmagErr_g le maxPhotErr and pri.psfmagErr_r le maxPhotErr and $
                   pri.psfmagErr_i le maxPhotErr and pri.psfmagErr_z le maxPhotErr and $
                   sec.psfmagErr_g le maxPhotErr and sec.psfmagErr_r le maxPhotErr and $
                   sec.psfmagErr_i le maxPhotErr and sec.psfmagErr_z le maxPhotErr and $
                   dist1 gt 0 and dist2 gt 0 and $
                   ddist le sig_ddist and ddist le maxDdist,nCandidateBry)      
   END
;;
   'vlm': BEGIN
       dist1 = dist_ms(gmag1,rmag1,imag1,zmag1)
       dist2 = dist_ms(gmag2,rmag2,imag2,zmag2)

       ddist = abs(dist1-dist2)
       sig_ddist = 0.1383*sqrt(dist1^2+dist2^2) ; add the errors of the two distances in quadrature
      
       ind2 = where(theta ge minTheta and theta le maxTheta and $
                    (imag1-zmag1) ge 1.14 and (imag2-zmag2) ge 1.14 and $
                    dist1 gt 0 and dist1 le 450 and $
                    dist2 gt 0 and dist2 le 450 and $ ;; max dist based on abs_i for M6
                    ddist le sig_ddist and ddist le maxDdist and $
                    pri.psfmagErr_i le maxPhotErr and pri.psfmagErr_z le maxPhotErr and $
                    sec.psfmagErr_i le maxPhotErr and sec.psfmagErr_z le maxPhotErr and $
                    imag1 le 21.3 and imag2 le 21.3 and $ 
                    zmag1 le 20.5 and zmag2 le 20.5,nCandidateBry)
   END
ENDCASE

print,'No. of candidate pairs found: ',nCandidateBry

nComp = 2 ;; two components
slw = replicate({ID:' ',$
                 targetID:0LL,$ 
                 objID:lon64arr(nComp),$
                 specObjID:lon64arr(nComp),$
                 run:intarr(nComp),rerun:intarr(nComp),camcol:bytarr(nComp),field:intarr(nComp),obj:intarr(nComp),$
                 ra:dblarr(nComp),dec:dblarr(nComp),$
                 lgal:dblarr(nComp), bgal:dblarr(nComp),$
                 psfMag:fltarr(nComp,5),psfMagErr:fltarr(nComp,5),extinction:fltarr(nComp,5),$
                 mag:fltarr(nComp,5),$
                 flags:lon64arr(nComp,5),flag_check:bytarr(nComp),$
                 MJD:fltarr(nComp),$
                 nObserve:bytarr(nComp),nDetect:bytarr(nComp),nEdge:bytarr(nComp),$
                 dist:fltarr(nComp), avg_dist:0.,$
                 theta:0., posang:0.,$
                 P_f:fltarr(2),$
                 rz:fltarr(nComp), iz:fltarr(nComp), spType:strarr(nComp), class:' ',$
                 nMulti:0B},nCandidateBry)
                 ;; pmRA_sdss:fltarr(nComp),pmDEC_sdss:fltarr(nComp),$
                 ;; ID_2mass:lon64arr(nComp),$  ;; 2MASS parameters
                 ;; ra_2mass:dblarr(nComp),dec_2mass:dblarr(nComp),$
                 ;; mag_2mass:fltarr(nComp,3),MagErr_2mass:fltarr(nComp,3),$
                 ;; ph_qual:strarr(nComp),rd_flg:strarr(nComp),cc_flg:strarr(nComp),$
                 ;; dist_2mass:fltarr(nComp),$
                 ;; JD_2mass:fltarr(nComp),$
                 ;; ID_ukidss:lon64arr(nComp),$ ;; UKIDSS parameters
                 ;; ra_ukidss:dblarr(nComp),dec_ukidss:dblarr(nComp),$
                 ;; distance_ukidss:fltarr(nComp),$
                 ;; pmRA_ukidss:fltarr(nComp),pmDEC_ukidss:fltarr(nComp),chi2_ukidss:fltarr(nComp),nFrames_ukidss:intarr(nComp),$
                 ;; aperMag_ukidss:fltarr(nComp,5),aperMagErr_ukidss:fltarr(nComp,5),A_ukidss:fltarr(nComp,5),E_BV:fltarr(nComp),$
                 ;; epoch_ukidss:fltarr(nComp),$
                 ;; priORsec_ukidss:lon64arr(nComp)

slw.ra =  form_2Darray(pri[ind2].ra,sec[ind2].ra)
slw.dec = form_2Darray(pri[ind2].dec,sec[ind2].dec)

name = adstring(slw.ra[0],slw.dec[0])
slw.ID = 'SLW'+strmid(name,1,2)+strmid(name,4,2)+strmid(name,13,3)+strmid(name,17,2)
slw.class = fname

slw.targetID = sec[ind2].targetID
slw.objID = form_2Darray(pri[ind2].objID,sec[ind2].objID)

slw.dist = form_2Darray(dist1[ind2],dist2[ind2])
slw.avg_dist = 0.5 * (dist1[ind2]+dist2[ind2])

euler,slw.ra,slw.dec,ll,bb,1
slw.lgal = ll & slw.bgal = bb

slw.theta = angdist(slw.ra[0],slw.dec[0],slw.ra[1],slw.dec[1])
posang,1,slw.ra[0]/15,slw.dec[0],slw.ra[1]/15,slw.dec[1],pa
slw.posang = pa

slw.psfMag = form_3Darray([[pri[ind2].psfmag_u],[pri[ind2].psfmag_g],[pri[ind2].psfmag_r],[pri[ind2].psfmag_i],[pri[ind2].psfmag_z]],$
                          [[sec[ind2].psfmag_u],[sec[ind2].psfmag_g],[sec[ind2].psfmag_r],[sec[ind2].psfmag_i],[sec[ind2].psfmag_z]])    

slw.psfMagErr = form_3Darray([[pri[ind2].psfMagErr_u],[pri[ind2].psfMagErr_g],[pri[ind2].psfMagErr_r],[pri[ind2].psfMagErr_i],[pri[ind2].psfMagErr_z]],$
                             [[sec[ind2].psfMagErr_u],[sec[ind2].psfMagErr_g],[sec[ind2].psfMagErr_r],[sec[ind2].psfMagErr_i],[sec[ind2].psfMagErr_z]])    

f = transpose(cmreplicate(f,nCandidateBry))
Av1 = cmreplicate(pri[ind2].extinction_r,5)
Av2 = cmreplicate(sec[ind2].extinction_r,5)
slw.extinction = form_3Darray(f*Av1,f*Av2)

;; extinction-corrected mags
umag1 = slw.psfmag[0,0] - slw.extinction[0,0] & umag2 = slw.psfmag[1,0] - slw.extinction[1,0] 
gmag1 = slw.psfmag[0,1] - slw.extinction[0,1] & gmag2 = slw.psfmag[1,1] - slw.extinction[1,1] 
rmag1 = slw.psfmag[0,2] - slw.extinction[0,2] & rmag2 = slw.psfmag[1,2] - slw.extinction[1,2] 
imag1 = slw.psfmag[0,3] - slw.extinction[0,3] & imag2 = slw.psfmag[1,3] - slw.extinction[1,3]
zmag1 = slw.psfmag[0,4] - slw.extinction[0,4] & zmag2 = slw.psfmag[1,4] - slw.extinction[1,4]   

slw.mag = form_3Darray([[umag1],[gmag1],[rmag1],[imag1],[zmag1]],$
                       [[umag2],[gmag2],[rmag2],[imag2],[zmag2]])

;; calculate the photometric spectral type
slw.rz[0] = rmag1 - zmag1 & slw.rz[1] = rmag2 - zmag2
slw.iz[0] = imag1 - zmag1 & slw.iz[1] = imag2 - zmag2

slw.sptype[0] = sptype_dm(slw.rz[0],/nosubtype)
slw.sptype[1] = sptype_dm(slw.rz[1],/nosubtype)

mwrfits,slw,'../bry2.'+fname+'.fits',/CREATE 

stop
end

FUNCTION reform_array,a1,a2

if N_params() eq 1 then aout =  a1 else begin   
   sz = size(a1,/N_dimensions)

   if sz gt 2 then begin
      print,'I cannot combine 3-D (or more) arrays. Either fix me or do it yourself'
      return,0
   endif else $
      aout = (sz eq 1) ? form_2Darray(a1,a2):form_3Darray(a1,a2)
endelse

return, aout
END

;; combine two 2-d arrays ([N1,N2])  such that the output is [2,N1,N2]
FUNCTION form_3Darray,a1,a2
return, transpose([[[a1]],[[a2]]])
END

;; combine two 1-d arrays ([N])  such that the output is [2,N]
FUNCTION form_2Darray,a1,a2
return,transpose([[a1],[a2]])
END

;; delete doubles when binary AB/BA exist
FUNCTION delete_doubles,bry

lun = -1

printf,lun,'Deleting the duplicate pairs ...'

bry = bry[bsort(bry.theta)]
n = n_elements(bry)

for i = 0l,n-2 do begin
    if ((bry[i].objID[0] eq bry[i+1].objID[0]) and $
        (bry[i].objID[1] eq bry[i+1].objID[1])) OR $
       ((bry[i].objID[0] eq bry[i+1].objID[1]) and $
        (bry[i].objID[1] eq bry[i+1].objID[0])) then begin
       bry[i+1].theta = 0
    endif
 endfor

ind = where(bry.theta ne 0,count)
printf,lun,'Doubles Deleted: ',n-count
printf,lun,'Pairs left     : ',count

return, bry[ind]
END
