-- creates a sample of low-mass stars in SDSS.

declare @PEAKCENTER bigint set @PEAKCENTER=dbo.fPhotoFlags('PEAKCENTER')

declare @NOTCHECKED bigint set @NOTCHECKED=dbo.fPhotoFlags('NOTCHECKED')
declare @PSF_FLUX_INTERP bigint set @PSF_FLUX_INTERP=dbo.fPhotoFlags('PSF_FLUX_INTERP')
declare @INTERP_CENTER bigint set @INTERP_CENTER=dbo.fPhotoFlags('INTERP_CENTER')
declare @BAD_COUNTS_ERROR bigint set @BAD_COUNTS_ERROR=dbo.fPhotoFlags('BAD_COUNTS_ERROR')
declare @SATURATED bigint set @SATURATED=dbo.fPhotoFlags('SATURATED')
declare @BRIGHT bigint set @BRIGHT=dbo.fPhotoFlags('BRIGHT')
declare @NODEBLEND bigint set @NODEBLEND=dbo.fPhotoFlags('NODEBLEND')
declare @DEBLEND_NOPEAK bigint set @DEBLEND_NOPEAK=dbo.fPhotoFlags('DEBLEND_NOPEAK')

declare @bad_flags bigint set
@bad_flags=(@PEAKCENTER|@NOTCHECKED|@PSF_FLUX_INTERP|@INTERP_CENTER|@BAD_COUNTS_ERROR|@SATURATED|@BRIGHT|@NODEBLEND|@DEBLEND_NOPEAK)

SELECT objID, ra, dec,
  psfmag_g,psfmag_r,psfmag_i,psfmag_z,
  psfmagerr_g,psfmagerr_r,psfmagerr_i,psfmagerr_z,
  flags_g,flags_r,flags_i,flags_z INTO LMS
FROM Star S
WHERE (flags_r & @bad_flags) = 0 and (flags_i & @bad_flags) = 0 and (flags_z & @bad_flags) = 0
  and nChild = 0
  and (psfmag_r-psfmag_i) >= 0.30
  and (psfmag_i-psfmag_z) >= 0.20
  and psfmagerr_r > 0 and psfmagerr_i > 0 and psfmagerr_z > 0
  and abs(psfmagerr_r) <= 0.10 and abs(psfmagerr_i) <= 0.10 and abs(psfmagerr_z) <= 0.10
  and extinction_r <= 0.50
