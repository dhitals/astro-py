-- query to search around cataloged stars within 30" in the SDSS STAR database
-- assumes a output table has already been created with a separate query.

-- Clean the photometry. This can be omitted and done later in IDL for speed.
@PEAKCENTER=dbo.fPhotoFlags('PEAKCENTER')declare @NOTCHECKED bigint set @NOTCHECKED=dbo.fPhotoFlags('NOTCHECKED')

declare @PSF_FLUX_INTERP bigint set @PSF_FLUX_INTERP=dbo.fPhotoFlags('PSF_FLUX_INTERP')
declare @INTERP_CENTER bigint set @INTERP_CENTER=dbo.fPhotoFlags('INTERP_CENTER')
declare @BAD_COUNTS_ERROR bigint set @BAD_COUNTS_ERROR=dbo.fPhotoFlags('BAD_COUNTS_ERROR')
declare @SATURATED bigint set @SATURATED=dbo.fPhotoFlags('SATURATED')
declare @BRIGHT bigint set @BRIGHT=dbo.fPhotoFlags('BRIGHT')
declare @NODEBLEND bigint set @NODEBLEND=dbo.fPhotoFlags('NODEBLEND')
declare @DEBLEND_NOPEAK bigint set @DEBLEND_NOPEAK=dbo.fPhotoFlags('DEBLEND_NOPEAK')

declare @bad_flags bigint set @bad_flags=(@PEAKCENTER|@NOTCHECKED|@PSF_FLUX_INTERP|@INTERP_CENTER|@BAD_COUNTS_ERROR|@SATURATED|@BRIGHT|@NODEBLEND|@DEBLEND_NOPEAK)

SELECT
  N.objID as targetID, N.NeighborObjID as objID,
  S.ra, S.dec,
  S.psfmag_u,S.psfmag_g,S.psfmag_r,S.psfmag_i,S.psfmag_z,
  S.psfmagerr_u,S.psfmagerr_g,S.psfmagerr_r,S.psfmagerr_i,S.psfmagerr_z,
  S.extinction_r
  INTO MYDB.N01
FROM MYDB.L01 as L, Neighbors as N JOIN STAR as S on S.objID = N.NeighborObjID
WHERE N.objID = L.objID
  and N.mode = 1 and N.neighborMode = 1
  and N.type = 6 and N.neighborType = 6
  and (S.flags_r & @bad_flags) = 0 and (S.flags_i & @bad_flags) = 0 and (S.flags_z & @bad_flags) = 0
  and S.nChild = 0
  and S.psfmagerr_r > 0 and S.psfmagerr_i > 0 and S.psfmagerr_z > 0
  and abs(S.psfmagerr_i) <= 0.10 and abs(S.psfmagerr_z) <= 0.10
  and S.extinction_r <= 0.50
ORDER BY N.objID
