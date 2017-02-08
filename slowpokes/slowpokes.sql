-- Examples of a query used in the original slowpokes work.
-- queries the PhotoObjAll and the ProperMotions table in SDSS

CREATE TABLE #upload (up_ra FLOAT, up_dec FLOAT, up_id BIGINT );
INSERT INTO #upload SELECT TOP 2 ra as up_ra, dec as up_dec, objID as up_id FROM MYDB.BooHighPM
CREATE TABLE #tmp (up_id BIGINT, objID BIGINT);
INSERT INTO #tmp EXEC spGetNeighbors 1

INSERT INTO MYDB.HighPMSDSS SELECT
t.up_id as target_id, t.objID,
S.ra, S.dec, pm.pmra, pm.pmdec, pm.pmraerr, pm.pmdecerr,
S.run, S.rerun, S.camcol, S.field, S.obj,
S.psfMag_u, S.psfMag_g, S.psfMag_r, S.psfMag_i, S.psfMag_z,
S.psfMagErr_u, S.psfMagErr_g, S.psfMagErr_r, S.psfMagErr_i, S.psfMagErr_z,
S.extinction_u, S.extinction_g, S.extinction_r, S.extinction_i, S.extinction_z,
S.flags, S.status, S.type, S.SpecObjID,
pm.pmL, pm.pmB, pm.delta, pm.match, pm.nFit, pm.dist22,
pm.sigRa, pm.sigDec, pm.O as usnoO, pm.J as usnoJ

FROM #tmp t, dr6.PhotoObjAll S, dr6.ProperMotions pm, MYDB.BooHighPM Boo
WHERE t.up_id = S.objID and t.up_id = pm.objID and S.objID = pm.objID
ORDER BY Boo.up_id
