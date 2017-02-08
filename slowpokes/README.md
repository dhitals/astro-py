### SLoWPoKES (Dhital et al. 2010, 2015)

This is a short guide to using the slowpokes code.

* SQL: `slw2.sql` is the main query. It was used in the CasJObs
  interface on the SDSS website. this searches for all stellar objects
  within 30 arcseconds of your in out table. Currently set up for
  cleaning the photometry. The `slowpokes.sql` was the query used in
  slowpokes-I that matched the SDSS `PhotoObjAll` and `ProperMotions`
  tables.

* Matching the candidate binaries (`slw_match.pro`): IDL script that
  matches a list of low-mass stars with another list of candidate
  companions (the result of the above SQL query). Assumes FITS format
  for both files. Uses sorted lists (and counters) to speed up the
  matching.

  Uses other custom IDL scripts, so let me know which ones are missing.
  
* Galactic Model: IDL code to calculate the probability that a
  candidate binary is a chance alignment. See the SLoWPoKES paper for
  a detailed discussion.
