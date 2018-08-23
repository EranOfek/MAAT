%
% Contents file for package: Catalogue
% Created: 29-Dec-2015
%---------
% astrolinks.m :  Given a list of coordinates, return the URL links to the following webpages: SDSS chart, SDSS object, NED, FIRST, NVSS, and DSS. Moreover, generate an html page with the catalog data and links to the verious web pages.
% build_htm4cat.m :  Given an astronomical catalog (in a structure format, e.g., 'FIRST.mat') construct a directory containing catalogs per HTM (Hierarchical Triangular Mesh) region, and an index catalog of all HTM regions.
% cat_match.m :  Match two catalogs by object positions. If the catalogs are not sorted then the program will sort them.
% cat_search.m :  Search a sorted astronomical catalog for all objects found within a given distance from the search position.
% get_adata.m :  Get 'everything' for a given celestial poistion.
% get_apass.m :  Search the APASS (AAVOS photometric all sky survey) catalog around a given coordinate.
% get_batse_lc.m :  Get a 4 chanel BATSE light curve in 64ms bins (only PREB sample) from a local catalog.
% get_cat.m :  Search an astronomical catalog structure by coordinates.
% get_dss.m :  Get link to and the FITS image of a digital sky survey image (POSS-I/II, UKST survey). Read the FITS file into matlab.
% get_first.m :  Search the FIRST (21 cm radio survey) catalog around a given coordinate.
% get_hstsrc.m :  Search a retrieve the HST source catalog around given coordinate.
% get_nvss.m :  Search the NVSS (21 cm radio survey) catalog around a given coordinate.
% get_orbit_files.m :  Get asteroids and comets orbital elements from JPL and read into a matlab structure.
% name_server_ned.m :  Resolve an astronomical object name into coordinates using NASA Extragalactic Database (NED).
% name_server_simbad.m :  Resolve an astronomical object name into coordinates using SIMBAD database.
% search_cat.m :  Given a catalog with Long,Lat coordinates position, search for lines near a list of reference positions. This function can be used to search for a near(est) position in a catalog or to match two catalogs. This function replaces cat_search.m and cat_match.m
% stellar_tracks.m :  Given an initial mass and metllicity return the Geneva stellar tracks as a function of time.
% tile_the_sky.m :  Tiling the celestial sphere with approximately equal area tiles.
% vizquery_path.m :  Return the path of the vizquery/cdsclient directory Edit this program before using wget_2mass.m, wget_usnob1.m and wget_ucac4.m.
% wget_2mass.m :  Query the 2MASS catalog using the VizieR web service. Installation: 1. install cdsclient (instructions can be found in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html) in $USER/matlab/fun/bin/vizquery/cdsclient-3.71/ 2. If you installed the cdsclient in a different location, then edit the first few lines of the code accordingly. This program is replacing search2mass.m
% wget_sdss.m :  Query SDSS PhotoPrimary table around specific coordinate. See run_sdss_sql.m for a more general queries.
% wget_ucac4.m :  Query the UCAC4 catalog using the VizieR web service. Installation: 1. install cdsclient (instructions can be found in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html) in $USER/matlab/fun/bin/vizquery/cdsclient-3.71/ 2. If you installed the cdsclient in a different location, then edit the first few lines of the code accordingly.
% wget_usnob1.m :  Query the USNO-B1 catalog using the VizieR web service. Installation: 1. install cdsclient (instructions can be found in: http://cdsarc.u-strasbg.fr/doc/cdsclient.html) in $USER/matlab/fun/bin/vizquery/cdsclient-3.4/ 2. If you installed the cdsclient in a different location, then edit the first few lines of the code accordingly. This program is replacing search2mass.m
