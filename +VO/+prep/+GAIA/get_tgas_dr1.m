function get_tgas_dr1
% Retrieve and format the GAIA-DR1 TGAS catalog
% Package: VO.prep.GAIA
% Description: Retrieve the GAIA DR1 TGAS catalog files and reformat it


LocationURL = 'http://cdn.gea.esac.esa.int/Gaia/tgas_source/fits/';

% get the file names
Links = www.find_urls(LocationURL,'match','.*?\.fits');

% Get files
Names = www.pwget(Links);

% gunzip 
%Names = AstroUtil.files.gunzip(Names);

%hip,tycho2_id,solution_id,source_id,random_index,ref_epoch,ra,ra_error,dec,dec_error,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,astrometric_n_obs_al,astrometric_n_obs_ac,astrometric_n_good_obs_al,astrometric_n_good_obs_ac,astrometric_n_bad_obs_al,astrometric_n_bad_obs_ac,astrometric_delta_q,astrometric_excess_noise,astrometric_excess_noise_sig,astrometric_primary_flag,astrometric_relegation_factor,astrometric_weight_al,astrometric_weight_ac,astrometric_priors_used,matched_observations,duplicated_source,scan_direction_strength_k1,scan_direction_strength_k2,scan_direction_strength_k3,scan_direction_strength_k4,scan_direction_mean_k1,scan_direction_mean_k2,scan_direction_mean_k3,scan_direction_mean_k4,phot_g_n_obs,phot_g_mean_flux,phot_g_mean_flux_error,phot_g_mean_mag,phot_variable_flag,l,b,ecl_lon,ecl_lat

SelectedCol = {'hip','source_id','ref_epoch','ra','dec','ra_error','dec_error','parallax','parallax_error',...
               'pmra','pmra_error','pmdec','pmdec_error','ra_dec_corr',...
               'phot_g_mean_mag'};

% read the files
Nfile = numel(Names);
for Ifile=1:1:Nfile
    Table=FITS.read_table(Names{Ifile});
    TS(Ifile)   = col_select(Table,SelectedCol);
end

TM = merge(TS);
TM = sortrows(TM,'dec');
AstCat.saveh(TM,'GAIA_TGAS_DR1.hdf5');
Util.files.delete_cell(Names);


