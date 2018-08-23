function [FilterData,ZP_AB,ZP_Vega]=get_filter(FilterName);
%-------------------------------------------------------------
% get_filter function      Get a filter transmission curves
%                        and zero points.
% Input  : - string containing Filter Name or
%            transmition curve [wavelength (\AA), transmition].
%            Available filters are:
%            'hH' - HST/NICMOS F160W (H-band).
%            'hV' - HST/WFC F555W (V-band).
%            'hR' - HST/WFC2 F675W (R-band).
%            'hI' - HST/WFC F814W (I-band). 
%            'Gf' - GALEX FUV band.
%            'Gn' - GALEX NUV band.
%            'jU' - jhonson U band.
%            'jB' - jhonson B band.
%            'jV' - jhonson V band.
%            'jR' - jhonson R band.
%            'jI' - jhonson I band.
%            'jK' - jhonson K band.
%            'Su' - SDSS u band.
%            'Sg' - SDSS g band.
%            'Sr' - SDSS r band.
%            'Si' - SDSS i band.
%            'Sip'- partial SDSS i band (cutted at about 7600A).
%            'Sz' - SDSS z band.
% Output : - Filter transmission curve [Wavelength (\AA), Trans].
%          - AB mag zero point - assuming flux is in erg/cm^2/s.
%            Add ZP to -2.5*log10(Flux) to get mag
%          - Vega mag zero point - assuming flux is in erg/cm^2/s.
%            Add ZP to -2.5*log10(Flux) to get mag
%-------------------------------------------------------------


switch FilterName
 case {'hH'}
    FilterData = f160w_nicmos_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'hV'}
    FilterData = f555w_wfc_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'hR'}
    FilterData = f675w_wfc2_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'hI'}
    FilterData = f814w_wfc_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;

 case {'Gf'}
    FilterData = galex_fuv_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;

 case {'Gn'}
    FilterData = galex_nuv_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;

 case {'jU'}
    FilterData = jhonson_u_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'jB'}
    FilterData = jhonson_b_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;

    VegaSpec = vega_spec;
    InterpMethod  = 'linear';
    [List1,List2]=eq_sampling(FilterData,VegaSpec,[],InterpMethod);    
    ZP_Vega = 2.5.*log10(trapz(List1(:,1),List1(:,2).*List2(:,2)));

 case {'jV'}
    FilterData = jhonson_v_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'jR'}
    FilterData = jhonson_r_band;

    VegaSpec = vega_spec;
    InterpMethod  = 'linear';
    [List1,List2]=eq_sampling(FilterData,VegaSpec,[],InterpMethod);    
    ZP_Vega = 2.5.*log10(trapz(List1(:,1),List1(:,2).*List2(:,2)));

    ZP_AB         = NaN;   % 2.5.*log10(trapz(FilterData(:,1),FilterData(:,2)))+ 48.6
    %ZP_Vega       = NaN;
 case {'jI'}
    FilterData = jhonson_i_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'jK'}
    FilterData = jhonson_k_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Su'}
    FilterData = sdss_u_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Sg'}
    FilterData = sdss_g_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Sr'}
    FilterData = sdss_r_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Si'}
    FilterData = sdss_i_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Sip'}
    FilterData = sdss_ipart_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 case {'Sz'}
    FilterData = sdss_z_new_band;
    ZP_AB         = NaN;
    ZP_Vega       = NaN;
 otherwise
    error('Unknown FilterName Option');
end
