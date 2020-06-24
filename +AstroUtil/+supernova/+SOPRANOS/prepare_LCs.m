function [] = prepare_LCs(sn_name,model,VecVs,VecRs,VecMs,VecFrho,VecEbv, Prog)
% Prepare Model light curves for interpolation
% Package: AstroUtil.supernove.SOPRANOS
% Description: Calculate model light curves grid for interpolation  by
%              calcGrid using Sapir & Waxman 2017 shock cooling model.
%              for each band in its input it calculates 1000x1001 grid 
%              spanned by bolometric luminosity and photosphere
%              temperature of the observed flux in the band, take into
%              consdireation the distance, redshift, extinction and filter
%              transmission curve.
% Input  : - The name of the supernova (The observations, the SN redshift
%              and coordinats  are expected to be in a file named <sn_name>_data.mat)
%          - Shock cooling model to use.
%          - Vector of Shock velocity parameter values [cm s^-1]
%          - Vector of Progenitor radius values [R_\sun]
%          - Vector of Ejecta Mass values [M_\sun]
%          - Vector of progenitor envelope parameter values
%          - Vector of Extinction values (E_{B-V}).
%          - Progenitor type ['RSG' or 'BSG']
% Output : - A file named <sn_name>_<model>_LCs.mat
%               
% See also: AstroUtil.supernova.SOPRANOS.calcGrid
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.calcGrid('PTF12gnt',
% Reliable: 2
%--------------------------------------------------------------------------

data = matfile(sprintf('%s_data.mat',sn_name));

filename = sprintf('%s_%s_LCs.mat', sn_name,model);
bands    = data.bands;
redshift = data.redshift;
VecRs   = VecRs(:);
VecVs   = VecVs(:);
VecMs   = VecMs(:);
VecFrho = VecFrho(:);
VecEbv  = VecEbv(:);

save(filename,'VecVs','VecRs','VecMs','VecFrho','VecEbv', 'bands','redshift', 'Prog', 'model', '-v7.3');


switch model
    case {'rw','sw'}
        [~,~,~,~,~,t_BO,t_max,t_opac,t_tr,t_delta]=...
            AstroUtil.supernova.sn_cooling_sw(1,'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);

        [L_BO,Tc_BO]=...
            AstroUtil.supernova.sn_cooling_sw(t_BO*(1+redshift),'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);

        if strcmpi(model,'sw')
            t_end = min(t_opac,t_max);
        else %'rw'
            t_end = min(t_opac,t_delta);
        end

        [L_end,Tc_end]=...
            AstroUtil.supernova.sn_cooling_sw(t_end*(1+redshift),'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);
    case 'msw'
        [~,~,~,~,~,t_min,t_max,times]=...
            AstroUtil.supernova.sn_cooling_msw(1,'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);

        % MSW t_min and t_max are already redshifted
        [L_BO,Tc_BO]=...
            AstroUtil.supernova.sn_cooling_msw(t_min,'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);

        [L_end,Tc_end]=...
            AstroUtil.supernova.sn_cooling_msw(t_max,'Type',Prog,'Rs',VecRs,'Vs',shiftdim(VecVs,-1),'Ms',shiftdim(VecMs,-2),'f_rho',shiftdim(VecFrho,-3),...
            'Model',model);
        t_BO = times.t_min_MSW; % t_min_msw is identical to t_min, but in the SN frame. We use this time to be backward compatible to sn_cooling_sw results.
        t_opac = times.t_opac;
        t_delta = times.t_delta;
        t_tr   = times.t_tr;
end
        

save(filename,'-append', 't_BO', 't_max', 't_opac', 't_delta', 't_tr', 'model');

vecTc = logspace(log10(min(Tc_end(:))),log10(max(Tc_BO(:))),1000).';
vecL  = logspace(log10(min(L_end(:))),log10(max(L_BO(:))),1001);

Dist = AstroUtil.cosmo.lum_dist(redshift);
vecL = vecL./(4.*pi.*(Dist.*constant.pc).^2);

LCsfile = matfile(filename,'Writable',true);
LCsfile.vecTc = vecTc;
LCsfile.vecL  = vecL;
LCsfile.Mag = cell(1,length(bands));
LCsfile.Flux = cell(1,length(bands));

for iband = 1:length(bands)
    fprintf('Calculating Mag table for %s %s band.\n', bands{iband}.instrumentname, bands{iband}.filtername);
    filter    = bands{iband}.filterObj;
    
    parfor iL = 1:length(vecL)
        tmpFile = matfile(sprintf('%stmp_spec_%d.mat',tempdir,iL),'Writable',true);
        L         = vecL(iL);           
        WaveRange = linspace(filter.min_wl,filter.max_wl,101);

        Spec=cell(1,2);
        Spec{1} = WaveRange;
        [~,~,Spec{2}] = AstroUtil.spec.black_body(vecTc./(1+redshift),WaveRange);

        F  = constant.sigma .* (vecTc./(1+redshift)).^4;
        Spec{2}=Spec{2}./F.*L;

        tmpFile.Mag = AstroUtil.spec.synphot(Spec,filter,[],'AB',[],shiftdim(VecEbv,-1),[],'photon');
        tmpFile.Flux = convert.flux(tmpFile.Mag,'AB','cgs/A',filter.pivot_wl_photon,'A');
    end
    
    for iL=1:length(vecL)
        tmpFile = matfile(sprintf('%stmp_spec_%d.mat',tempdir,iL));
        Mag(:,iL,:) = tmpFile.Mag;
        Flux(:,iL,:) = tmpFile.Flux;
        delete(sprintf('%stmp_spec_%d.mat',tempdir,iL));
    end
    
    LCsfile.Mag(1,iband) = num2cell(Mag,1:ndims(Mag));
    LCsfile.Flux(1,iband) = num2cell(Flux,1:ndims(Flux));
    
    clear Mag Flux
end


end