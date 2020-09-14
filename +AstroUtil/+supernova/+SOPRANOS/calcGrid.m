function calcGrid(sn_name,Vect0,model,ModelError,BGbands,bg)
% Calculate chi2/dof grid
% Package: AstroUtil.supernove.SOPRANOS
% Description: Calculate chi2/dof grid for a SN light curve against a
%              Sapir & Waxman 2017 shock cooling model.
% Input  : - The name of the supernova (The observations are expected to be
%                   in a file named <sn_name>_data.mat and model light
%                   curves for interpolation in a file named
%                   <sn_name>_<model>_LCs.mat)
%          - Vector of reference times.
%          - Vector of Extinction values (E_{B-V}) - An obselete parameter,
%                               VecEbv is now a part of the LCs file.
%          - Shock cooling model to use. 
%          - The model accuracy - default is 0.
%          - Ordinal number of the bands for which to calculate the BG out
%            of the data
% Output : - A file named <sn_name>_<model>_grid.mat
%               the output file includes the vector which span the grid
%               chi2 and number of valid data points grids
%               cell array with chi2 and number of valid data points grids
%               for each band.
%               cell array with chi2 and number of outliers data points grids
%               for each band, in case outlier points were icnluded.
%               
% See also: AstroUtil.supernova.SOPRANOS.prepare_LCs
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.calcGrid('PTF12gnt',
% Reliable: 2
%--------------------------------------------------------------------------
if nargin==3
    BGbands = [];
    ModelError = 0;
elseif nargin==4
    BGbands = [];
end

data = matfile(sprintf('%s_data.mat',sn_name));
LCs  = matfile(sprintf('%s_%s_LCs.mat', sn_name, model));
gridfile = matfile(sprintf('%s_%s_grid', sn_name, model),'Writable',true);

Ebv = AstroUtil.spec.sky_ebv(data.RArad ,data.decRad);
if (LCs.VecEbv(1,1)>Ebv)||(LCs.VecEbv(end,1)<Ebv)
    warning('The local galaxy Ebv is outside VecEbv');
end

if ~strcmpi(LCs.model, model)
    warning('The input model and the one in %s_%s_LCs.mat are different', sn_name, model);
end

data_vars = who(data);
if any(strcmp(data_vars,'TransientStart'))&&any(strcmp(data_vars,'TransientEnd'))
    transientStart = data.TransientStart;
    transientEnd   = data.TransientEnd;
    gridfile.transientStart = transientStart;
    gridfile.transientEnd   = transientEnd;
    transient      = true;
else
    transient      = false;
end

gridfile.Vect0   = Vect0(:);
gridfile.VecVs   = LCs.VecVs;
gridfile.VecRs   = LCs.VecRs;
gridfile.VecMs   = LCs.VecMs;
gridfile.VecFrho = LCs.VecFrho;
gridfile.VecEbv  = LCs.VecEbv;
gridfile.ModelError = ModelError;
gridfile.BGbands = BGbands;
gridfile.redshift = LCs.redshift;
gridfile.ProgType = LCs.Prog;
gridfile.model    = LCs.model;

Nv  = length(LCs.VecVs);
Nr  = length(LCs.VecRs);
Nm  = length(LCs.VecMs);
Nfrho= length(LCs.VecFrho);
NEbv= length(LCs.VecEbv);
Nt0 = length(Vect0);
t_BO   = LCs.t_BO.*(1+data.redshift);

switch LCs.model
    case 'sw'
        t_max  = min(LCs.t_max,LCs.t_opac).*(1+data.redshift);
    case 'rw'
        t_max  = min(LCs.t_delta,LCs.t_opac).*(1+data.redshift);
    case 'msw'
        t_max = LCs.t_max;
end

bands = data.bands;
band = cell(1,1);
gridfile.bandschi2 = cell(size(bands));
gridfile.bandspnts = cell(size(bands));
gridfile.bandsOutliersChi2 = cell(size(bands));
gridfile.bandsOutliersPnts = cell(size(bands));
gridfile.bg        = cell(size(bands));
gridfile.bands = LCs.bands;

gridfile.chi2   = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'single');
gridfile.points = zeros(Nt0,Nr,Nv,Nm,Nfrho,'uint16');

if transient
    gridfile.NTransient = cell(size(bands));
end
t_BO = shiftdim(t_BO,-1);
t_max = shiftdim(t_max,-1);


for iband=1:length(bands)
    fprintf('iband=%d\n',iband)
    LC = sortrows(unique(bands{iband}.LC),'MJD');
    Nobs = length(LC.MJD);
    [MJDobs, t0obs] = ndgrid(LC.MJD, Vect0);
    W.tobs = MJDobs - t0obs;
    
    if transient
        W.TransientObs = (LC.MJD>=transientStart)&(LC.MJD<=transientEnd);
    end
    clear MJDobs t0obs
   
    W.Prog = LCs.Prog;
    W.VecRs = LCs.VecRs;
    W.VecVs = LCs.VecVs;
    W.VecMs = LCs.VecMs;
    W.VecFrho = LCs.VecFrho;
    W.redshift = LCs.redshift;
    W.model = LCs.model;
    W.dist = AstroUtil.cosmo.lum_dist(LCs.redshift);
    W.WaveRange = (bands{iband}.filterObj.min_wl: (bands{iband}.filterObj.max_wl - bands{iband}.filterObj.min_wl)./100 :bands{iband}.filterObj.max_wl);
    W.filter = bands{iband}.filterObj;
    W.model = LCs.model;

    chi2  = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'single');
    points = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'uint16');
    chi2outliers  = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'single');
    pointsOutliers = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'uint16');
    if transient
        transientPoints = zeros(Nt0,Nr,Nv,Nm,Nfrho,'uint16');
    end
    if ismember(iband,BGbands)
        BG     = zeros(Nt0,Nr,Nv,Nm,Nfrho,NEbv,'single');
    elseif isfield(bands{iband},'background')
        BG     = bands{iband}.background;
    else
        BG     = 0;
    end
    
    fprintf('calculating L & T:');
    for it0 = 1:Nt0
        fprintf('\nit0=%d, iEbv=',it0);
        for iEbv=1:NEbv
            maxdiff=zeros(length(LC.MJD),1);
            fprintf('%d,',iEbv);
            loopfileName = sprintf('%s_iEbv%d.mat',sn_name,iEbv);
%            loopfile = matfile(loopfileName,'Writable',true);
            LTcgrid = cell2mat(LCs.Flux(1,iband));
            LTcgrid = LTcgrid(:,:,iEbv).';
            W.F = griddedInterpolant({LCs.vecL,LCs.vecTc},LTcgrid);
            clear LTcgrid;
            W.Ebv = LCs.VecEbv(iEbv,1);
            c=parallel.pool.Constant(W);
            Flux = zeros([length(LC.MJD), Nr, Nv, Nm, Nfrho]);
            if transient
                NTransient = zeros(Nr, Nv, Nm, Nfrho);
            end
            parfor iObs=1:length(LC.MJD)
                validObs = shiftdim((c.Value.tobs(iObs,it0)>=t_BO) & (c.Value.tobs(iObs,it0)<=t_max),1);
                validTransObs = validObs & c.Value.TransientObs(iObs);
                tmpLC  = zeros(Nr, Nv, Nm, Nfrho);
                Spec=cell(1,2);
                if any(validObs(:))
                    switch c.Value.model
                        case {'rw','sw'}
                            [L, Tc]= AstroUtil.supernova.sn_cooling_sw(c.Value.tobs(iObs,it0),'Type',c.Value.Prog,...
                                'Rs',c.Value.VecRs,'Vs',shiftdim(c.Value.VecVs,-1),'Ms',shiftdim(c.Value.VecMs,-2),...
                                'f_rho',shiftdim(c.Value.VecFrho,-3), 'redshift', c.Value.redshift, 'Model', c.Value.model);
                        case {'msw'}
                            [L, Tc]= AstroUtil.supernova.sn_cooling_msw(c.Value.tobs(iObs,it0),'Type',c.Value.Prog,...
                                'Rs',c.Value.VecRs,'Vs',shiftdim(c.Value.VecVs,-1),'Ms',shiftdim(c.Value.VecMs,-2),...
                                'f_rho',shiftdim(c.Value.VecFrho,-3), 'redshift', c.Value.redshift, 'Model', c.Value.model);
                    end
                    L = L(validObs(:))./(4.*pi.*(c.Value.dist.*constant.pc).^2);
                    Tc = Tc(validObs(:));

                    tmpLC(validObs) = c.Value.F(L,Tc);
                    tmpLC(tmpLC<0) = 0;
                end
                Flux(iObs,:,:,:,:) = tmpLC;
                if transient
                    NTransient = NTransient+validTransObs;
                end
            end
            FluxPnt = LC.flux;
            FluxErr = LC.fluxerr;
            chi2t   = zeros(Nr,Nv,Nm,Nfrho,'single');
            pointst = zeros(Nr,Nv,Nm,Nfrho,'uint16');

            if any(ismember(LC.Properties.VariableNames,'outliers'))
                outliers = LC.outliers;
                chi2to   = zeros(Nr,Nv,Nm,Nfrho,'single');
                pointsto = zeros(Nr,Nv,Nm,Nfrho,'uint16');
            else
                outliers = false;
            end
            
            if ismember(iband,BGbands)
                BGt = zeros(Nr,Nv,Nm,Nfrho,'single');
            else
                BGt = 0;
            end
            
            parfor iRs = 1:Nr
%             for iRs = 1:Nr
%                 valid = (c.Value.tobs(:,it0)>=t_BO(:,iRs,:,:,:)) & (c.Value.tobs(:,it0)<=t_max(:,iRs,:,:,:));
                valid = (W.tobs(:,it0)>=t_BO(:,iRs,:,:,:)) & (W.tobs(:,it0)<=t_max(:,iRs,:,:,:));
                if ismember(iband,BGbands)
%                     BGi = single(sum(~outliers.*valid.*(FluxPnt-Flux(:,iRs,:,:,:))./FluxErr.^2,1)./sum(~outliers.*valid.*FluxErr.^-2,1));
                    k=6; % number of model parameters to estimate
                    sigma2 = sum(~outliers.*valid./FluxErr.^2,1);                                  %$$\Sum_i 1/\sigma_i^2$$
                    OminE  = sum(~outliers.*valid.*(FluxPnt-Flux(:,iRs,:,:,:))./FluxErr.^2,1);      %$$\Sum_i \frac{O_i - E_i}{sigma_i^2}$$  
                    OminE2 = sum(~outliers.*valid.*(FluxPnt-Flux(:,iRs,:,:,:)).^2./FluxErr.^2,1);   %$$\Sum_i \frac{(O_i-E_i)^2}{sigma_i^2}$$
                    
                    a = -sigma2.^2-sigma2./bg{iband}.sigma.^2;
                    b = 3.*sigma2.*OminE + 2./bg{iband}.sigma.^2.*OminE + bg{iband}.value/bg{iband}.sigma.^2.*sigma2;
                    c = (k-2).*sigma2 - 2.*OminE.^2 - OminE2.*(sigma2 - 1./bg{iband}.sigma.^2) - 2.*bg{iband}.value/bg{iband}.sigma.^2.*OminE;                  
                    d = (2-k).*OminE + OminE2.*(OminE + bg{iband}.value./bg{iband}.sigma.^2);
                    
                    sigma2 = []; OminE = []; OminE2 = [];
                    
                    BGi = cubicZeros(a,b,c,d);
                    a = []; b = []; c = []; d = [];
                else
                    BGi = BG;
                end

                OminE = FluxPnt-BGi-Flux(:,iRs,:,:,:);        

                if ModelError>0
                    ModelError2 = (ModelError.*Flux(:,iRs,:,:,:)).^2;
                    chi2i = single(shiftdim(sum(~outliers.*valid.*OminE.^2./(FluxErr.^2+ModelError2),1),1));
                    chi2o = single(shiftdim(sum(outliers.*valid.*OminE.^2./(FluxErr.^2+ModelError2),1),1));
                else
                    chi2i = single(shiftdim(sum(~outliers.*valid.*(OminE./FluxErr).^2,1),1));
                    chi2o = single(shiftdim(sum(outliers.*valid.*(OminE./FluxErr).^2,1),1));
                end
                chi2t(iRs,:,:,:,:) = chi2i;                
                pointst(iRs,:,:,:,:) = uint16(shiftdim(sum(~outliers.*valid,1),1));
                if length(outliers)>1
                    chi2to(iRs,:,:,:,:) = chi2o;
                    pointsto(iRs,:,:,:,:) = uint16(shiftdim(sum(outliers.*valid,1),1));
                end
                if ismember(iband,BGbands)
                    BGt(:,iRs,:,:,:)=shiftdim(BGi,1);
                else
                    BGt = BG;
                end
            end
            loopfile.chi2 = chi2t;
            loopfile.points = pointst;
            if length(outliers)>1
                loopfile.chi2outliers = chi2to;
                loopfile.pointsOutliers = pointsto;
            end
            if transient
                loopfile.NTransient = NTransient;
            end
            loopfile.BG = BGt;
            save(loopfileName,'-struct','loopfile','-v7.3');
            clear loopfile
        end

        fprintf('\nMerging Ebv files.\n');

        for iEbv=1:NEbv
            loopfile = matfile(sprintf('%s_iEbv%d.mat',sn_name,iEbv));
            chi2(it0,:,:,:,:,iEbv) = loopfile.chi2;
            points(it0,:,:,:,:,iEbv) = loopfile.points;
            if ismember(iband,BGbands)
                BG(it0,:,:,:,:,iEbv) = loopfile.BG;            
            end
            if length(outliers)>1
                chi2outliers(it0,:,:,:,:,iEbv) = loopfile.chi2outliers;
                pointsOutliers(it0,:,:,:,:,iEbv) = loopfile.pointsOutliers;
            end
            clear loopfile
            delete(sprintf('%s_iEbv%d.mat',sn_name,iEbv));        
        end
        if transient
            transientPoints(it0,:,:,:,:)=NTransient;
        end
    end
    gridfile.bandschi2(1,iband) = mat2cell(chi2,Nt0);
    gridfile.bandspnts(1,iband) = mat2cell(points,Nt0);
    
    if length(outliers)>1
        gridfile.bandsOutliersChi2(1,iband) = mat2cell(chi2outliers,Nt0);
        gridfile.bandsOutliersPnts(1,iband) = mat2cell(pointsOutliers,Nt0);
    end
    gridfile.bg(1,iband) = mat2cell(BG,size(BG,1));
    if transient
        gridfile.NTransient(1,iband) = {transientPoints};
    end
    gridfile.chi2    = gridfile.chi2 + chi2;
    gridfile.points  = gridfile.points + points;   
end
