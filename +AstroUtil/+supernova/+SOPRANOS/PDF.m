function [PDF, chi2, dof, chi2bands, pntsbands, transbands] = ...
    PDF(data, redshift, Vect0, Rs, Vs, VecMs, VecFrho, VecEbv, model, n_params, ProgType, priors, bg, transient, mintpoints, c)
% calculate Probabiliy Distribution Function (PDF) for a given SW17 shock cooling model
% Package: AstroUtil.supernove.SOPRANOS
% Description: calculate Probabiliy Distribution Function (PDF) for a given SW17 shock cooling model
% Input  : - cell array of bands, each on containts the observations and
%               AstFilterObj
%          - redshift
%          - Vector of refernce time values.
%          - Progenitor radius [R_\sun]
%          - Shock velocity paraneter [cm s^-1]
%          - Vector of ejecta masses values [M_\sun]
%          - Vector of envelope density parameter values
%          - Vector of extinction values
%          - Models type [SW or RW]
%          - Number of parameters estimated (for calculation of likelihood)
%          - Prgoenitor type [RSG or BSG]
%          - Parameter priors range
%               (Model outside the prior will result 0 likelihood)
%          - Background values for each one of the bands for sustraction.
%          - The transient range for transinet points conatraint
%          - The number of minimal number of transient points for each band
%               (0 means no constraint)
%          - parallel.pool.Constant which includes the data input to reduce
%            communication between workers.
% Output : - Likelihood value for the given model. In case Vect0, VecMs,
%               VecFrho or VecEbv includes more than 1 value they are
%               marginalized out resulting a scalar likelihood.
%          - A matrix spanned by Vect0, VecMs, VecFrho and VecEbv contains
%               the chi2 sum for each model 
%          - A matrix spanned by Vect0, VecMs, VecFrho and VecEbv contains
%               the degrees of freeddom for each model 
%          - A cell array with matrices spanned by Vect0, VecMs, VecFrho
%               and VecEbv contains the chi2 sum for each band & each model 
%          - A cell array with matrices spanned by Vect0, VecMs, VecFrho
%               and VecEbv contains the number of valid data points for
%               each band & each model.
%          - A cell array with matrices spanned by Vect0, VecMs, VecFrho
%               and VecEbv containts the number of valid data 
%               points for each band & each model 
%          - A cell array with matrices spanned by Vect0, VecMs, VecFrho
%               and VecEbv containts the number of valid data 
%               points with the transient for each band & each model 
%           %               
% See also: AstroUtil.supernova.SOPRANOS.findMaximum
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

    if nargin<16
        c=parallel.pool.Constant(data);
    end
    
    if ( Rs<priors.Rs(1)||Rs>priors.Rs(end)||Vs<priors.Vs(1)||Vs>priors.Vs(end)||...
            Vect0(1)<priors.t0(1)||Vect0(end)>priors.t0(end)||VecMs(1)<priors.Ms(1)||VecMs(end)>priors.Ms(end)||...
            VecFrho(1)<priors.frho(1)||VecFrho(end)>priors.frho(end)||VecEbv(1)<priors.Ebv(1)||VecEbv(end)>priors.Ebv(end) )
        PDF  = 0;
        chi2 = inf(size(Vect0));
        dof  = zeros(size(Vect0));
    else       
        
        chi2bands = zeros(length(data),length(Vect0),length(VecMs),length(VecFrho),length(VecEbv));
        pntsbands = zeros(length(data),length(Vect0),length(VecMs),length(VecFrho));
        transbands = zeros(length(data),length(Vect0),length(VecMs),length(VecFrho));
        
        for iband = 1:length(data)
            nMs   = length(VecMs);
            nFrho = length(VecFrho);
            chi2 = zeros(length(Vect0),nMs*nFrho,length(VecEbv));
            pnts = zeros(length(Vect0),nMs*nFrho);
            transpnts = zeros(length(Vect0),nMs*nFrho);
            parfor i = 1:nMs*nFrho
                [iMs, iFrho] = ind2sub([nMs nFrho], i);
                LC = sortrows(unique(c.Value{iband}.LC),'MJD');
                if any(contains(LC.Properties.VariableNames,'outliers'))
                    LC = LC(~LC.outliers,:);
                end
                [t, t0] = ndgrid(LC.MJD,Vect0);
                trans = (t>=transient(1)&t<=transient(2));
                t = t-t0;

                switch model
                    case {'sw','rw'}
                        [~, ~, ~, ~, ~, t_min, t_max, t_opac,~,t_delta] = ...
                            AstroUtil.supernova.sn_cooling_sw(1, 'Type', ProgType, 'Vs', Vs, 'Rs', Rs, ...
                            'Ms', VecMs(iMs), 'f_rho',VecFrho(iFrho), 'redshift', redshift, 'model',model);
                        valid = [];
                        if strcmpi(model,'sw')
                            valid = (t >= t_min.*(1+redshift)) & (t <= min(t_max,t_opac).*(1+redshift));
                        elseif strcmpi(model,'rw')
                            valid = (t >= t_min.*(1+redshift)) & (t <= min(t_delta,t_opac).*(1+redshift));
                        else
                            error('unknown model %s',model);
                        end
                    case 'msw'
                        [~, ~, ~, ~, ~, t_min, t_max] = ...
                            AstroUtil.supernova.sn_cooling_msw(1, 'Type', ProgType, 'Vs', Vs, 'Rs', Rs, ...
                            'Ms', VecMs(iMs), 'f_rho',VecFrho(iFrho), 'redshift', redshift, 'model',model);
                        valid = (t>=t_min) & (t<=t_max);
                end
                
               if any(valid(:))
                   switch model
                       case {'sw','rw'}
                           [~,~,~,~, Mag] = AstroUtil.supernova.sn_cooling_sw(t(valid), 'Type', ProgType, 'Vs', Vs, 'Rs', Rs, ...
                                           'Ms', VecMs(iMs),  'f_rho', VecFrho(iFrho), 'Ebv', shiftdim(VecEbv,-2),'redshift', redshift, ...
                                           'FiltFam', c.Value{iband}.filterObj,'model',model);
                       case 'msw'
                           [~,~,~,~, Mag] = AstroUtil.supernova.sn_cooling_msw(t(valid), 'Type', ProgType, 'Vs', Vs, 'Rs', Rs, ...
                                           'Ms', VecMs(iMs),  'f_rho', VecFrho(iFrho), 'Ebv', shiftdim(VecEbv,-2),'redshift', redshift, ...
                                           'FiltFam', c.Value{iband}.filterObj,'model',model);
                   end
                   F = zeros([size(t), length(VecEbv)]);
                   F(valid&shiftdim(VecEbv==VecEbv,-2))=convert.flux(Mag,'AB','cgs/A', c.Value{iband}.filterObj.pivot_wl_photon, 'A');
                   chi2(:,i,:) = squeeze(sum(valid.*((LC.flux - bg(iband) - F)./LC.fluxerr).^2,1));
                   pnts(:,i) = squeeze(sum(valid,1));
                   transpnts(:,i) = squeeze(sum(valid&trans,1));
               else
                   chi2(:,i,:)=zeros(length(Vect0),1,length(VecEbv));
                   pnts(:,i)=zeros(length(Vect0),1);
                   transpnts(:,i) = zeros(length(Vect0),1);                   
               end
            end           
            
            chi2bands(iband,:,:,:,:) = reshape(chi2,length(Vect0),nMs,nFrho,length(VecEbv));
            pntsbands(iband,:,:,:) = reshape(pnts,length(Vect0),nMs,nFrho);
            transbands(iband,:,:,:) = reshape(transpnts,length(Vect0),nMs,nFrho);
                        
        end
        
        % check which data grids has less than the minimum required number
        % of data points:
        valid = (transbands>=mintpoints(:))&shiftdim(true(length(VecEbv),1),-4);
        valid = all(valid,1);

        chi2 = shiftdim(sum(chi2bands,1),1);
        dof  = shiftdim(sum(pntsbands,1),1)-n_params;
        dof  = dof.*(shiftdim(ones(length(VecEbv),1),-3));
        
        chi2(~valid) = 0;
        
        if any(any(any(any(valid,2),3),4),5)
            PDF = chi2pdf(chi2,dof);
            PDF(~valid) = 0;
            PDF(dof<=0) = 0;        % models with less data points than number of params - dof is defined as valid data points - n_params
            if length(VecEbv)>1
                PDF = trapz(VecEbv,PDF,4);
            end           
            if length(VecFrho)>1
                PDF     = trapz(log10(VecFrho),PDF,3);
            end
            if length(VecMs)>1
                PDF     = trapz(VecMs,PDF,2);
            end
            if length(Vect0)>1
                PDF     = trapz(Vect0,PDF,1);
            end
        else
            % there are no valid phase space points
            PDF  = 0;
        end
    end    
end
