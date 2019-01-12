function [Mag,Flag,FiltEffWave]=synphot(Spec,FiltFam,FiltName,MagSys,Algo,Ebv,R)
% Synthetic photometry of spectrum
% Package: AstroUtil.spec
% Description: Calculate synthetic photometry of a spectrum.
% Input  : - Spectrum [wavelength(Ang), Flux(F_{\lambda})].
%          - Filter normalized transmission curve,
%            or a string containing filter familiy name,
%            or AstFilter class object.
%            See AstFilter.get for details.
%            Filter transmission curve override filter name.
%          - Filter name, see AstFilter.get.m for details.
%          - Magnitude system: {'Vega' | 'AB'}
%          - Algorithm used:
%            'Poz' - D. Poznanski, basic_synthetic_photometry.m
%            'cos' - transmission curve interpolated on spectrum.
%                    Default.
%            'soc' - spectrum interpolated on transmission curve.
%            If empty matrix, then use default.
%          - E_{B-V} extinction to apply to spectrum before calculating
%            the synthetic photometry. Default is 0.
%            This function is using the Cardelli et al. model implemented
%            in extinction.m
%          - R_v of extinction. Default is 3.08.
% Output : - Synthetic magnitude.
%          - The fraction of flux that was extrapolated in case of
%            partial coverage between spectrum and filter.
%            0 - means no extrapolation.
%          - Filter effective wavelength [Ang].
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    May 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Comments: The 'Poz' option requires basic_synthetic_photometry.m 
%           by: Dovi Poznanski.
% Example: [Mag,Flag]=synphot(Spec,'SDSS','r','AB');
%          [Mag,Flag]=synphot(Spec,'SDSS','r','AB',[],0.1); % apply extinction
% Reliable: 1
%------------------------------------------------------------------------------
InterpMethod = 'linear';
Def.Algo = 'cos';
Def.Ebv  = 0;
Def.R    = 3.08;
if (nargin==4)
   Algo   = Def.Algo;
   Ebv    = Def.Ebv;
   R      = Def.R;
elseif (nargin==5)
   Ebv    = Def.Ebv;
   R      = Def.R;
elseif (nargin==6)
   R      = Def.R;
elseif (nargin==7)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (isempty(Algo))
   Algo = Def.Algo;
end

if (ischar(FiltFam)==1)
   %Filter = get_filter(FiltFam,FiltName);
   Filter = AstFilter.get(FiltFam,FiltName);
   Tran   = Filter(1).nT;
   FiltEffWave = Filter(1).eff_wl;
elseif (AstFilter.isAstFilter(FiltFam))
    Tran = FiltFam(1).nT;
    FiltEffWave = FiltFam(1).eff_wl;
else
   Tran   = FiltFam;
   FiltEffWave = trapz(Tran(:,1),Tran(:,1).*Tran(:,2))./trapz(Tran(:,1),Tran(:,2));
end
%FiltEffWave = Filter.eff_wl{1};

if (iscell(Spec))
    if ~isvector(Spec{1})
        error('MAAT:AstroUtils:spec:synphot:spec_wave_no_vector',...
            'Error. The first component of the Spec cell must be a vector.');
    else
        if ismatrix(Spec{2})
            if size(Spec{2},1)~=length(Spec{1})
                if size(Spec{2},2)==length(Spec{1})
                    Spec{2} = Spec{2}.';
                else
                    error('MAAT:AstroUtils:spec:synphot:spec_wrong_elements',...
                        'Error. One of the two dimensions of the second component of the Spec cell must be in agreement with the number of elements of the first component.');
                end
            end
        elseif isvector(Spec{2})
            if length(Spec{2})~=length(Spec{1})
                error('MAAT:AstroUtils:spec:synphot:spec_wrong_elements',...
                    'Error. The number of elements in the two components of the Spec cell must be in agreement with eachother.');
            end
        else
            error('MAAT:AstroUtils:spec:synphot:spec_wrong_elements',...
                'Error. Only matrix and vector are supported in the second component of the Spec cell.');
        end
    end
end             
                  
if (Ebv>0)
   % apply extinction
   A = extinction(Ebv,Spec(:,1)./10000,[],R);
   Spec(:,2) = Spec(:,2).*10.^(-0.4.*A);
end

TranNorm = trapz(Tran(:,1),Tran(:,2));

switch lower(Algo)
 case 'poz'
    [Mag,Dmag,Flag] = basic_synthetic_photometry(Spec,Tran,MagSys,[],[0 1]);
 otherwise

    Direction = 'cos';
    switch lower(Direction)
     case {'curve_on_spec','cos'}
        % Interp transminssion curve on Spec
        if iscell(Spec)
            [Spec,Tran]     = AstroUtil.spec.eq_sampling(Spec,Tran,Spec{1},InterpMethod);
        else
            [Spec,Tran]     = AstroUtil.spec.eq_sampling(Spec,Tran,Spec(:,1),InterpMethod);
%             I = find(~isnan(Tran(:,2)));
%             Spec = Spec(I,:);
%             Tran = Tran(I,:);
        end
     case {'spec_on_curve','soc'}
        % Interp Spec on transminssion curve
        Spec     = AstroUtil.spec.eq_sampling(Spec,Tran,Tran(:,1));
     otherwise
        error('Unknown Direction option');
    end
    
    if iscell(Spec)
        [~,MinI] = min(Spec{1});
        [~,MaxI] = max(Spec{1});
    else
        [~,MinI] = min(Spec(:,1));
        [~,MaxI] = max(Spec(:,1));
    end
    if (nargout>=2)
        if (MaxI-MinI)<2
            Flag = 1;
        else
            Flag = 1 - trapz(Tran(MinI:MaxI,1),Tran(MinI:MaxI,2))./TranNorm;
        end
        if (Flag<1e-5)
            Flag = 0;
        end
    end
    %Min = min(min(Tran(:,1))-min(Spec(:,1)),0);
    %Max = max(max(Spec(:,1))-max(Tran(:,1)),0);
    %Flag = (Min+Max)./range(Tran(:,1));

    switch lower(MagSys)
     case 'ab'
        % convert F_lambda to F_nu
        Freq     = convert.energy('A','Hz',Tran(:,1));
        if iscell(Spec)
            SpecFnu  = convert.flux(Spec{2},'cgs/A','cgs/Hz',Tran(:,1),'A');
        else
            SpecFnu  = convert.flux(Spec(:,2),'cgs/A','cgs/Hz',Tran(:,1),'A');
        end

        if (size(Tran,1)<2)
 %           NormTran = 1;
            Fnu      = NaN;
        else
            NormTran = trapz(Freq,Tran(:,2));
            Fnu      = trapz(Freq,SpecFnu.*Tran(:,2))./NormTran;
            %Flam     = convert_flux(Fnu,'cgs/Hz','cgs/A',FiltEffWave,'A');
        end
        Mag      = -48.6 - 2.5.*log10(Fnu);
        
     case 'vega'
        load vega_spec.mat;
        VegaF   = AstroUtil.spec.eq_sampling(vega_spec,Tran,Tran(:,1));
%        Freq    = convert.energy('A','Hz',Tran(:,1));
        Fvega   = trapz(Tran(:,1),Spec(:,2).*Tran(:,2)./VegaF(:,2));
        Mag     = -2.5.*log10(Fvega);

     otherwise
        error('Unknown MagSys option');
    end
end
