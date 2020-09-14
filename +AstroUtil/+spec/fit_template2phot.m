function [Res,SynMag,ExtrapFlag]=fit_template2phot(PhotData,varargin)
% Fit a set of spectral templates to photometric observations of a source.
% Package: AstroUtil
% Description: Fit a set of spectral templates to photometric observations
%              of a single source. Return the spectral template that best
%              fit the photometric data.
% Input  : - Photometric data of a single source to which to fit spectral
%            templates.
%            This is one of the following:
%            Cell array of {Mag, Err, BandFamily, BandName, MagType}
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Template' - An AstSpec object of spectral templates to fit.
%                       Default is AstSpec.get_pickles([],'V').
%            'SelectionMethod' - Method by which to select best fit
%                       template. {'rms' | 'chi2'}. Default is 'rms'.
%            'IgnorePartial' - A flag indicating if to ignore synthetic
%                       magnitude which is based on extrapolation.
%                       Default is true.
%            'SynMag' - A matrix of synthetic magnitude.
%                       Default is empty.
%                       This is useful in order to avoid re-calculating the
%                       synthetic magnideu multiple times.
%            'ExtrapFlag' - Like SynMag but for the extrapolation flag.
% Output : - A structure with the following fields:
%            'BestRMS' - RMS for best template.
%            'BestChi2' - Chi^2 for best template.
%            'NfiltUse' - Number of filter used for fitting each template.
%            'BestTempInd' - Index of best fitted template.
%          - Matrix of synthetic magnitude [template, filter]
%          - Matrix of syn.mag interpolation flag (0 if ok)
%            [template, filter].
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.Template             = AstSpec.get_pickles([],'V');
DefV.SelectionMethod      = 'rms';  % 'rms' | 'chi2'
DefV.IgnorePartial        = true; % ignore interpolated syn phot
DefV.SynMag               = [];
DefV.ExtrapFlag           = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (iscell(PhotData))
    Mag     = [PhotData{:,1}]';
    Err     = [PhotData{:,2}]';
    Family  = PhotData(:,3);
    Band    = PhotData(:,4);
    MagType = PhotData(:,5);
else
    error('Unknown PhotData type');
end

% number of template to fit
Ntemplate = numel(InPar.Template);
% number of phot data
Nphot     = numel(Mag);


if isempty(InPar.SynMag)
    SynMag     = nan(Ntemplate,Nphot);
    ExtrapFlag = nan(Ntemplate,Nphot);
    for Iphot=1:1:Nphot
        [SynMag(:,Iphot), ExtrapFlag(:,Iphot)] = synphot(InPar.Template,Family{Iphot},Band{Iphot},MagType{Iphot});
    end
else
    SynMag = InPar.SynMag;
    ExtrapFlag = InPar.ExtrapFlag;
    if (isempty(ExtrapFlag))
        ExtrapFlag = ones(size(SynMag));
    end
end
    
Diff = (Mag.' - SynMag);
if (InPar.IgnorePartial)
    FlagOK = ExtrapFlag==0;
    
    Diff(~FlagOK) = NaN;
    
else
    % use interpolated values
end
Chi2 = nansum((Diff./Err.').^2,2);
RMS  = std(Diff,[],2);
NfiltUse = sum(ExtrapFlag==0,2);  % number of filters fitted to each template

% select best match
switch lower(InPar.SelectionMethod)
    case 'rms'
        [BestRMS,MinInd] = min(RMS);
        BestChi2 = Chi2(MinInd);
    case 'chi2'
        [BestChi2,MinInd] = min(Chi2);
        BestRMS = RMS(MinInd);
    otherwise
        error('Unknown SelectionMethod option');
end


Res.BestRMS     = BestRMS;
Res.BestChi2    = BestChi2;
Res.NfiltUse    = NfiltUse;
Res.BestTempInd = MinInd;








    