function Sim=populate_sim(Sim,varargin)
% Populate SIM object with background, error, talog and PSF
% Package: @SIM
% Description: Given a SIM object with an image, populate all the other
%              fields (BackIm, ErrIm, Cat, PSF).
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'RePop'  - Repopulate field even if exist. Default is false.
%            'BackPar'- Cell array of parameters to pass to background.m
%                       for the background and std image calculation.
%                       Default is {}.
%            'CatFun' - Function to call for generating the catalog.
%                       Default is @mextractor.m
%            'CatFunPar' - Cell array of parameters to pass to 'CatFun'
%                       function.
%                       Default is {}.
%            'PSFPar' - Cell array of parameters to pass to
%                       psf_estimator.m. Default is {}.
% Output : - A SIM object with the BackIm, ErrIm, Cat, PSF fields
%            populated.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=populate_sim(S);
% Reliable: 
%--------------------------------------------------------------------------

BackField    = 'BackIm';
ErrField     = 'ErrIm';
CatField     = 'Cat';
PsfField     = 'PSF';

DefV.RePop              = false;
DefV.BackPar            = {};
DefV.CatFun             = @mextractor;
DefV.CatFunPar          = {};
DefV.PSFPar             = {};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);


% Check what to populate
if (InPar.RePop)
    PopBack = true(Nsim,1);
    PopCat  = true(Nsim,1);
    PopPSF  = true(Nsim,1);
else
    PopBack = ~(isfield_populated(Sim,BackField) && isfield_populated(Sim,ErrField));
    PopCat  = ~(isfield_populated(Sim,CatField) && isfield_populated(Sim,CatField));
    PopPSF  = ~(isfield_populated(Sim,PsfField) && isfield_populated(Sim,PsfField));
end

% populate background and error images
if (any(PopBack))
    Sim(PopBack) = background(Sim(PopBack),InPar.BackPar);
end

% populate catalog
if (any(PopCat))
    Sim(PopCat) = InPar.CatFun(Sim(PopCat),InPar.CatFunPar{:});
end

% populate PSF
if (any(PopPSF))
    Sim(PopPSF) = psf_estimator(Sim(PopPSF),InPar.PSFPar{:},'OutType','origsim');
end

warning('NEED to implement: WCS, 2nd iter catalog')

