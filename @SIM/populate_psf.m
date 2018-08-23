function Sim=populate_psf(Sim,varargin)
% Populate the PSF SIM object
% Package: class/@SIM
% Description: Populate the PSF SIM object
% Input  : - A SIM object
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FilterFindThresh' - Threshold, in units of sigmas, for
%                           finding bright stars for PSF extraction.
%                           Default is 10.
%            'DefFilter'    - Default filter to use for PSF stars search.
%                           Default is @Kernel2.gauss;
%            'DefFilterFunPar' - Default filters parameters.
%                           Default is {1.5,1.5,0,15,15}.
%            'RePopulate' - Repopulate PSF. Default is false.
%            'GetPsfPar'  - Cell array of additional arguments to pass to
%                           ClassPSF/getmpsf.m.
%            'PSF_EstimatorPar' - Cell array of additional arguments to
%                           pass to @psf_cat_selector.m.
%            'PSF_SelectPar' - Cell array of additional arguments to pass
%                           to @psf_estimator.m.
% Output : - A SIM object with populated PSF.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=populate_psf(Sim)
% Reliable: 2
%--------------------------------------------------------------------------

DefV.FilterFindThresh   = 10;
DefV.DefFilter          = @Kernel2.gauss;
DefV.DefFilterFunPar    = {1.5,1.5,0,15,15};
DefV.RePopulate         = false;
DefV.GetPsfPar          = {};
% DefV.PSF_EstimatorPar   = {'NoiseFactor',3,'Rad0',6}; %120,'StampSize',120};
% DefV.PSF_SelectPar      = {'ColSN','SN',...
%                            'ColMaxFlux','PEAKF_VALTOT',...
%                            'ColBitFlag','FLAGS',...
%                            'MaskVal',0,...
%                            'PosCol',{'XWIN_IMAGE','YWIN_IMAGE'},...
%                            'BoundryDist',10,...
%                            'MinSN',10,...
%                            'SatLevel',40000};
                       
DefV.HalfSizePSF          = 10;
DefV.PSFSelectorFun       = @psf_cat_selector;
DefV.PSFSelectorPar       = {'ColSN','SN',...
                           'ColMaxFlux','PEAKF_VALTOT',...
                           'ColBitFlag','FLAGS',...
                           'MaskVal',0,...
                           'PosCol',{'XWIN_IMAGE','YWIN_IMAGE'},...
                           'BoundryDist',10,...
                           'MinSN',15,...
                           'SatLevel',40000};
DefV.PSFCombinerPar       = {};                       
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Nsim    = numel(Sim);


if (InPar.RePopulate)
    FlagPop = true(Nsim,1);
else
    Psf     = getmpsf(Sim,InPar.GetPsfPar{:});
    FlagPop = Util.cell.isempty_cell(Psf);
end

SimPass1          = Sim;
SimPass1(FlagPop) = mextractor(Sim(FlagPop),'FilterFind',false,...
                           'Thresh',InPar.FilterFindThresh,...
                           'ColCell',{'XWIN_IMAGE','YWIN_IMAGE','SN','SN_UNF',...
                                      'FLAGS','PEAKF_VALTOT','X2WIN_IMAGE','Y2WIN_IMAGE','XYWIN_IMAGE',...
                                      'NEAREST_SRCDIST'},...
                           'AddFilter',{},...
                           'CleanByBackGrad',false,...
                           'DefFilter',InPar.DefFilter,...
                           'DefFilterFunPar',InPar.DefFilterFunPar);
  
%-----------------------------
%--- Execute psf_estimator ---
%-----------------------------
% Note that 'SelectFunPar' should be a cell array and hence
% the use of (:)
% [SimPass1(FlagPop),MetaData.ResPSF] = psf_estimator(SimPass1(FlagPop),'OutType','OrigSim',...
%                                       'SelectFun',@psf_cat_selector,...
%                                       InPar.PSF_EstimatorPar{:},...
%                                       'SelectFunPar',InPar.PSF_SelectPar(:));

[SimPass1(FlagPop),~]=psf_extractor(SimPass1(FlagPop),'StampHalfSize',InPar.HalfSizePSF,...
                               'PSFSelectorFun',InPar.PSFSelectorFun,...
                               'PSFSelectorPar',InPar.PSFSelectorPar,...
                               'PSFCombinerPar',InPar.PSFCombinerPar);

                                  
% Copy PSF from SimPass1 to Sim
Psf     = getmpsf(SimPass1,InPar.GetPsfPar{:});
Sim     = insert_psf(Sim,Psf);
     
