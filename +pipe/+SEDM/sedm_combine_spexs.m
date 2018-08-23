function [SpecMat,WaveGrid]=sedm_combine_spexs(SI,varargin)
%--------------------------------------------------------------------------
% sedm_combine_spexs function                                         SEDM
% Description: Given a list of spexcell indices, interpolate all the
%              spectra to a common grid and coadd.
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

'need to finish this function...'

DefV.Ind          = [];
DefV.WaveGrid     = logspace(log10(360),log10(950),500).';
DefV.FluxField    = 'SpexSpecFit';
DefV.InterpMethod = 'linear';
DefV.GoodRangeX                = [150 1900];
DefV.GoodRangeY                = [100 1950];

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

Nseg = numel(SI);
if (isempty(InPar.Ind)),
    InPar.Ind = (1:1:Nseg).';
end
    
Nind = length(InPar.Ind);
SpecMat = zeros(Nind,length(InPar.WaveGrid));
for Iind=1:1:Nind,
    
    if (length(SI(InPar.Ind(Iind)).WaveCalib)>2 && SI(InPar.Ind(Iind)).MeanX>InPar.GoodRangeX(1) && ...
                                                   SI(InPar.Ind(Iind)).MeanX<InPar.GoodRangeX(2) && ...
                                                   SI(InPar.Ind(Iind)).MeanY>InPar.GoodRangeY(1) && ...
                                                   SI(InPar.Ind(Iind)).MeanY<InPar.GoodRangeY(2)),
                
        NN = ~isnan(SI(InPar.Ind(Iind)).WaveCalib);
        FW = interp1(SI(InPar.Ind(Iind)).WaveCalib(NN),SI(InPar.Ind(Iind)).(InPar.FluxField)(NN),...
                     InPar.WaveGrid,InPar.InterpMethod);
    else
        FW = zeros(length(InPar.WaveGrid),1).*NaN;
    end
    SpecMat(Iind,:) = FW.';
end

WaveGrid = InPar.WaveGrid;             
             