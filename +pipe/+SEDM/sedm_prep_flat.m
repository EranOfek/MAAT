function [Res,SI]=sedm_prep_flat(varargin)
%--------------------------------------------------------------------------
% sedm_prep_flat function                                             SEDM
% Description: 
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


DefV.FlatKey       = 'Calib: twilight flat';
DefV.SI            = [];
DefV.SIflat        = [];   % SIflat is the Flat SI after efine_slopepos, extract_spexcell and copy_wavecalib
DefV.SIarc         = [];
DefV.FlatImageName = [];
DefV.SpecField     = 'SpexSpecFit';
DefV.WaveRange     = [450 850];
DefV.WaveGrid      = logspace(log10(360),log10(950),500).';
DefV.InterpMethod  = 'linear';
DefV.MeanFun       = @nanmedian;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

SI     = InPar.SI;
SIarc  = InPar.SIarc;
SIflat = InPar.SIflat;

if (~isempty(SIflat)),
    % SIflat is provided
    % continue
else
    if (~isempty(SI) && ~isempty(SIarc)),
        % SI and SIarc arc provided
        
        if (isempty(InPar.FlatImageName)),
            error('FlatImageName must be provided');
        end

        SIflat = sedm_refine_slopepos('Input',FlatImageName,'SI',SI);
        SIflat = sedm_extract_spexcell('SI',SIflat,'ScienceImage',FlatImageName,'SubSpexcellBack',...
                                       'SpexFit','SpexFitMF',20); 
        SIflat = sedm_copy_wavecalib(SIarc,SIflat);
    else
        error('either SIflat or SI and SIarc must be provided');
    end
end
    
Nseg = numel(SIflat);
Res.SpecMat = zeros(Nseg,length(InPar.WaveGrid));

for Iseg=1:1:Nseg,
    
    if (isempty(SIflat(Iseg).WaveCalib)),
        SpecMat(Iseg,:) = NaN;
    else
        FlatSpecSeg = [SIflat(Iseg).WaveCalib.', SIflat(Iseg).(InPar.SpecField).'];
        
        NN = ~isnan(FlatSpecSeg(:,1));
        FW = interp1(FlatSpecSeg(NN,1),FlatSpecSeg(NN,2),...
                     InPar.WaveGrid,InPar.InterpMethod);

        Res.SpecMat(Iseg,:) = FW.';
    end
end

Res.WaveGrid = InPar.WaveGrid;
Res.Xpos     = [SIflat.MeanX].';
Res.Ypos     = [SIflat.MeanY].';


FlagWave = InPar.WaveGrid>InPar.WaveRange(1) & InPar.WaveGrid<InPar.WaveRange(2);
%Res.SumFlat  = sum(Res.SpecMat(:,FlagWave),2);
%Res.NormSpecFlat  = bsxfun(@times,Res.SpecMat(:,FlagWave), 1./InPar.MeanFun(Res.SpecMat(:,FlagWave),1));
Res.NormFlat  = sum(Res.SpecMat(:,FlagWave),2);
Res.NormFlat  = Res.NormFlat./InPar.MeanFun(Res.NormFlat);

