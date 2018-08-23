function [Sim,SimPeaks]=local_maxima(Sim,ColXY,ExecField,RegionMaxConn)
%--------------------------------------------------------------------------
% local_maxima function                                         class/@SIM
% Description: Use imregionalmax.m to locate all the local maxim in a SIM
%              image and produce an AstCat object of these local maxima.
%              Note that in a noisy image, all the local maxima
%              corresponding to noise will be found too.
% Input  : - A SIM object.
%          - A cell array of the X and Y position columns in which to write
%            the X and Y peak positions in the AstCat object.
%            Default is {'XPIX_PEAK','YPIX_PEAK'}.
%          - A string containing the field name in the SIM
%            object on which to perform the thresholding.
%            Default is 'Im'.
%          - Connectivity for imregionalmax. Default is 8.
% Output : - The input SIM object with an AstCat object indicating the X
%            and Y positions of the pixel (integer values) in which local
%            maxima were found.
%          - A SIM object with the PeaksIm.
% See also: threshold.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=local_maxima(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField     = 'Im';
CatField       = 'Cat';
ColCellField   = 'ColCell';

Def.ColXY              = {'XPIX_PEAK','YPIX_PEAK'};
Def.ExecField          = ImageField;
Def.RegionMaxConn      = 8;  % [4 | 8]

if (nargin==1)
    ColXY         = Def.ColXY;
    ExecField     = Def.ExecField;
    RegionMaxConn = Def.RegionMaxConn;
elseif (nargin==2)
    ExecField     = Def.ExecField;
    RegionMaxConn = Def.RegionMaxConn;
elseif (nargin==3)
    RegionMaxConn = Def.RegionMaxConn;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments: local_maxima(Sim,ColXY,ExecField,RegionMaxConn)');
end


% allocate output
if (nargout>1)
    SimPeaks = SIM(size(Sim));
end

Nsim = numel(Sim);
for Isim=1:1:Nsim
    % locate peaks
    
    if (sum(Sim(Isim).(ExecField))==0)
        PeaksIm = Sim(Isim).(ExecField);
    else
        PeaksIm  = imregionalmax(Sim(Isim).(ExecField),RegionMaxConn);
    end

    %PeaksInd = find(PeaksIm>0);
    PeaksInd = find(PeaksIm);
    % Read peaks coordinates
    [PeakY,PeakX]    = ind2sub(size(PeaksIm),PeaksInd);

    Sim(Isim).(CatField)     = [PeakX, PeakY];
    Sim(Isim).(ColCellField) = ColXY;
    Sim(Isim)                = colcell2col(Sim(Isim));
    
    if (nargout>1)
        SimPeaks(Isim).(ExecField) = PeaksIm;
    end
end
