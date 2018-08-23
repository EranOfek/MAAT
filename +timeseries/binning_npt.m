function [Bin,ColCell]=binning_npt(TS,Npt,ColCell)
% Binning a time series by equal number of sucssesive points
% Package: timeseries
% Description: Binning a time series by equal number of sucssesive points
% Input  : - A two column [Time property] series.
%          - Number of points in each bin. Default is 10.
%          - Cell array of columns to calculate. Options:
%            'meanbin'   - mean time of bin.
%            'medianbin' - median time of bin.
%            'stdbin'    - std time of bin.
%            'rangebin'  - range time of bin.
%            'minbin'    - min time of bin.
%            'maxbin'    - max time of bin.
%            and any function handle that works on the first dimension.
%            Default is: {'meanbin',@nanmean,@nanstd}.
% Output : - A matrix with the binned data.
%          - Cell array of the requested columns.
%     By : - Eran O. Ofek                  Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [B,C]=timeseries.binning_npt(rand(102,2))
%          [B,C]=timeseries.binning_npt(rand(102,2),30,{'medianbin',@range,@max})
% Reliable : 2

ColTime = 1;
ColProp = 2;

Def.Npt     = 10;
Def.ColCell = {'meanbin',@nanmean,@nanstd};
if (nargin==1)
    Npt      = Def.Npt;
    ColCell  = Def.ColCell;
elseif (nargin==2)
    ColCell  = Def.ColCell;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments: binning_npt(TS,[Npt],[ColCell])');
end

TS = sortrows(TS,ColTime);
% number of points in series
N  = size(TS,ColTime);

% number of bins - remove extra points
Nbin = floor(N./Npt);
TS   = TS(1:Nbin.*Npt,:);

ReshapeTSt = reshape(TS(:,ColTime),Npt,Nbin);
ReshapeTSp = reshape(TS(:,ColProp),Npt,Nbin);

Ncol = numel(ColCell);
Bin  = zeros(Nbin,Ncol);
for Icol=1:1:Ncol
    % for each requested column
    if (ischar(ColCell{Icol}))
        switch lower(ColCell{Icol})
            case 'meanbin'
	        Bin(:,Icol) = nanmean(ReshapeTSt).';
            case 'medianbin'
                Bin(:,Icol) = nanmedian(ReshapeTSt).';
            case 'stdbin'
                Bin(:,Icol) = nanstd(ReshapeTSt).';
            case 'rangebin'
                Bin(:,Icol) = range(ReshapeTSt).';
            case 'minbin'
                Bin(:,Icol) = min(ReshapeTSt).';
            case 'maxbin'
                Bin(:,Icol) = max(ReshapeTSt).';
            otherwise
	        error('Unknown ColCell string option');
        end
     else
         % assume a function on property column
         Bin(:,Icol) = ColCell{Icol}(ReshapeTSp).';
     end
end
