function [Res,OutCol]=binning(Data,BinSize,StartEnd,OutCol,OutType)
% Binning function. Calculate various functions in data bins.
% Package: timeseries
% Description: Binning a timeseries.
%              In each bin, calculate various functions (e.g., mean) or
%              any user supplied functions.
%              Note that this function was changed on June 2016, and it
%              is not backward compatible. For the old binning.m function
%              use binning_old.m.
%              Note that NaN values are removed prior to binning.
% Input  : - Matrix, in which the first column is the depandent variable
%            by which to bin (e.g., time) and the second column is for the
%            observations. An optional third column may contains the errors
%            in the observations (default is equal weight errors of 1).
%          - Bin size, in units of the "time" column.
%            Alternaively this can be a vector of edges (override the
%            StartEnd input argument).
%            Default is to divide the "time" range into 10 bins.
%          - A two element vector of [Start End] "time".
%            A NaN value for either Start or End will be replaced with the
%            min(time) and max(time), respectively.
%            Default is [NaN NaN].
%          - Cell array of column names to calculate. Possible column
%            names are:
%            'MidBin'
%            'StartBin'
%            'EndBin'
%            'MeanBin'
%            'MedianBin'
%            'StdBin'
%            Or any function that returns a scalar to apply to the
%            observations (e.g., @numel, @mean, @median).
%            There are several special functions that will be applied also
%            for the errors column. These are @wmean, @wmedian.
%            Default is {'MidBin', @numel, @mean, @median, @std}.
%          - Output type: 'mat'|'astcat'. Default is 'mat'.
% Output : - Matrix or an AstCat object with with the requested columns.
%          - Output column names.
% See also: bin_by_eye.m, runmean.m, bin_en.m
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: B=timeseries.binning(rand(1000,2),0.1,[0 1]);
%          B=timeseries.binning(rand(1000,2),0.1,[0 1],{'MeanBin',@rstd,@numel});
% Reliable: 2
%--------------------------------------------------------------------------

Col.X = 1;
Col.Y = 2;
Col.E = 3;

% Clean NaNs
Flag = all(~isnan(Data),2);
Data = Data(Flag,:);

Def.Nbin    = 10;
Def.OutCol  = {'MidBin', @numel, @mean, @median, @std};
Def.OutType = 'mat';
if (nargin==1)
    StartEnd = [min(Data(:,Col.X))-eps, max(Data(:,Col.X))+eps];
    BinSize  = diff(StartEnd)./Def.Nbin;
    OutCol   = Def.OutCol;
    OutType  = Def.OutType;
elseif (nargin==2)
    StartEnd = [min(Data(:,Col.X)), max(Data(:,Col.X))];
    OutCol   = Def.OutCol;
    OutType  = Def.OutType;
elseif (nargin==3)
    OutCol   = Def.OutCol;
    OutType  = Def.OutType;
elseif (nargin==4)
    OutType  = Def.OutType;
elseif (nargin==5)
    % do nothing
else
    error('Illegal number of input arguments: binning(Data,BinSize,StartEnd,OutCol,OutType)');
end
    
   

if (numel(BinSize)>1)
    % BinSize is a vector of edges
    % override StartEnd
    Edges = BinSize;
    Nbin  = numel(Edges)-1;
else
    % BinSize is the bin size
    % deal with NaNs in StartEnd
    if (isnan(StartEnd(1)))
        StartEnd(1) = min(Data(:,Col.X))-eps;
    end
    if (isnan(StartEnd(2)))
        StartEnd(2) = max(Data(:,Col.X))+eps;
    end

    Nbin    = floor(diff(StartEnd)./BinSize);
    BinSize = diff(StartEnd)./Nbin;
    Edges = (StartEnd(1):BinSize:StartEnd(2))';
end

Ncol = numel(OutCol);


Res   = zeros(Nbin,Ncol);
for Ibin=1:1:Nbin
    Flag = Data(:,Col.X)>Edges(Ibin) & Data(:,Col.X)<=Edges(Ibin+1);
    
    % calculate columns
    for Icol=1:1:Ncol
        if (ischar(OutCol{Icol}))
            % Col is a pre-defined string
            switch lower(OutCol{Icol})
                case 'midbin'
                    Res(Ibin,Icol) = (Edges(Ibin) + Edges(Ibin+1)).*0.5;
                case 'startbin'
                    Res(Ibin,Icol) = Edges(Ibin);
                case 'endbin'
                    Res(Ibin,Icol) = Edges(Ibin+1);
                case 'meanbin'
                    Res(Ibin,Icol) = mean(Data(Flag,Col.X));
                case 'medianbin'
                    Res(Ibin,Icol) = median(Data(Flag,Col.X));
                case 'stdbin'
                    Res(Ibin,Icol) = std(Data(Flag,Col.X));
                case 'wmean'
                    Res(Ibin,Icol) = Util.stat.wmean(Data(Flag,[Col.X,Col.E]));
                otherwise
                    error('Unknown Column option');
            end
        else
            % Col is a function handle
            switch lower(func2str(OutCol{Icol}))
                case {'util.stat.wmean','util.stat.wmedian'}
                    Res(Ibin,Icol) = OutCol{Icol}(Data(Flag,Col.Y), Data(Flag,Col.E));
                case {'werr', 'wmeanerr'} 
                    [~,Res(Ibin,Icol)] = Util.stat.wmean(Data(Flag,Col.Y), Data(Flag,Col.E));
                otherwise
                    Res(Ibin,Icol) = OutCol{Icol}(Data(Flag,Col.Y));
            end
        end
    end
end

% Convert function names in columns to strings
for Icol=1:1:Ncol
    if (isa(OutCol{Icol},'function_handle'))
        OutCol{Icol} = func2str(OutCol{Icol});
    end
end
            
switch lower(OutType)
    case 'astcat'
        % Convert matrix to AstCat object
        Res = AstCat.array2astcat(Res,OutCol);
    otherwise
        % do nothing - output is a matrix
end

