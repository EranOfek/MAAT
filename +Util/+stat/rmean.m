function Mean=rmean(Mat,Dim,Range)
% Calculate the rubust mean over one of the dimensions.
% Package: Util.stat
% Description: Calculate the rubust mean over one of the dimensions.
% Input  : - A matrix.
%          - Dimension along to calculate the ribust mean.
%            If empty calc over all dimensions.
%            Deafult is 'all'.
%          - [Low, High] fraction of the values to remove prior to the
%            mean calculation. Default is [0.05 0.05].
%            I.e., remove the loewer and upper 5 percentile and
%            calculate the mean.
% Output : - The mean.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: A=rand(10000,4); A(1,:)=1e7; rmaen(A)
% Reliable: 2
%--------------------------------------------------------------------------

Def.Dim   = [];
Def.Range = [0.05 0.05];
if (nargin==1)
    Dim   = Def.Dim;
    Range = Def.Range;
elseif (nargin==2)
    Range = Def.Range;
elseif (nargin==3)
    % do nothing
else
    error('Illegal number of input arguments');
end

if ischar(Dim)
    % assume all
    Mat = Mat(:);
    SortedMat = sort(Mat);
    N         = numel(Mat);

else
    SortedMat = sort(Mat,Dim);
    N         = size(Mat,Dim);
end

Low  = ceil(N.*Range(1));
High = floor(N.*(1-Range(2)));

Low  = min(N,Low+1);
if (Low>=High)
    warning('Robust mean return empty as the sample is too small');
    Mean = [];
else
    if ischar(Dim)
        Mean = nanmean(SortedMat(Low:High,:),'all');
    else
        
        if (Dim==1)
            Mean = nanmean(SortedMat(Low:High,:),1);
        elseif (Dim==2)
            Mean = nanmean(SortedMat(:,Low:High),2);
        elseif (Dim==3)
            Mean = nanmean(SortedMat(:,:,Low:High),3);
        else
            error('Dim must be 1,2 or 3');
        end
    end
end


