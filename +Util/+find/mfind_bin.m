function Im=mfind_bin(X,Vals)
% Binary search on a vector running simolutnously on multiple values.
% Package: Util.find
% Description: Binary search on a vector running simolutnously on
%              multiple values. A feature of this program is that it
%              you need to add 1 to the index in order to make sure 
%              the found value is larger than the searched value.
% Input  : - Sorted column vector.
%          - Row vector of values to search.
% Output : - Indices of nearest values.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=sort(rand(1e6,1)); Vals=rand(1,1e5);
%          Im=Util.find.mfind_bin(X,Vals)
% Reliable: 2
%--------------------------------------------------------------------------

Nvals = length(Vals);
N     = length(X);
I1    = ones(1,Nvals);
I2    = N.*ones(1,Nvals);
Im    = uint32(floor(0.5.*(I1+I2)));
PrevIm= zeros(size(Im),'uint32');

if (numel(X)<2)
    if (isempty(X))
        Im = uint32([]);
    else
        Im = ones(1,Nvals);
    end
else
    %Niter = 0;
    while (~all(Im==PrevIm))
        %Niter = Niter+ 1;
        FlagU = Vals>X(Im).';
        FlagD = ~FlagU; %Vals<X(Im).';
        I1(FlagU) = Im(FlagU);
        I1(FlagD) = I1(FlagD);
        I2(FlagU) = I2(FlagU);
        I2(FlagD) = Im(FlagD);
        PrevIm    = Im;
        Im        = uint32(floor(0.5.*(I1+I2)));
    end
end