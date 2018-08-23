function Matched=mat2matched_array(ColNames,varargin)
% Store the columns of matched arrays in specific matrices.
% Package: ImUtil
% Description: Given a list of matrices with the same number of columns and
%              rows, and in which the rows are matched, store the i'th
%              column of each matrix in a single matrix. If the matrices
%              represent different epochs, then there is an output matrix
%              (stored in a structure) per each column and the columns in
%              the output matrix represent different epochs.
% Input  : - Cell array of column names.
%          * Arbitrary number of matrices.
% Output : - A structure, with field names equal to the column names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Matched=ImUtil.pattern.mat2matched_array(ColNames,M1,M2)
% Reliable: 2
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Nep  = numel(varargin);
if (Nep<1)
    error('Must supply at least one epoch');
end
Ncol = numel(ColNames);
for Icol=1:1:Ncol
    
    Nsrc = size(varargin{1},1);
    Matched.(ColNames{Icol}) = nan(Nsrc,Nep);
    for Iep=1:1:Nep
        Matched.(ColNames{Icol})(:,Iep) = varargin{Iep}(:,Icol);
    end
end
