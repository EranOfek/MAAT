function [Res,Summary,N_Ep]=astcat2matched_array(AstC,Cols,OutType)
% AstCat object to a matrix of matched sources.
% Package: @AstCat
% Description: Given an AstCat object containing multiple element, in
%              which each element containing the same number of rows
%              (e.g., the output of AstCat/match.m), return a matrix
%              that contains, for selected columns, the requested column
%              in all the AstCat elements. E.g., the number of columns
%              in the matrix is equal to the number of AstCat elements
%              (column per element) and the number of rows is equal to
%              the number of rows in each AstCat element.
% Input  : - An AstCat object containing multiple element, in
%            which each element containing the same number of rows
%            (e.g., the output of AstCat/match.m)
%          - Column indices or column names (string or a cell array of
%            strings) for which to prepare the array.
%          - A flag indicating if the ouput matrix is
%            [Nepoch X Nsrc] (true), or
%            [Nsrc X Nepoch] (false).
%            Default is false.
% Output : - A structure in which each field name corresponds to a
%            requested column name and contains the matrix of all the
%            column entries in all the AstCat object elements.
%          - A structure array containing a summary.
%            The following fields are available:
%            .Nsrc - Number of sources (size of the 1st dimension of the
%                    output matrix).
%            .Nepoch - Number of epochs (size of the 2nd dimension of the
%                    output matrix).
%            .Nnn - number of not NaNs for each source in the first Column.
%          - Vector of length equal to the number of epochs. The value in
%            each element is the number of stars that appears in exactly
%            I epochs.
% See also: match.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: M=match(Ref,AstC);
%          Res=astcat2matched_array(M,{'XWIN_IMAGE','YWIN_IMAGE','MAG_APER'});
%          % step by step application:

%          S=images2sim('PTF_2015*.fits');
%          S=gain_correct(S);
%          S=background(S,'Block',[64 64]);
%          S=mextractor(S);
%          [~,I]=max(sizecat(S)); 
%          Sref = S(I);
%          [M,UM]=match(S,Sref);
%          [Res,Summary]=astcat2matched_array(M,{'MAG_PSF'});
%          II=find(Summary.Nnn==Summary.Nepoch); % find sources that appears in all epochs
%          Res.MAG_PSF(:,II)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2)
    OutType = false;
end

CatField   = 'Cat';

Ncat        = numel(AstC);
[CatRow,~]  = sizecat(AstC);
if all(CatRow(1)==CatRow)
    Nrow = CatRow(1);
else
    error('Number of rows in all AstCat elemnts must be equal');
end
    

if (~iscell(Cols))
    Cols = {Cols};
end

% get column names and indices
ColInd  = colname2ind(AstC(1),Cols);
ColName = ind2colname(AstC(1),Cols);
Ncol    = numel(ColInd);

% Initiate Res
for Icol=1:1:Ncol
    Res.(ColName{Icol}) = zeros(Nrow,Ncat);
end

if (OutType)
    DimSrc   = 2;
    DimEpoch = 1;
else
    DimSrc   = 1;
    DimEpoch = 2;
end

% For each selected column
for Icol=1:1:Ncol
    % for each catalog
    for Icat=1:1:Ncat
        % E.g., Res.XWIN_IMAGE(:,10) = AstC(10).Cat(:,Column)
        Res.(ColName{Icol})(:,Icat) = AstC(Icat).(CatField)(:,ColInd(Icol));
    end
    if (OutType)
        % transpose the matrix
        Res.(ColName{Icol}) = Res.(ColName{Icol}).';
    end
end


if (nargout>1)
    % calculate summary
    Icol = 1;
    Summary.Nsrc   = size(Res.(ColName{Icol}),DimSrc);
    Summary.Nepoch = size(Res.(ColName{Icol}),DimEpoch);
    Summary.Nnn    = sum(~isnan(Res.(ColName{Icol})),DimEpoch);
    
    if (nargout>2)
        % calculate the number of stars that appears in N images
        N_Ep = zeros(Summary.Nepoch,1);
        for Iep=1:1:Summary.Nepoch
            % Summary.Nnn is the number of epochs in each stars appears
            % N_Ep is the number of stars that appears in exactly Iep
            % epochs.
            N_Ep(Iep) = sum(Summary.Nnn==Iep);
        end
    end
end

    
    
