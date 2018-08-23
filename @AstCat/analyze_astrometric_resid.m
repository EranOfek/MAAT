function [ResAnalysis]=analyze_astrometric_resid(Cat,Ref,ImSize,BlockSize)
%--------------------------------------------------------------------------
% analyze_astrometric_resid function                         class/@AstCat
% Description: 
% Input  : - AstCat object in which the first three columns are [X,Y,Mag].
%          - AstCat object for the reference catalog in which the first
%            fourcolumns are [X,Y,Mag,Color].

%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ColXc  = 1;
ColYc  = 2;
ColMc  = 3;
ColXr  = 1;
ColYr  = 2;
ColMr  = 3;
ColCr  = 4;

% general correlations and trends
DeltaX = Cat.Cat(:,ColXc) - Ref.Cat(:,ColXr);
DeltaY = Cat.Cat(:,ColYc) - Ref.Cat(:,ColYr);
N      = numel(DeltaX);

% correlations of the residuals
ResAnalysis.Corr.X_dX = corr(Ref.Cat(:,ColXr),DeltaX);
ResAnalysis.Corr.Y_dX = corr(Ref.Cat(:,ColYr),DeltaX);
ResAnalysis.Corr.X_dY = corr(Ref.Cat(:,ColXr),DeltaY);
ResAnalysis.Corr.Y_dY = corr(Ref.Cat(:,ColYr),DeltaY);
ResAnalysis.Corr.M_M  = corr(Ref.Cat(:,ColMr),Cat.Cat(:,ColMc));
ResAnalysis.Corr.M_dX = corr(Ref.Cat(:,ColMr),DeltaX);
ResAnalysis.Corr.M_dY = corr(Ref.Cat(:,ColMr),DeltaY);
ResAnalysis.Corr.C_dX = corr(Ref.Cat(:,ColCr),DeltaX);
ResAnalysis.Corr.C_dY = corr(Ref.Cat(:,ColCr),DeltaY);


% rms in blocks

Merged = col_concat(Cat,Ref,[1 2 3],[1 2 3 4]);

[SubCat,ListEdge,ListCenter]=subcat_regional(Merged,ImSize,BlockSize,0,[1 2]);
Nsub = numel(SubCat);
for Isub=1:1:Nsub,
    DeltaX = SubCat(Isub).Cat(:,1) - SubCat(Isub).Cat(:,4);
    DeltaY = SubCat(Isub).Cat(:,2) - SubCat(Isub).Cat(:,5);
    
    ResAnalysis.Block(Isub).X    = ListCenter(Isub,1);
    ResAnalysis.Block(Isub).Y    = ListCenter(Isub,2);
    ResAnalysis.Block(Isub).Sec  = ListEdge(Isub,:);
    ResAnalysis.Block(Isub).N    = numel(DeltaX);
    ResAnalysis.Block(Isub).rmsX = std(DeltaX);
    ResAnalysis.Block(Isub).rmsY = std(DeltaY);
    ResAnalysis.Block(Isub).rms  = std(sqrt(DeltaX.^2 + DeltaY.^2));
    ResAnalysis.Block(Isub).meanX= mean(DeltaX);
    ResAnalysis.Block(Isub).meanY=mean(DeltaY);
    ResAnalysis.Block(Isub).maxX = max(abs(DeltaX));
    ResAnalysis.Block(Isub).maxY = max(abs(DeltaY));
end

    
ResAnalysis.WorseBlock.rmsX   = max([ResAnalysis.Block.rmsX]);
ResAnalysis.WorseBlock.rmsY   = max([ResAnalysis.Block.rmsY]);
ResAnalysis.WorseBlock.rms    = max([ResAnalysis.Block.rms]);
ResAnalysis.WorseBlock.meanX  = max(abs([ResAnalysis.Block.meanX]));
ResAnalysis.WorseBlock.meanY  = max(abs([ResAnalysis.Block.meanY]));


    
    


