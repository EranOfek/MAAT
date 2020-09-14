function B=braek_astcat2zones(AC,varargin)
% Break an AstCat object to multiple files, sorted by one of the columns.
% Package: @AstCat
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;
CatField = AstCat.CatField;

DefV.BreakCol             = 'Dec';
DefV.BreakStep            = 1./RAD;   % deg
DefV.StepBuffer           = 0./RAD;
DefV.CatUnits             = 'rad';
DefV.SaveName             = 'GALEX2mEpochs_%s%02d_%s%02d.mat';
DefV.SignM                = 'm';
DefV.SignP                = '';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

AC = sortrows(AC,InPar.BreakCol);
Col = AC.SortedByCol;

Start = min(AC.(CatField)(:,Col))-10.*eps;
End   = max(AC.(CatField)(:,Col))+10.*eps;

ConvF = convert.angular(InPar.CatUnits,'deg');

VecZone = (Start:InPar.BreakStep:End).';  % rad
VecZone = floor(VecZone.*ConvF);
VecZone = [VecZone; VecZone(end)+InPar.BreakStep.*ConvF];
VecZone = VecZone./ConvF;

Nzone   = numel(VecZone);

if (nargout>0)
    B = AstCat(Nzone-1,1);
end
for Izone=1:1:Nzone-1
    Izone
    V1 = VecZone(Izone);
    V2 = VecZone(Izone+1);
    
    Flag = AC.(CatField)(:,Col)>(V1-InPar.StepBuffer) & ...
           AC.(CatField)(:,Col)<(V2+InPar.StepBuffer);
       
    Tmp = row_select(AC,Flag);
%     if (nargout>0)
%         B(Izone) = Tmp;
%     end
    % save the Tmp
    if ~isempty(InPar.SaveName)
        V1N = V1.*ConvF;
        V2N = V2.*ConvF;
        if (V1N)>0
            Sign1 = InPar.SignP;
        else
            Sign1 = InPar.SignM;
        end
        if (V2N)>0
            Sign2 = InPar.SignP;
        else
            Sign2 = InPar.SignM;
        end
        FileName = sprintf(InPar.SaveName,Sign1,round(abs(V1N)),Sign2,round(abs(V2N)));
        save(FileName,'Tmp','-v7.3');
        
    end
end
    
    