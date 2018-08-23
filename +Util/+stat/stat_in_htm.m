function Res=stat_in_htm(Data,varargin)
% SHORT DESCRIPTION HERE
% Package: Util
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Util.stat.stat_in_htm(Table(:,[1 2 3]));
% Reliable: 
%--------------------------------------------------------------------------

DefV.Fun                  = {@nanmean,@nanmedian,@nanstd,@numel};
DefV.ColRA                = 1;
DefV.ColDec               = 2;
DefV.ColProp              = 3;
DefV.HTM_Level            = 4;
DefV.CooUnits             = 'deg';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Data(:,InPar.ColRA)  = convert.angular(InPar.CooUnits,'rad',Data(:,InPar.ColRA));
Data(:,InPar.ColDec) = convert.angular(InPar.CooUnits,'rad',Data(:,InPar.ColDec));


[HTM,LevList]=celestial.htm.htm_build(InPar.HTM_Level);

Ndata = size(Data,1);

DataInd = zeros(Ndata,1);
for Idata=1:1:Ndata
    DataInd(Idata) = celestial.htm.htm_search_point(HTM,Data(Idata,[InPar.ColRA, InPar.ColDec]));
    
end


Nfun    = numel(InPar.Fun);
FunName = cell(1,Nfun);
for Ifun=1:1:Nfun
    FunName{Ifun} = func2str(InPar.Fun{Ifun});
end
    

Nptr = numel(LevList(end).ptr);
Res.MeanRA  = zeros(Nptr,1);
Res.MeanDec = zeros(Nptr,1);
Res.FunName = FunName;
Res.Prop    = nan(Nptr,Nfun);

for Ihtm=1:1:Nptr
    Ptr = LevList(end).ptr(Ihtm);
    Res.MeanRA(Ihtm)  = mean(HTM(Ptr).coo(:,1));
    Res.MeanDec(Ihtm) = mean(HTM(Ptr).coo(:,2));
    
    If = find(Ptr == DataInd);
    
    for Ifun=1:1:Nfun
        Res.Prop(Ihtm,Ifun) = InPar.Fun{Ifun}(Data(If,InPar.ColProp));
    end
end
        
    
    
    


