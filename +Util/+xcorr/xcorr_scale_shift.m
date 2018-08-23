function Info=xcorr_scale_shift(Data,Template,varargin)
%--------------------------------------------------------------------------
% xcorr_scale_shift function                                       General
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Outout : - A structure containing the best shift and scale between
%            the two serieses.
%            .DataTemp - A vector with the length of the first input
%                        argument (Data), with coordinates corresponding
%                        to the template coordinates.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Template = SpecArcs(9).Spec;
%          Data = interp1(Template(:,1),Template(:,2),(3300:0.7:5400).');
%          Data = [(1:1:length(Data)).'.*0.75+10,Data];
%          Info=xcorr_scale_shift(Data,Template)
%          X1=(1:1.3:100).'; Y1=zeros(size(X1)); Y1(10)=1; Y1(20)=1;
%          X2=(1:0.9:100).'; Y2=zeros(size(X2)); Y2(12)=1; Y2(30)=1;
%          Data = [X1,Y1]; Template = [X2, Y2];
%          Info=xcorr_scale_shift(Data,Template)
% Reliable: 2
%--------------------------------------------------------------------------
Col.X = 1;
Col.Y = 2;


DefV.MaxRange      = 1e-8;
DefV.InterpMethod  = 'linear';
DefV.ScaleVec      = logspace(-1,1,1000);
DefV.Back          = 'median';
DefV.BackPar       = {};
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% make ScaleVec a row vector
if (size(InPar.ScaleVec,1)>1),
    InPar.ScaleVec = InPar.ScaleVec.';
end

% check if Data series is equally space
if (range(unique(diff(Data(:,Col.X))))<InPar.MaxRange)
    % equally spaced
else
    % make the series equally spaced
    X = (min(Data(:,Col.X)):min(diff(Data(:,Col.X))):max(Data(:,Col.X))).';
    Y = interp1(Data(:,Col.X),Data(:,Col.Y),X,InPar.InterpMethod);
    Data = [X,Y];
end
% check if Template series is equally space
if (range(unique(diff(Template(:,Col.X))))<InPar.MaxRange)
    % equally spaced
else
    % make the series equally spaced
    X = (min(Template(:,Col.X)):min(diff(Template(:,Col.X))):max(Template(:,Col.X))).';
    Y = interp1(Template(:,Col.X),Template(:,Col.Y),X,InPar.InterpMethod);
    Template = [X,Y];
end


% build a version of the template with all possible stretches
Ndata = size(Data,1);
Ntemp = size(Template,1);
DataI = (1:1:Ndata).';
TempI = (1:1:Ntemp).';

StTemplateX = bsxfun(@times,TempI,InPar.ScaleVec);
StTemplateY = interp1(TempI,Template(:,Col.Y),StTemplateX,InPar.InterpMethod,'extrap');


%Data(:,2)     = Data(:,2).*DataI;
%Template(:,2) = Template(:,2).*TempI; 

[~,~,Info]         = xcorr_fft_multi(StTemplateY,Data(:,Col.Y),InPar.Back,InPar.BackPar);

Info.BestScale     = InPar.ScaleVec(Info.BestCol);
Info.ScaleData     = mean(diff(Data(:,Col.X)));
Info.ScaleTemplate = mean(diff(Template(:,Col.X)));


%DataItemp = (Info.BestScale*DataI+Info.BestShift.*Info.BestScale);

Info.DataTemp = (DataI + Info.BestShift).*Info.BestScale.*Info.ScaleTemplate+min(Template(:,Col.X));


% plot(Info.DataTemp,Data(:,2));
% hold on;
% plot(Template(:,1),Template(:,2),'r-');


