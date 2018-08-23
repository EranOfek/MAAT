function TaperVal=taper(X,varargin)
% Generate a taper function
% Package: timeseries
% Description: Generate a taper function for a timeseries. Taper function is
%              a weight/window function in the time domain.
%              This can be used in order to give reduce weight to data
%              found near the edges of a timeseries in order to reduce
%              aliasing with the series length.
% Input  : - A column vector of X (e.g., time). If multiple columns are
%            given, then the program will use the first column only.
%          * Arbitrary number of pairs of ...,key,val,...
%            The following keywords are available:
%            'TaperFun' - Name of taper function
%                         'cosbell' - cosine bell (default).
%                         'trapz'   - a trapzoid.
%                         Alternatively, this can be a nuerical function
%                         [X,Y] norzmlized that such X is in the range 0 to 1.
%                         Alternatively, this can be a function handle.
%                         Y=@fun(X), where X is in the rane 0 to 1.
%            'TaperPar' - Parameters of the taper function.
%                         For 'cosbell' and 'trapz' this is the precentage
%                         of flat part of the cosine bell (in the range 0..1).
%                         Where the total range of the taper function is
%                         the total range of X.
%                         Default is 0.9.
%                         If TaperFun is a function handle, this can be
%                         a cell array of additional parameters that will be
%                         passed to the function.
%            'Norm'     - {'y'|'n'}. Normalize the taper function such that
%                         the maximum height of the taper is 1.
%                         This is useful if TaperFun is an unnormalized
%                         function. Default is 'n'.
%            'Y'        - Y value of the time series. If this parameter is
%                         provided, then the output argument is the taper
%                         function multiplied by the Y value.
%                         Default is Y=ones(size(X));
%            'InterpMethod' - Interpolation method to be used if the taper
%                         function is numeric. See interp1.m for options.
%                         Default is 'linear'.
% Output - The value of the taper function for each X value in the input
%          vector. If The keyword 'Y' is provided then this will be the
%          taper function multiply by the Y value.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X = (100:1:200)';
%          TaperVal=timeseries.taper(X);
%          TaperVal=timeseries.taper(X,'TaperPar',0.0);
%          TaperVal=timeseries.taper(X,'TaperPar',0.5,'TaperFun','trapz','Y',ones(size(X)).*2);
%          TaperVal=timeseries.taper(X,'TaperFun',@sin);
% Reliable: 2
%------------------------------------------------------------------------------

Xdata = X(:,1);

DefV.TaperFun     = 'cosbell';
DefV.TaperPar     = 0.9;
DefV.Norm         = 'n';
DefV.Y            = ones(size(Xdata));
DefV.InterpMethod = 'linear';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);                                                  

RangeX = range(Xdata);
Xs     = min(Xdata);
Xe     = max(Xdata);
Xfrac  = (Xdata-Xs)./RangeX;

TaperVal = ones(size(Xdata));
if (ischar(InPar.TaperFun)),
   switch lower(InPar.TaperFun)
    case 'cosbell'
       Xstart = 0.5.*(1-InPar.TaperPar);      % start of flat part
       Xend   = 1 - 0.5.*(1-InPar.TaperPar);  % end of flat part

       Irise = find(Xfrac<Xstart);
       Iset  = find(Xfrac>Xend);
       Iflat = find(Xfrac>=Xstart & Xfrac<=Xend);

       TaperVal(Irise) = sin(Xfrac(Irise)./Xstart .*pi./2);
       TaperVal(Iset)  = cos((Xfrac(Iset)-Xend)./Xstart .*pi./2);
    case 'trapz'
       Xstart = 0.5.*(1-InPar.TaperPar);      % start of flat part
       Xend   = 1 - 0.5.*(1-InPar.TaperPar);  % end of flat part

       Irise = find(Xfrac<Xstart);
       Iset  = find(Xfrac>Xend);
       Iflat = find(Xfrac>=Xstart & Xfrac<=Xend);

       TaperVal(Irise) = Xfrac(Irise)./Xstart;
       TaperVal(Iset)  = 1-(Xfrac(Iset)-Xend)./Xstart;
    otherwise
       error('Unknown TaperFun option');
   end
elseif (isnumeric(InPar.TaperFun)),
   TaperVal = interp1(InPar.TaperFun(:,1),InPar.TaperFun(:,2),Xfrac,InPar.InterpMethod);
elseif (isa(InPar.TaperFun,'function_handle')),
   if (iscell(InPar.TaperPar)),
      TaperVal = feval(InPar.TaperFun,Xfrac,InPar.TaperPar{:});
   else
      TaperVal = feval(InPar.TaperFun,Xfrac);
   end
else
   error('Illegal TaperFun option');
end

switch lower(InPar.Norm)
 case 'y'
    TaperVal = TaperVal./max(TaperVal);
 otherwise
    % do nothing
end   

	   

% applay taper
TaperVal = TaperVal.*InPar.Y;

