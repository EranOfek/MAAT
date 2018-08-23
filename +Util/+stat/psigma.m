function P=psigma(Sigma,Sides)
%--------------------------------------------------------------------------
% psigma function                                                AstroStat
% Description: Return the two sided or one sided probability for a given
%              sigma level.
% Input  : - Sigma level.
%          - Number of sides {1|2}, default is 2.
%            For example, 1 sigma for 2 sides is 0.6827
%            and for 1 side is 0.8413.
% Output : - Probability.
% Tested : Matlab 7.6
%     By : Eran O. Ofek                    Mar 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: P=psigma(1,2);
% Reliable: 1
%--------------------------------------------------------------------------
DefSides = 2;
if (nargin==1),
   Sides = DefSides;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

switch Sides
 case 2
    P = 1-(1-normcdf(Sigma,0,1)).*2;
 case 1
    P = normcdf(Sigma,0,1);
 otherwise
    error('Unknown Sides option');
end
