function [NewRA,NewDec]=interp_coo(Time,RA,Dec,NewTime,InterpMethod)
% Interpolate celestial coordinates as a function of time
% Package: celestial.coo
% Description: Interpolate on celestial ccordinates as a function of time.
%              Use the built in matlab interpolation functions.
% Input  : - Vector of Times (e.g., JD).
%          - Vector of longitudes [radians].
%          - Vector of latitudes [radians].
%          - Vector of new times in which to interpolate position.
%          - Algorithm (see interp1.m for details), default is 'pchip'.
% Output : - Vector of longitudes interpolated to the to the vector
%            of new times.
%          - Vector of latitudes interpolated to the to the vector
%            of new times.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Mar 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Time=(1:1:11)';
%          RA = (0.1:0.01:0.2)'; Dec = (0.1:0.01:0.2).';
%          NewTime = [2.3;3.4];
%          [NewRA,NewDec]=celestial.coo.interp_coo(Time,RA,Dec,NewTime);
% Reliable: 2
%--------------------------------------------------------------------------

DefInterpMethod = 'pchip';
if (nargin==4)
   InterpMethod = DefInterpMethod;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end

CD = celestial.coo.cosined([RA, Dec]);
NewCD = zeros(length(NewTime),3);

NewCD(:,1) = interp1(Time,CD(:,1),NewTime,InterpMethod);
NewCD(:,2) = interp1(Time,CD(:,2),NewTime,InterpMethod);
NewCD(:,3) = interp1(Time,CD(:,3),NewTime,InterpMethod);

NewCoo = celestial.coo.cosined(NewCD);
NewRA  = NewCoo(:,1);
NewDec = NewCoo(:,2);
