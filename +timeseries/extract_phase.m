function OutMat=extract_phase(InMat,T0,Period,PhaseRange)
% Extract observations take at some phase range
% Package: timeseries
% Description: Given a time series, extract observation made in a given
%              phases range. The phase and epoch is defined by the user.
% Input  : - observations matrix in which first column is time and
%            second column is observed value.
%          - ephemeris start time, T0.
%          - ephemeris period.
%          - Vector of phase range: [From To].
% Output : - Observations matrix, containing only the observations
%            for which the phase restriction are fulfilled.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OutMat=timeseries.extract_phase(rand(100,2),0,0.1,[0.1 0.2]);
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin~=4),
   error('Illegal number of input arguments');
end


PhTemp = (InMat(:,1) - T0)./Period;
Phase  = PhTemp - floor(PhTemp);

if (PhaseRange(2)>PhaseRange(1)),
   I = (Phase>=PhaseRange(1) & Phase<=PhaseRange(2));
else
   I = (Phase>=PhaseRange(1) | Phase<=PhaseRange(2));
end

OutMat = InMat(I,1:end);

