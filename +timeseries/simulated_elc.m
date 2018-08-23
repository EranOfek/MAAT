function [Events,Nbe]=simulated_elc(LC,Window,Back)
% Simulated photons light curve
% Package: timeseries
% Description: Given a model light curve, generate a list of time-tagged
%              events that their rate follows the model light curve.
%              The events are generated between the first and last time
%              in the model light curve.
%              Note that the first event is always at the first time
%              in the model light curve.
%              The rate at each time is calculated using a second
%              order approximation: Tau0*(1 + dTau/dt + [dTau/dt]^2).
% Input  : - Light curve model [time, Lambda], where Lambda is the mean
%            rate (>0) at each time (e.g., events per unit time).
%          - Matrix, containing in each line the start and end times of
%            active windows. Events which thier time tags is outside
%            these windows are deleted.
%            If empty matrix (e.g., []) then, all times are used.
%            Default is [].
%          - Rate (per unit time) of optional background.
%            If empty matrix (e.g., []), then no background. Default is [].
% Output : - Sorted time tagged events.
%          - Number of backround events.
% Tested : Matlab 7.3
%     By : Eran O. Ofek                    Jan 2008
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Events,Nbe]=timeseries.simulated_elc([1 10; 2 10; 3 500; 4 10; 5 10],[1 5],5);
% Reliable: 2
%------------------------------------------------------------------------------

InterpMethod = 'linear';

Col.Time   = 1;
Col.Lambda = 2;

Def.Back   = [];
Def.Window = [];
if (nargin==1),
   Window    = Def.Window;
   Back      = Def.Back;
elseif (nargin==2),
   Back      = Def.Back;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

Events = LC(1,Col.Time);
Tau    = 1./LC(:,2);
while (Events(end)<LC(end,Col.Time)),
   J =  (LC(:,Col.Time) - Events(end))>0;
   [NearestDelta] = min(LC(J,Col.Time) - Events(end));
   InterpTau = interp1(LC(:,Col.Time),Tau,Events(end)+[0; NearestDelta],InterpMethod);

   DTauDt = (InterpTau(2)-InterpTau(1))./NearestDelta;
   NextTau = InterpTau(1).*(1 + DTauDt + DTauDt.^2);

   Events = [Events; Events(end)+exprnd(NextTau,1,1)];
end

Events = Events(1:end-1);


if (isempty(Back)==1),
   % do nothing
   Nbe = 0;
else
   TimeRange  = range(LC(:,Col.Time));
   Nbe        = ceil(TimeRange.*Back);  % number of background events
   BackEvents = rand(Nbe,1).*TimeRange + min(LC(:,Col.Time));

   Events  = sort([Events; BackEvents]);
end

if (isempty(Window)==1),
   % do nothing
else
   ListGood = [];
   for I=1:1:size(Window,1),
      ListGood = [ListGood; find(Events>=Window(I,1) & Events<=Window(I,2))];
   end
   Events = sort(Events(ListGood));
end
