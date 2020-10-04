function [W,ResW]=weight_cadence(JD,LastObs,NightCounter,varargin)
% Calculate the cadence weight for a list of targets.
% Package: +celestial.scheduling
% Description: The cadence weight is a combination of two weights, based on
%              the last time the target was observed during the night, and
%              the last time it was observed before the current night.
%              This is calculated using the current JD, the JD of last
%              observation, and the nightly counter.
% Input  : - The current JD.
%          - A vector of JD on which each target was last observed.
%          - A vector of night counter. The number of times each target was
%            observed on the current night.
%          * Arbitrary number of pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'MainCadence' - The main cadence of the survey [day]
%                       Default is 2.4.
%            'NightCadence' - The nightly cadence [day].
%                       Default is 40./1440.
%            'Nfast' - Number of epochs per night in the nightly cadence.
%                       Default is 2.
%            'MainWFun' - Weight function for main cadence.
%                       Default is @celestial.scheduling.fermiexp
%            'MainWFunPar' - Parameters for weight function for main
%                       cadence.
%                       Default is {2.4, 1, 0.03, 1, 0.5}.
%            'NightWFun' - Weight function for nightly cadence.
%                       Default is @celestial.scheduling.fermiexp
%            'NightWFunPar' - Parameters for weight function for nightly
%                       cadence.
%                       Default is {40./1440, 1, 0.003, 1.5, 0.5}.
% Output : - Vector of weights for each target.
%          - A structure with additional information.
%            The following fields are available:
%            .TimeSinceLastObs
%      By: Eran Ofek                      Oct 2020
% Example : 




RAD = 180./pi;

InPar = inputParser;
addOptional(InPar,'MainCadence',2.4);  % [day]
addOptional(InPar,'NightCadence',40./1440); % [day]
addOptional(InPar,'Nfast',2); % [day]
addOptional(InPar,'MainWFun',@celestial.scheduling.fermiexp); %@(t) 1.0+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
addOptional(InPar,'MainWFunPar', {2.4, 1, 0.03, 1, 0.5} );  %t0, DecayExp, SoftFermi, BaseW, ExtraW
addOptional(InPar,'NightWFun',@celestial.scheduling.fermiexp); %@(t) 1.5+0.5.*exp(-t./1) ); % weight as a function of time since it is allowed to observe the target
addOptional(InPar,'NightWFunPar',{40./1440, 1, 0.003, 1.5, 0.5});  %t0, DecayExp, SoftFermi, BaseW, ExtraW

parse(InPar,varargin{:});
InPar = InPar.Results;

% time since last observation
TimeSinceLastObs = JD - LastObs;
% flag indicating objects that requires the next observation on the same
% night
FlagNight = NightCounter>0 & NightCounter<InPar.Nfast & TimeSinceLastObs>InPar.NightCadence;
W_Night   = InPar.NightWFun(TimeSinceLastObs, InPar.NightWFunPar{:}) .* FlagNight;


% flag indicating objects that requires the first observation on the night
FlagMain    = NightCounter==0 & TimeSinceLastObs>InPar.MainCadence;
W_Main      = InPar.MainWFun(TimeSinceLastObs, InPar.MainWFunPar{:}) .* FlagMain;

W = max(W_Main, W_Night);
ResW.TimeSinceLastObs = TimeSinceLastObs;
