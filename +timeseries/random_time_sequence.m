function T=random_time_sequence(TotalT,Cadence,LenAnnual,RandStd,MonthlyFrac)
% Generate random times for an astronomical time series.
% Package: timeseries
% Description: Generate random times for an astronomical time series,
%              including daily and annual gaps.
% Input  : - Total time span [days]. Default is 365.*3.
%          - Mean cadence [days]. Default is 2.
%          - Length of annual obsrving period [days]. Default is 240.
%          - Random std for observing times [days]. Default is 0.05.
%          - Montly observing fraction. Default is 0.8.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=timeseries.random_time_sequence
% Reliable: 
%--------------------------------------------------------------------------
Year = 365.25;
Month = 29.53;

if (nargin==0)
    TotalT      = 365.*3;
    Cadence     = 2;
    LenAnnual   = 240;
    RandStd     = 0.05;
    MonthlyFrac = 0.8;
end

T = (1:Cadence:TotalT).';
ModYear = mod(T,365.25);
FlagY = ModYear<LenAnnual;

ModMonth = mod(T,Month)./Month;
FlagM = ModMonth<MonthlyFrac;

T = T(FlagY & FlagM);

T = T + randn(size(T)).*RandStd;




