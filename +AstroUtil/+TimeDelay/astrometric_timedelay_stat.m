function []=astrometric_timedelay_stat(Time,TotF,ErrF,X,Y,ErrX,ErrY,TimeDelay,FluxRtaio,PosX,PosY,varargin)
% Calculate the astrometric time delay statistics for a give time delay
% Package: AstroUtil.TimeDelay
% Description: Given total flux and center of mass (astrometric)
%              observations, calculate the astrometric time delay
%              statistics for a give time delay and flux ratio.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



DefV.Resample             = false;
DefV.ResamplePar          = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% sort by time
[Time,SI] = sort(Time);
TotF      = TotF(SI);
ErrF      = ErrF(SI);
X         = X(SI);
Y         = Y(SI);
ErrX      = ErrX(SI);
ErrY      = ErrY(SI);

% resample LC to an equally spaced LC
if InPar,Resample
    % time series is not equally spaced
    % resample
    [Prop,Time] = timeseries.resample_uniform(Time,[TotF,ErrF,X,Y,ErrX,ErrY],InPar.ResamplePar{:});
    TotF = Prop(:,1);
    ErrF = Prop(:,2);
    X    = X(:,2);
    Y    = Y(:,2);
    ErrX = ErrX(:,1);
    ErrY = ErrY(:,2);
    
end

% recover source light curve
[~,fftSourceF,SourceF,ShiftFactor]=AStroUtil.TimeDelay.recover_source_flux([Time,TotF,ErrF],TimeDelay,FluxRatio,'Resample',[])

PosX(1)./ShiftFactor 
fftA = 



