function [Interp]=interpolant_mag(Data,varargin)
% Interpolant from a time series of photometric observations in one band 
% Package: AstroUtil.spec
% Description: Given a time series of observations take at a single band
%              return an interpolant that allows to calculate the magnitude
%              in each time within observations range.
% Input  : - [Time, Mag, Err]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MaxTimeNoObs' - Maximum time without observations.
%                             If gap larger than this exist in the data the
%                             interpolant will not be valid in the gap.
%                             Default is 5.
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Interp=AstroUtil.spec.interpolant_mag(Data)
% Reliable: 
%--------------------------------------------------------------------------

Col.T = 1;
Col.M = 2;
Col.E = 3;


DefV.MaxTimeNoObs         = 5;   % units of time
DefV.ExtrapOutOfRange     = 0.2; % units of time
DefV.MaxDeg               = 4;   % max polynomial degree
DefV.PolyDegTable         = [0 0.2 0;0.2 1 1;1 2 2;2 5 3;5 10 4; 10 Inf 5];        % (timeSpan), (timeSpan), Deg
DefV.LogTime              = false;
DefV.TableStep            = 0.01;
DefV.Plot                 = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (size(Data,2)<3)
    Data = [Data, ones(size(Data,1),1)];
end

Data = sortrows(Data);
Nobs = size(Data,1);
DT   = diff(Data(:,Col.T));  % delta T

Igap = find(DT>InPar.MaxTimeNoObs);
Igap = [1;Igap+1];
Nsec = numel(Igap); % number of continous sections (with no gap)
for Isec=1:1:Nsec
   
    % Is is the indices of data in continous section
    if (Isec<Nsec)
        Is = (Igap(Isec):1:Igap(Isec+1)-1).';
    else
        Is = (Igap(Isec):1:Nobs).';
    end
    
    
    T = Data(Is,Col.T);
    M = Data(Is,Col.M);
    E = Data(Is,Col.E);
    
    Interp(Isec).MinT       = min(T);
    Interp(Isec).MaxT       = max(T);
    Interp(Isec).ValidRange = [Interp(Isec).MinT - InPar.ExtrapOutOfRange, Interp(Isec).MaxT + InPar.ExtrapOutOfRange];
    Interp(Isec).Data       = [T, M, E];
    
    % fit
    NobsSec  = numel(T);
    TimeSpan = range(T);
    TableI   = TimeSpan > InPar.PolyDegTable(:,1) & TimeSpan < InPar.PolyDegTable(:,2);
    MaxDeg   = InPar.PolyDegTable(TableI,3);
    Deg      = min(NobsSec-2,InPar.MaxDeg);
    if InPar.LogTime
        T = log10(T);
    end
    % subtract meanT from R
    MeanT    = mean(T);
    T        = T - MeanT;
    
    Par      = polyfit(T,M,Deg);
    Y        = polyval(Par,T);
    % residuals from fit
    Resid    = M - Y;
    Std      = std(Resid);

    Chi2     = sum(Resid./E).^2;
    Dof      = NobsSec - Deg;
    
    
    %Interp(Isec).DataOrig   = [T, Mmag, E];
    Interp(Isec).Par        = Par;
    Interp(Isec).MeanT      = MeanT;
    Interp(Isec).Deg        = Deg;
    Interp(Isec).Nobs       = NobsSec;
    Interp(Isec).Resid      = Resid;
    Interp(Isec).Std        = Std;
    Interp(Isec).Chi2       = Chi2;
    Interp(Isec).Dof        = Dof;
    
    tt = (Interp(Isec).ValidRange(1):InPar.TableStep:Interp(Isec).ValidRange(2)).';
    ttNL = tt;
    if InPar.LogTime
        tt = log10(tt);
    end
    tt = tt - MeanT;
    YY = polyval(Par,tt);
%     if (InPar.PolyFlux)
%         YY = -2.5.*log(YY);
%     end
    Interp(Isec).IntTable   = [ttNL, YY];
    
    if InPar.Plot
        plot.errorxy([T,Mmag,E])
        hold on;
        plot(ttNL,YY,'k-');
    end
end    
    
    
    
    
    


