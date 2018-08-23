function [FAP]=fab_counts(C,B,Ratio)
%-------------------------------------------------------------------------
% fab_counts function                                           AstroStat                              
% Description: Calculate the False Alarm Probability (FAP) that a source
%              is real rather than a background fluctutation,
%              given the source counts and the background counts. 
% Input  : - Source counts.
%          - Expectency of background counts within the source
%            extraction aperture.
%            Alternatively, if a third argument is given this is the
%            total measured background counts within a given area.
%            In this case the third parameter give the ratio between
%            the background area and the source extraction area.
%          - The ratio between the background area and the source
%            extraction area.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Oct 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%-------------------------------------------------------------------------


if (nargin==2),
    N = length(C);
    % calculate the probability that the observed count rate is
    % due to background noise
    FAP = zeros(N,1).*NaN;
    for I=1:1:N,
       FAP(I) = 1 - poisscdf(C(I)-1,B(I));
    end
elseif (nargin==3),
    
    Nsim = 10000;
    SimFAP = zeros(Nsim,1);
    RP_b = poissrnd(B,Nsim,1);   % random poisson counts for background
    SimFAP = 1 - poisscdf(C-1,RP_b./Ratio);
    FAP = mean(SimFAP);
end


