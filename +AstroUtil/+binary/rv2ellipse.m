function [RV_Even,RV_Odd,Phase,Time]=rv2ellipse(Time,RV,P)
% Convert radial velocity as a function of time to ellipse
% Package: AstroUtil.binary
% Description: Given radial velocity as a function of time and a trial
%              period convert it to a RV_Odd vs RV_Even (ellipse).
%              I.e., for a good period plot(RV_Even,RV_Odd,'.') should be
%              an ellipse.
% Input  : - Vector of time.
%          - Vector of radial velocity
%          - Trial period
% Output : - RV_Even
%          - RV_Odd
%          - Corresponding phase
%          - Corresponding time
% Reference: Bhattacharyya & Nityanada (2008; MNRAS 387, 273-278)
%     By : Eran O. Ofek            Dec 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2



% calculate phase [fraction]
Nt = numel(Time);

Phase = mod(Time,P)./P;

% add eps to phase to make sure its monotonic...
RV1 = [Phase+(eps:eps:Nt.*eps)', RV, Time];
RV1 = sortrows(RV1,1);
Phase = RV1(:,1);
Time  = RV1(:,3);

% interpolate
RV2 = interp1(RV1(:,1),RV1(:,2),1-Phase);
% remove extrapolated points
FlagOK = ~isnan(RV2);
RV1 = RV1(FlagOK,:);
RV2 = RV2(FlagOK);
Phase = Phase(FlagOK);
Time  = Time(FlagOK);


RV_Even = (RV1(:,2) + RV2).*0.5;
RV_Odd  = (RV1(:,2) - RV2).*0.5;
    