function Val=fermi_fun(X,X0,DX,Decay)
%--------------------------------------------------------------------------
% fermi_fun function                                               General
% Description: Calculate the Fermi function (Fermi-Dirac statistics)
%              of the form 1/(1+exp((x-x0)/DX)).
% Input  : - X.
%          - X0.
%          - DX.
%          - Decay (true) or rising (false) function. Default is true.
% Output : - Value.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: plot(fermi_fun((0:1:100),50,2))
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==3),
    Decay = true;
end
if (Decay)
    Sign = 1;
else
    Sign = -1;
end

Val = 1./(1+exp(Sign.*(X-X0)./DX));


    
