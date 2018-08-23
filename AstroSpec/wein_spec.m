function [I_nu,I_lam]=wein_spec(T,Lambda,R)


%--------------------------------------------------------------------------
% wein_spec function                                             AstroSpec
% Description: 
% Input  : - 
% Output : - 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

Lambda = Lambda.*1e-8;   % Ang->cm
sigmaB = get_constant('sigma');
c      = get_constant('c');
h      = get_constant('h');
kB     = get_constant('kB');

Nu     = c./Lambda;

I_nu   = 2./3 .*pi .*sigmaB .*R.^2 .*Nu.^4 .*exp(-h.*Nu./(kB.*T));
I_lam  = convert.flux(I_nu,'cgs/Hz','cgs/A',Nu,'Hz');
