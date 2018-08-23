function [OmL_EF,OmL_NS]=omega_m_lambda_lines(OmegaM)
% Selected lines in Omega_m-Omega_lambda diagram
% Package: AstroUtil.cosmo
% Description: Given a universe with \Omega_{m} and \Omega_{\Lambda}
%              contributions, and given \Omega_{m} vector, find for
%              each value of \Omega_{m}: (i) the value of
%              \Omega_{\Lambda} for which the universe will expand
%              forever; (ii) The \Omega_{\Lambda} criterion for which
%              there have been no singularity in the past (rather than
%              Big Bang its early history consisted of a period of
%              gradually slowing contraction to a minimum radius before
%              begining its current expansion).
% Input  : - Vector of \Omega_{m}.
% Output : - Vector of lower(>=) \Omega_{\Lambda} (for each \Omega_{m}),
%            for which the universe will expand forever.
%            (in the \Omega_{\Lambad} vs. \Omega_{m} plane, the
%            big-crunch, is to the left of this critical line).
%          - Vector of lower(>=) \Omega_{\Lambda} (for each \Omega_{m}),
%            for which there was no singularity in the past.
%            (in the \Omega_{\Lambad} vs. \Omega_{m} plane, the
%            no-big-bang, is to the right of this critical line).
% Tested : Matlab 5.1
%     By : Eran O. Ofek                    Jul 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% References : Carroll, S.M. 2000, astro-ph/0004075
% Example: OmegaM=[0:0.01:3]';  % define a vector of Omega matter
%          [OmL_EF,OmL_NS]=AstroUtil.cosmo.omega_m_lambda_lines(OmegaM);
%          plot(OmegaM,OmL_EF); hold on; plot(OmegaM,OmL_NS,'r');
%          xlabel('\Omega_{m}'); ylabel('\Omega_{\Lambda}');
%          axis([0 3 -1 4.5]);
%          text(1.5,-0.2,'Expand forever'); text(0.5,3.5,'No singularity');
% Reliable: 2
%------------------------------------------------------------------------------

if (nargout>0),
   OmL_EF = zeros(size(OmegaM));
   I0 = find(OmegaM<=1);
   I1 = find(OmegaM>1);
   OmL_EF(I0) = zeros(size(I0));
   OmL_EF(I1) = 4.*OmegaM(I1).*cos( acos( (1-OmegaM(I1))./(OmegaM(I1)))./3 + 4.*pi./3).^3;
end

if (nargout>1),
   OmL_NS = zeros(size(OmegaM));
   I0 = find(OmegaM<=0.5);
   I1 = find(OmegaM>0.5);
   OmL_NS(I0) = 4.*OmegaM(I0).*cosh( acosh((1-OmegaM(I0))./(OmegaM(I0)) )./3).^3;
   OmL_NS(I1) = 4.*OmegaM(I1).*cos( acos((1-OmegaM(I1))./(OmegaM(I1)) )./3).^3;
end
