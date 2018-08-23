function [Dist,Z]=dist_mod2dist(DM)
% Convert distance modulous to luminosity distance and redshift
% Package: AstroUtil.cosmo
% Description: Convert distance modulous to luminosity distance and
%              redshift.
% Input  : - Distance modulous [mag].
% Output : - Luminosity distance [pc].
%          - Redshift.
% See also: AstroUtil.cosmo.inv_lum_dist
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Dist,Z]=AstroUtil.cosmo.dist_mod2dist(35)
% Reliable: 2
%--------------------------------------------------------------------------


Dist = 10.^(0.2.*(DM+5));
if (nargout>1),
    Z    = AstroUtil.cosmo.inv_lum_dist(Dist,'LD');
end
