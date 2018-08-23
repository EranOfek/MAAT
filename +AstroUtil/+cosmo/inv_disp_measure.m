function Z=inv_disp_measure(DM,varargin)
% Convert cosmological line of sight dispersion measure to redshift
% Package: AstroUtil.cosmo
% Description: Convert cosmological line of sight dispersion measure to
%              redshift.
% Input  : - Vector of DM [cm^-3 pc].
%          - An optional structure of cosmological parameters (see
%            AstroUtil.cosmo.cosmo_pars). Alternatively a string of
%            parameters names. Default is 'planck'.
% Output : - Vector of redshifts.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Z=AstroUtil.cosmo.inv_disp_measure(1000);
% Reliable: 2
%--------------------------------------------------------------------------

Ndm = numel(DM);
Z   = zeros(Ndm,1);
for Idm=1:1:Ndm
    Z(Idm) = Util.find.fun_binsearch(@AstroUtil.cosmo.disp_measure,DM(Idm),[0 3],varargin{:});
end