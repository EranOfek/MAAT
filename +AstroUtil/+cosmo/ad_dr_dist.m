function Dist=ad_dr_dist(Z_vecA,CosmoPars)
% Calculates the Dyer-Roeder angular diameter distance for the empty beam 
% Package: AstroUtil.cosmo
% Description: Calculates the Dyer-Roeder angular diameter distance,
%              for the empty beam case (for filled beam use: ad_dist.m),
%              between two redshifts on the line of sight.
% Input  : - Two, or one columns matrix containing redshifts [z1 z2].
%            For each row, the angular diameter distance is
%            calculated between z1 and z2. If one column is given [z2],
%            then calculate the angular diameter distance from
%            z1=0 to z2.
%          - Cosmological parameters : [H0, \Omega_{m}, \Omega_{\Lambda}],
%            or cosmological parmeters structure, or a string containing
%            parameters source name, default is 'wmap9' (see cosmo_pars.m).
% Output : - Vector of Angular diameter distance, in parsecs.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Oct 1999
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Dist=AstroUtil.cosmo.ad_dr_dist(0.6);
%------------------------------------------------------------------------------

c = 29979245800; % cm/sec

if (nargin==1),
    CosmoPars = 'wmap9';
end
Pars    = cosmo_pars(CosmoPars);
H_0     = Pars.H0;
Omega_0 = Pars.OmegaM;
Lambda_0= Pars.OmegaL;

% convert all the units to cgs:
% c is allready in cgs
% pc = 3.08568e18 cm
H_0 = H_0.*1000.*100./(3.08568e18.*1e6); % 1/sec.


SizeZ = size(Z_vecA);
Dist = zeros(SizeZ(1),1);
for J=1:1:SizeZ(1),
   Z_vec = Z_vecA(J,:);


   if (length(Z_vec(1,:))==1),
      Z1 = 0;
      Z2 = Z_vec(:,1);
   elseif (length(Z_vec(1,:))==2),
      Z1 = Z_vec(:,1);
      Z2 = Z_vec(:,2);
   else
      error('Illigal number of columns in redshift matrix');
   end
   
   
   
   
   Lh = length(H_0);
   Lo = length(Omega_0);
   Ll = length(Lambda_0);
   MaxSizeVec = max([Lh,Lo,Ll]);
   if (MaxSizeVec>1),
      if (Lh==1),
         H_0 = H_0.*ones(MaxSizeVec,1);
      end
      if (Lo==1),
         Omega_0 = Omega_0.*ones(MaxSizeVec,1);
      end
      if (Ll==1),
         Lambda_0 = Lambda_0.*ones(MaxSizeVec,1);
      end
   end
   
   
   
   % calculate the DR distance.
   R0 = c./H_0;  % cm.
   
   
   Tol = 1e-4;   % integration tolerance
   
   % integrate ChiIntegrand from z1 to z2
   
   IntVal = quad('chiint_dr',Z1,Z2,Tol,[],Omega_0,Lambda_0);
   
   
   % calculate the DR angular diameter distance
   Dist(J) = R0.*(1+Z1).*IntVal;

end

% convert distance in cm to parsecs
Dist = Dist./3.08568e18;

