function [Mag,Res]=microlens_psfs(Pars,T,LD,Ngrid)
% Microlening light curve with finite source effect
% Package: AstroUtil.microlensing
% Description: Calculate the microlensing light curve for a point source.
% Input  : - Either vector of free parameters
%            [T0,Beta,V,Alpha,BaseMag, FS]
%            or a matrix of the same free parameters,
%            or a structure array with these fields and parameters.
%            T0 is time of min. impact parameters [days].
%            Beta is the impact parameter in units of the Einstein radius.
%            V is the relative velocity (e.g., Einstein radius per day).
%            Alpha is the blending parameter, 0<Alpha<=1, Alpha=1 if no blending.
%            BaseMag is the base line magnitude [mag].
%            FS is the finite source radius in units of the Einstein
%            radius.
%            Default is [0, 0.05, 0.1, 1, 18, 0.1].
%          - Vector of times. Default is (-20:0.1:20).'.
%          - Function handle for a limb darkning function.
%            I=Fun(Radius); where Radius is between 0 and 1, and I
%            is the intensity in the range of 0 and 1.
%            Default is flat function (@(r)1).
%          - Number of integration grid point along the radius.
%            Default is 10.
% Output : - Magnitude for each time.
%          - Structure array of additional parameters for each time.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mag,Res]=microlens_psfs; plot(Res.T,Mag);
% Reliable: 2
%--------------------------------------------------------------------------

Def.T     = (-20:0.1:20).';
%           T0,Beta,V,Alpha,BaseMag, FS
Def.Pars  = [0, 0.05, 0.1, 1, 18, 0.1];
Def.LD    = @(r) ones(size(r));
Def.Ngrid = 100;
if (nargin==0)
    Pars = Def.Pars;
    T    = Def.T;
    LD   = Def.LD;
    Ngrid = Def.Ngrid;
elseif (nargin==1)
    T    = Def.T;
    LD   = Def.LD;
    Ngrid = Def.Ngrid;    
elseif (nargin==2)
    LD   = Def.LD;
    Ngrid = Def.Ngrid;    
elseif (nargin==3)
    Ngrid = Def.Ngrid;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isstruct(Pars))
    In = Pars;
elseif (isnumeric(Pars))
    In.T0      = Pars(:,1);
    In.Beta    = Pars(:,2);
    In.V       = Pars(:,3);
    In.Alpha   = Pars(:,4);
    In.BaseMag = Pars(:,5);
    In.FS      = Pars(:,6);
else
    error('Unknown Pars input');
end
    


Xt = In.V.*(T-In.T0);
StepFS = In.FS./Ngrid;
RadFS = (-In.FS:StepFS:In.FS).';
[MatX,MatY] = meshgrid(RadFS,RadFS);
MatR = sqrt(MatX.^2+MatY.^2);
Lum  = LD(MatR);
%Lum(100,100)=0.5;
Lum(MatR>In.FS) = 0;
MatYp = MatY + In.Beta;

Nt = numel(T);
Res.Mu = zeros(Nt,1);
for It=1:1:Nt
   MatXp = MatX + Xt(It);
   
   U2 = MatXp.^2 + MatYp.^2 + (StepFS./2).^2; % last term to avoid singularity 
   Mu = (U2 + 2)./(sqrt(U2).*sqrt(U2+4));  % Total Magnification

   Res.Mu(It) = sum(Lum(:).*Mu(:))./sum(Lum(:));
end

BaseFlux = 1;
Res.Flux = (1-In.Alpha).*BaseFlux  + In.Alpha.*BaseFlux.*Res.Mu;
Mag = In.BaseMag-2.5.*log10(Res.Flux);
Res.T    = T;