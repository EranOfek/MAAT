function Obs=obstruction(D,R1,R2,Nstep,LimbFun,Pars)
% Stellar obstruction due to the eclipse
% Package: AstroUtil.binary
% Description: Calculate stellar obstruction due to the eclipse given the
%              stars radii, distance and limb darkening function of the
%              background star.
% Input  : - Vector of distances between primary and secondary
%            [consistent length unit].
%          - Radius of background star [length unit].
%          - Radius of forground star [length unit].
%          - Number of steps in integration.
%          - limb-darkening function, returning the luminosity per unit area
%            as function of radius.
%          - Cell array of additional optional parameters for the
%            LimbFun function. Default is {'Milne',1}.
% Output : - Vector of total obstraction of background star, for each distance
%            in distance vector [background-star luminosity].
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Aug 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Obs=AstroUtil.binary.obstruction(D,R1,R2,100,'limb_darkening',{'Kling',[1 0.5]});
% Reliable: 2
%---------------------------------------------------------------------------
import AstroUtil.EB.*

Def.Nstep   = 100;
Def.LimbFun = 'limb_darkening';
Def.Pars    = {'Milne',1};

if (nargin==3)
   Nstep   = Def.Nstep;
   LimbFun = Def.LimbFun;
   Pars    = Def.Pars;
elseif (nargin==4)
   LimbFun = Def.LimbFun;
   Pars    = Def.Pars;
elseif (nargin==5)
   Pars    = Def.Pars;
elseif (nargin==6)
   % do nothing
else
   error('Illigal number of input arguments');
end

N = length(D);

DelR = R1./Nstep;
Obs = zeros(size(D));
for I=1:1:N
   % for each distance in distances-vector
   if ((D(I)+R2)<R1)
      % total eclipse
      LimA = D(I) + R2;
      LimB = D(I) - R2;
   elseif (D(I)>(R1+R2))
      % no eclipse
      LimA = 0;
      LimB = 0;
   else
      % partial eclipse
      LimA = R1;
      LimB = D(I) - R2;
   end
   
   if (LimB<0)
      LimB = 0;
   end
   
   R      = [LimB:DelR:LimA].';
   if (length(R)>1)
      %R      = [LimB:DelR:LimA].';
      if ((R2>R1) && D(I)<(R2-R1))
         Obs(I) = total_light(R1,LimbFun,Nstep,Pars);
      else
         %L      = eval([LimbFun,'(R./R1,Pars)']);
         L      = feval(LimbFun,R./R1,Pars{:});
         I0     = find(D(I)==0 | R==0);
         In0    = find(D(I)~=0 & R~=0);
         Alpha  = zeros(size(R));
         Alpha(In0)  = acos((D(I).^2 + R(In0).^2 - R2.^2)./(2.*D(I).*R(In0)));
         % set special cases
         Ir0 = find(R==0);    % on center
         Alpha(Ir0) = pi;
         Iic = find((R+D(I))<R2);   % inner circle
         Alpha(Iic) = pi;
         %Ibc = find(R2>R1 & R<(R2-D(I))); % inner inverse size
         %Alpha(Ibc)
         %[R, abs(Alpha)]
         Alpha = abs(Alpha);
      
         Obs(I) = trapz(R,2.*R.*Alpha.*L);
      end
   else
      Obs(I) = 0;
   end
end
