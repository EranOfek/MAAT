function NewSpec=shift_spec(Spec,Z,Method,MethodF)
%--------------------------------------------------------------------------
% shift_spec function                                            AstroSpec
% Description: Transform a spectrum from the observed frame to the rest frame.
%              If the redshift is negative then transform from the rest
%              frame to the observed frame.
% Input  : - Spectrum, in which the first column is
%            the wavelength/frequency/energy, and the second
%            is the specific flux.
%          - The redshift Z.
%            If the redshift is negative then transform from the rest
%            frame to the observed frame.
%          - Working method in: wavelength/frequency/energy
%            'w' - wavelength [lam_rest = lam_obs /(1+Z)] -default.
%            'f' - frequency.
%            'e' - energy.
%          - Specific flux conversion method:
%            'w' - f_{lambda} -> f_rest=f_obs * (1+Z)^3 - default.
%            'f' - f_{nu}     -> f_rest=f_obs * (1+Z) 
%            'wi'- f_{lambda} -> f_obs=f_rest * (1+Z)^-3
%            'fi'- f_{nu}     -> f_obs=f_rest * (1+Z)^-1
%            'n' - do nothing.
% Output : - The new shifted spectrum.
% Tested : Matlab 5.3
%     By : Eran O. Ofek / Ofer Yaron       May 2000   
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==2),
   Method  = 'w';
   MethodF = 'w';
elseif (nargin==3),
   MethodF = 'w';
elseif (nargin==4),
   % no default.
else
   error('Illegal number of input arguments');
end

switch Method
 case 'w'
    % Wavelength    \lambda_{o} =  \lambda_{e} / (1+Z)
    W_o = Spec(:,1);
    if (Z>=0),
        W_n = W_o./(1+Z);
    else
        W_n = W_o.*(1+abs(Z));
    end
    NewSpec1 = W_n;
 case 'f'
    % Frequency    f_{o} = f_{e} (1+Z)
    F_o = Spec(:,1);
    if (Z>0),
        F_n = F_o.*(1+Z);
    else
        F_n = F_o./(1+abs(Z));
    end
    NewSpec1 = F_n;
 case 'e'
    % Energy    e_{o} = e_{e} (1+Z)
    E_o = Spec(:,1);
    if (Z>=0),
        E_n = E_o.*(1+Z);
    else
        E_n = E_o./(1+abs(Z));
    end
    NewSpec1 = E_n;
 otherwise
    error('Unknown method');
end

switch MethodF
 case 'w'
     if (Z>0),
         NewSpec2 = Spec(:,2) .* (1+Z).^3;
     else
         NewSpec2 = Spec(:,2) ./ (1+abs(Z)).^3;
     end
 case 'f'
     if (Z>0),
         NewSpec2 = Spec(:,2) .* (1+Z);
     end
 case 'wi'
     if (Z>0),
         NewSpec2 = Spec(:,2) .* (1+Z).^-3;
     else
         NewSpec2 = Spec(:,2) .* (1+abs(Z)).^3;
     end
 case 'fi'
     if (Z>0),
         NewSpec2 = Spec(:,2) ./ (1+Z);
     else
         NewSpec2 = Spec(:,2) .* (1+abs(Z));
     end
 case 'n'
     NewSpec2 = Spec(:,2);      
 otherwise
    error('Unknown method');
end


NewSpec = [NewSpec1, NewSpec2];

