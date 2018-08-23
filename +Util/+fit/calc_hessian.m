function [H]=calc_hessian(Fun,MinX,DeltaX,varargin)
% Calculate the Hessian matrix of a multivariable function.
% Package: Util.fit
% Description: Calculate the Hessian (second derivative) matrix of a
%              multivariable function.
% Input  : - Function: y=fun(X,varargin)
%          - Vector of parameters X around which to calculate the Hessian.
%          - Vector of DX to use in the differentiation process,
%            default (in case empty) is MinX.*1e-2.
%          - Additional parameters to pass to Fun.
% Output : - Hessian matrix
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Notes  : The covariance matrix is: inv(0.5.*H)
%          The errors in the parameters are: sqrt(diag(inv(0.5.*H)))
% Reliable: 2
%--------------------------------------------------------------------------
RelDiff  = 1e-2;
if (nargin==2)
   DeltaX = [];
end

if (isempty(DeltaX)==1)
   DeltaX    = MinX.*RelDiff;
   I         = find(abs(DeltaX)<RelDiff.*1e-3);
   DeltaX(I) = RelDiff;
end


FunMin = feval(Fun,MinX,varargin{:});
Nvar = length(MinX);
H    = zeros(Nvar,Nvar);
for I=1:1:Nvar
   for J=1:1:Nvar

      CurrentX = MinX;
      if (I==J)
         MinXn     = MinX;
         MinXn(I)  = MinX(I) - DeltaX(I);
         MinXp     = MinX;
         MinXp(I)  = MinX(I) + DeltaX(I);
         FunMin_n  = feval(Fun,MinXn,varargin{:});
         FunMin_p  = feval(Fun,MinXp,varargin{:});

         H(I,J)    = (-2.*FunMin + FunMin_n + FunMin_p)./(DeltaX(I).*DeltaX(J));
      else
         % -1,+1
         CurX      = MinX;
         CurX(J)   = CurX(J) - DeltaX(J);    % x coo
         CurX(I)   = CurX(I) + DeltaX(I);    % y coo
         FunMin_np = feval(Fun,CurX,varargin{:});
         % +1,-1
         CurX      = MinX;
         CurX(J)   = CurX(J) + DeltaX(J);    % x coo
         CurX(I)   = CurX(I) - DeltaX(I);    % y coo
         FunMin_pn = feval(Fun,CurX,varargin{:});
         % -1,-1
         CurX      = MinX;
         CurX(J)   = CurX(J) - DeltaX(J);    % x coo
         CurX(I)   = CurX(I) - DeltaX(I);    % y coo
         FunMin_nn = feval(Fun,CurX,varargin{:});
         % +1,+1
         CurX      = MinX;
         CurX(J)   = CurX(J) + DeltaX(J);    % x coo
         CurX(I)   = CurX(I) + DeltaX(I);    % y coo
         FunMin_pp = feval(Fun,CurX,varargin{:});
        
         %H(I,J)    = 0.25.*(2.*FunMin_nn - FunMin_np - FunMin_pn)./(DeltaX(I).*DeltaX(J));
         H(I,J)    = 0.25.*(FunMin_pp + FunMin_nn - FunMin_np - FunMin_pn)./(DeltaX(I).*DeltaX(J));
      end
   end
end


