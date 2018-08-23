function [D1,D2,D3]=numerical_diff(X,dt)
% Numerical differentiation of an equally spaced a vector.
% Package: Util
% Description: Numerical differentiation of an equally spaced a vector.
% Input  : - A vector representing an equally spaced evaluated function.
%          - function evaluation step. Default is 1.
% Output : - Column vector of first derivative.
%            The vector length is identical to the input vector, and the
%            derivative are evaluated at the original points.
%          - Column vector of second derivative (same length as input).
%          - Column vector of third derivative (same length as input).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    May 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [D1,D2,D3]=Util.fun.numerical_diff(sin((0:0.01:1)),0.01)
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    dt = 1;
end

X = X(:);

% first derivative
D1a = (X(2:end-1) - X(1:end-2))./dt;
D1b = (X(3:end)   - X(2:end-1))./dt;
D1  = [NaN; 0.5.*(D1a + D1b); NaN];

if (nargout>1)
    D2  = [NaN; (D1b - D1a)./dt; NaN];

    if (nargout>2)
        D3a = (X(4:end-1) - 3.*X(3:end-2) + 3.*X(2:end-3) - X(1:end-4))./(dt.^3);
        D3b = (X(5:end)   - 3.*X(4:end-1) + 3.*X(3:end-2) - X(2:end-3))./(dt.^3);
        D3  = [NaN; NaN; 0.5.*(D3a + D3b); NaN; NaN];
    end
end


