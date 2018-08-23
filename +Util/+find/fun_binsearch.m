function Mid=fun_binsearch(Fun,Y,Range,Tol,varargin)
%--------------------------------------------------------------------------
% fun_binsearch function                                            FitFun
% Description: Given a monotonic function, Y=Fun(X), and Y, search for
%              X that satisfy Y=F(X). The search is done using a binary
%              search between the values stated at X range.
% Input  : - Function [i.e., Y=Fun(X)].
%          - Y value to search
%          - X range [min max] in which to search for a solution.
%          - Relative tolerance, default is 1e-3.
%          - Additional optional parameters of Fun
% Output : - X value corresponds to the input Y value.
% Tested : Matlab 6.5
%     By : Eran O. Ofek                    Feb 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fun=@(x,a) a+x
%          Mid=fun_binsearch(Fun,6,[0 100],1e-6,1)
% Notes: Previously called: binsear_f.m
% Reliable: 1
%--------------------------------------------------------------------------
if (nargin==3),
   Tol = 1e-3;
else
   % do nothing
end

Mid = mean(Range);
Y1  = feval(Fun,Range(1),varargin{:});
Y2  = feval(Fun,Range(2),varargin{:});
Ym  = feval(Fun,Mid,varargin{:});

if (Y2>Y1),
   % ascending function
   Type = 'a';
else
   % descending function
   Type = 'd';
end


while (diff(Range)>(Mid.*Tol)),
   switch Type
    case 'a'
       if (Y>Ym),
          Range = [Mid, Range(2)];
       else
          Range = [Range(1), Mid];
       end
    case 'd'
       if (Y>Ym),
          Range = [Range(1), Mid];
       else
          Range = [Mid, Range(2)];
       end
    otherwise
       error('Unknwon Type Option');
   end
   Mid = mean(Range);
   Y1  = feval(Fun,Range(1),varargin{:});
   Y2  = feval(Fun,Range(2),varargin{:});
   Ym  = feval(Fun,Mid,varargin{:});
end
