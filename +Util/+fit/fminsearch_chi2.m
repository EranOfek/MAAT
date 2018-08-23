function [varargout]=fminsearch_chi2(X,Y,Err,Fun,varargin)
%------------------------------------------------------------------------------
% fminsearch_chi2 function                                              FitFun
% Description: \chi^2 fitting using fminsearch.m and given a model function.
% Input  : - Vector of X.
%          - Vector of Y.
%          - A model function handle to fit, Y=Fun(OptimizedPars,X).
%            Alternatively, a cell vector in which the first
%            element is the function handle, and the rest of the array
%            contains additional parameters to pass to Fun
%            {@Fun,Par1,Par2,...}.
%            In this case Y=Fun(OptimizedPars,X,Par1,Par2,...).
%          * Additional parameters as defined in fminsearch.m
%            See fminsaerch.m for details.
%            E.g., the first argument here is the guess parameters.
% Output : * The output parameters returned by fminsearch.m 
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fun = @(P,X,A)(P(1).*X+P(2).*sin(X)-A).^2
%          X = [1:1:100]'; Y=Fun([1.5,5],X,2)+randn(size(X)).*0.01;
%          [Beta,Chi2,E,O]=fminsearch_chi2(X,Y,0.01,{Fun,2},[0.7,1]); % returns Beta=[1.5,5]
% Reliable: 2
%------------------------------------------------------------------------------

if (~iscell(Fun)),
   Fun = {Fun};
end

if (numel(Err)==1),
   Err = Err.*ones(size(X));
end

CallFun = Fun{1};
FunPar  = Fun(2:end);
 
varargout = cell(nargout,1);

[varargout{:}] = fminsearch(@call_fun,varargin{:});


   function Chi2=call_fun(Pars);
   ModelY = feval(CallFun,Pars,X,FunPar{:});
   Chi2 = sum(((Y-ModelY)./Err).^2);
   end

end
