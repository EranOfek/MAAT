function [varargout]=nlinfit_my(X,Y,Fun,varargin)
%------------------------------------------------------------------------------
% nlinfit_my function                                                  General
% Description: A version of the built in nlinfit.m function in which it
%              is possible to pass additional parameters to the function,
%              with no need for nested functions.
% Input  : - Vector of X data.
%          - Vector of Y data.
%          - A function handle to fit, Y=Fun(OptimizedPars,X).
%            Alternatively, a cell vector in which the first
%            element is the function handle, and the rest of the array
%            contains additional parameters to pass to Fun
%            {@Fun,Par1,Par2,...}.
%            In this case Y=Fun(OptimizedPars,X,Par1,Par2,...).
%          * Additional parameters as defined in nlinfit.m
%            See fminsaerch.m for details.
% Output : * The output parameters returned by nlinfit.m 
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Fun = @(P,X,A)(P(1).*X-A).^2
%          X = [1:1:100]'; Y=Fun(1.5,X,2)+randn(size(X)).*0.01;
%          [Beta,R]=nlinfit_my(X,Y,{Fun,2},0.7); % returns Beta=1.5
% Reliable: 2
%------------------------------------------------------------------------------

if (~iscell(Fun)),
   Fun = {Fun};
end

CallFun = Fun{1};
FunPar  = Fun(2:end);
 
varargout = cell(nargout,1);

[varargout{:}] = nlinfit(X,Y,@call_fun,varargin{:});


   function ModelY=call_fun(Pars,X);
   ModelY = feval(CallFun,Pars,X,FunPar{:});
   end

end
