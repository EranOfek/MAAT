function [FWHM,Center,MaxX]=get_fwhm(X,Y)
%--------------------------------------------------------------------------
% get_fwhm function                                                 FitFun
% Description: Given a 1-D vector estimate the FWHM by measuring the
%              distance between the two extreme points with height
%              of 0.5 of the maximum height.
% Input  : - Equally spaced (and asendingly sorted) X.
%          - Y vector.
% Output : - Full Width at Half Maximum.
%          - Center of FWHM.
%          - Maximum height position.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    May 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:100)'; Y=fun_gauss([1,50,5],X); [FWHM,Center,MaxX]=get_fwhm(X,Y);
% Reliable: 2
%--------------------------------------------------------------------------
if (isempty(X)),
    X = (1:1:length(Y))';
end

Y = Y./max(Y);  % normalize Y


[~,Ind] = max(Y);
MaxX    = X(Ind);

X1      = interp1(Y(1:Ind),X(1:Ind),0.5,'cubic');
X2      = interp1(Y(Ind:end),X(Ind:end),0.5,'cubic');
FWHM    = X2 - X1;
Center  = 0.5.*(X1+X2);
%ListZ   = find_local_zeros(X,Y-0.5);
%FWHM    = max(ListZ(:,1))-min(ListZ(:,1));
%Center  = 0.5.*(max(ListZ(:,1))+min(ListZ(:,1)));


