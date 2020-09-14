function [x1,x2,x3]=cubicZeros(a,b,c,d)
% this functions solves the cubic equation of the form:
% ax^3+bx^2+cx+d=0
% a,b,c,d should by either scalar or matrices of the same form.

a2=b./a;
a1=c./a;
a0=d./a;

% solve for x=y-a2/3=y-b/3a;
% the equation now becomes:
% y^3+y(3a1-a2^2)/9-(9a1a2-27a0-2a2^3)/27 = 0;

% or in the general form
% y^3+py=q
% where p&q are 
%p=a1-a2.^2/3;
%q=(9.*a1.*a2-27.*a0-2.*a2.^3)/27;

%substitute y=w-p/3w and the equation becomes:
%Cardanos method:
%(introduce v,w such as y=v+w
% v^3+w^3+(3vw+p)(v+w)=q
% require 3vw=p =>v=-p/3w which lead to y=w-p/3w )
%
% w^3-p^3/(27w^3)-q = 0;
% or (w^3)^2-q(w^3) - p^3/27 = 0
% which its solution is
% w^3 = 0.5(q+-sqrt(q^2+4/27*p^3) = 0.5q+-sqrt(0.25q^2+1/27p^3) =
% =R+-sqrt(R^2+Q^3)
% where
Q = a1./3-a2.^2/9;
R = (9.*a1.*a2-27.*a0-2.*a2.^3)/54;
D = Q.^3+R.^2;

x1 = zeros(size(D));
x2 = zeros(size(D));
x3 = zeros(size(D));

% D>0 - there is only one real solution
% x1(D>0)=nthroot(R(D>0)+sqrt(D(D>0)),3)+nthroot(R(D>0)-sqrt(D(D>0)),3)-a2(D>0)./3;
% nthroot uses one iteration of Newton's method to deal with numerical
% errors. This doubles this computation time and is not required here.
wp = R(D>0)+sqrt(D(D>0));
wm = R(D>0)-sqrt(D(D>0));
x1(D>0)=sign(wp).*abs(wp).^(1/3)+sign(wm).*abs(wm).^(1/3)-a2(D>0)./3;
x2(D>0)=nan;
x3(D>0)=nan;

% D=0 (Q^3+R^2 = 0)
% y^3+py=q <==> y^3+3Qy-2R=0
% Q=0 ==> R=0 ==> trivial case
% otherwise (Q~=0)
% (y-2R/Q)(y+R/Q)^2=(y-2R/Q)(y^2+2yR/Q+R^2/Q^2) =
%                  =y^3 -y^2*2R/Q +y^2*2R/Q -4yR^2/Q^2 +yR^2/Q^2 +2R^3/Q^3=
%                  =y^3 -3R^2/Q^2+ 2R^3/Q^3 =
% [Q^3=-R^2]       =y^3 +3yQ^3/Q^2- 2Q^3R/Q^3 = y^3+ 3Qy - 2R
%
% y^3+3Qy-2R = 0 = (y-2R/Q)(y+R/Q)^2
x1(D==0&Q==0)=-a2(D==0&Q==0)./3;
x1(D==0&Q~=0)=2.*R(D==0&Q~=0)./Q(D==0&Q~=0)-a2(D==0&Q~=0)./3;
x2(D==0&Q==0)=nan;
x2(D==0&Q~=0)=-R(D==0&Q~=0)./Q(D==0&Q~=0)-a2(D==0&Q~=0)./3;
x3(D==0)=nan;

% solution for D<0
% y=u*cos(theta)
% y^3+py-q=0=u^3cos^3(theta)+pucos(theta)-q    / *4/u^3
%          0=4cos^3(theta)+4/u^2*cos(theta)-4q/^u^3
%
% substitute u=2*sqrt(-p/3)
%          0=4cos^3(theta)+3cos(theta)-3q/2p*sqrt(-3/p)
%
% trigonometrical identity:
%          0=4cos^3(theta)+3cos(theta)-cos(3*theta)
%
%  ==>     cos(3*theta) = 3q/2p*sqrt(-3/p) = R/Q*sqrt(-Q) = sqrt(-R^2/Q^3)
%
% y=u*cos(theta) = -2*sqrt(-Q)*cos(theta)
theta = zeros(size(D));
theta(D<=0)=acos(R(D<=0)./sqrt(-Q(D<=0).^3));
sqrtmQ = 2*sqrt(-Q(D<=0));
x1(D<=0)=sqrtmQ.*cos(theta(D<=0)/3)-a2(D<=0)./3;
x2(D<=0)=sqrtmQ.*cos((theta(D<=0)+2*pi)/3)-a2(D<=0)./3;
x3(D<=0)=sqrtmQ.*cos((theta(D<=0)+4*pi)/3)-a2(D<=0)./3;

