function pass = test_poisson( pref )
% Test Poisson solver in Chebfun2

if ( nargin == 0 )
    pref = chebfunpref;
end

tol = 1e6*pref.cheb2Prefs.chebfun2eps;

% One-argument syntax (RHS):
v = chebfun2( @(x,y) (1-x.^2).*(1-y.^2).*cos(2*pi*x.*y) );
f = lap(v);
u = chebfun2.poisson(f);
pass(1) = ( norm(u-v) < tol );

% Two-argument syntax (RHS + Dirichlet data):
v = chebfun2( @(x,y) cos(2*pi*x.*y) );
f = lap(v);
u = chebfun2.poisson(f, v);
pass(2) = ( norm(u-v) < tol );

% Full syntax:
m = 101; n = 101;
v = chebfun2( @(x,y) (1-x.^2).*(1-y.^2).*cos(2*pi*x.*y) );
f = lap(v);
u = chebfun2.poisson(f, 0, m, n);
pass(3) = ( norm(u-v) < tol );

% Rectangular discretization sizes:
m = 87; n = 101;
v = chebfun2( @(x,y) (1-x.^2).*(1-y.^2).*cos(2*pi*x.*y) );
f = lap(v);
u = chebfun2.poisson(f, 0, m, n);
pass(4) = ( norm(u-v) < tol );

% Nonstandard rectangular domain:
m = 99; n = 101;
a = -3; b = 1; c = 4; d = 4.2;
v = chebfun2( @(x,y) (x-a).*(x-b).*(y-c).*(y-d).*cos(2*pi*x.*y), [a b c d] );
f = lap(v);
u = chebfun2.poisson(f, 0, m, n);
pass(5) = ( norm(u-v) < tol );

% Nonzero, but constant, Dirichlet data:
m = 99; n = 101;
a = -3; b = 1; c = 4; d = 4.2;
v = chebfun2( @(x,y) (x-a).*(x-b).*(y-c).*(y-d).*cos(2*pi*(x+y))+1, [a b c d] );
f = lap(v);
u = chebfun2.poisson(f, 1, m, n);
pass(6) = ( norm(u-v) < tol );

% General Dirichlet data, given as function handle:
m = 201; n = 201;
a = -3; b = 1; c = 4; d = 4.2;
p = @(x,y) x.*y + cos(3*x.^2.*(y-.2));
v = chebfun2( p, [a b c d] );
f = lap(v);
u = chebfun2.poisson(f, p, m, n);
pass(7) = ( norm(u-v) < tol );

% General Dirichlet data, given as chebfun2:
m = 201; n = 201;
a = -3; b = 1; c = 4; d = 4.2;
p = @(x,y) x.*y + cos(3*x.^2.*(y-.2));
v = chebfun2( p, [a b c d]);
f = lap(v);
u = chebfun2.poisson(f, v, m, n);
pass(8) = ( norm(u-v) < tol );

end
