% Test file for chebtech1 constructor.
% Here, we check populate().  (This function is not user-facing.)

function pass = test_constructor(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Initialize with default data:
data = chebtech.parseDataInputs(struct());

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'nested';
f = @(x) sin(x);
g = populate(chebtech1, f, data, pref);
x = chebtech1.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(1) = norm(f(x) - values, inf) < 10*vscale(g)*eps;

% Test on an array-valued function:
pref.refinementFunction = 'nested';
f = @(x) [sin(x) cos(x) exp(x)];
g = populate(chebtech1, f, data, pref);
x = chebtech1.chebpts(length(g.coeffs));
values = g.coeffs2vals(g.coeffs);
pass(2) = norm(f(x) - values, inf) < 10*max(vscale(g)*eps);

%%
% Test on a scalar-valued function:
pref.refinementFunction = 'resampling';
f = @(x) sin(x);
[g, values] = populate(chebtech1, f, data, pref);
x = chebtech1.chebpts(length(values));
pass(3) = norm(f(x) - values, inf) < 10*vscale(g)*eps;

% Test on an array-valued function:
pref.refinementFunction = 'resampling';
f = @(x) [sin(x) cos(x) exp(x)];
[g, values] = populate(chebtech1, f, data, pref);
x = chebtech1.chebpts(length(values));
pass(4) = norm(f(x) - values, inf) < 10*max(vscale(g)*eps);

%%
% Some other tests:

% This should fail with an error:
try
    f = @(x) x + NaN;
    populate(chebtech1, f, data, pref);
    pass(5) = false;
catch ME
    pass(5) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% As should this:
try
    f = @(x) x + Inf;
    populate(chebtech1, f, data, pref);
    pass(6) = false;
catch ME
    pass(6) = strcmp(ME.message, 'Too many NaNs/Infs to handle.');
end

% Check that things don't crash if pref.minSamples and pref.maxLength are equal.
try
    pref.minSamples = 8;
    pref.maxLength = 8;
    populate(chebtech1, @sin, data, pref);
    pass(7) = true;
catch
    pass(7) = false;
end

% Test logical-valued functions:
f = chebtech1(@(x) x > -2);
g = chebtech1(1);
pass(8) = normest(f - g) < eps;

f = chebtech1(@(x) x < -2);
g = chebtech1(0);
pass(9) = normest(f - g) < eps;

end
