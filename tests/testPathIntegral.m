% Tests parameters
absTol = 1e-12;

% Parametric cases
%% Straight parametric path
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t];
dcurve = @(t) [1, 1];
tmin = -2;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 3;
assert(abs(s - expected) < absTol);

%% Column shaped input
field = @(x) [-x(2); -x(1)];
curve = @(t) [t; t];
dcurve = @(t) [1; 1];
tmin = -2;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 3;
assert(abs(s - expected) < absTol);

%% Curvilinear parametric path
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t.^2];
dcurve = @(t) [1, 2.*t];
tmin = 0;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = -1;
assert(abs(s - expected) < absTol);

%% Closed cycle on gradient field
field = @(x) [-x(2).^2, -2.*x(1).*x(2)];
curve = @(t) [cos(t), sin(t)];
dcurve = @(t) [-sin(t), cos(t)];
tmin = -pi;
tmax = pi;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 0;
assert(abs(s - expected) < absTol);

%% 3D field
field = @(x) [-x(2).*x(3), -x(1).*x(3), -x(2).*x(3)];
curve = @(t) [t, t, t];
dcurve = @(t)[1, 1, 1];
tmin = 0;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = -1;
assert(abs(s - expected) < absTol);

% Integration over gradient field
%% Cartesian
field = @(x) [-2.*x(1).*x(2), -x(1).^2];
x0 = [1 2];
x = [3 4];
s = PathIntegral(field, x0, x);

expected = -34;
assert(abs(s - expected) < absTol);

%% Origin as default
field = @(x) [-x(2), -x(1)];
x = [1 1];
s = PathIntegral(field, x);

expected = -1;
assert(abs(s - expected) < absTol);

%% 1D integral
field = @(x) x;
x0 = 0;
x = 1;
s = PathIntegral(field, x0, x);

expected = 0.5;
assert(abs(s - expected) < absTol);

%% 1D symbolic
syms x t;
field(x) = x;
curve(t) = t;
tmin = 0;
tmax = 1;
s = double(PathIntegral(field, curve, tmin, tmax));

expected = 0.5;
assert(abs(s - expected) < absTol);

%% 2D symbolic
syms x y t;
field(x,y) = [-y, -x];
curve(t) = [t, t.^2];
tmin = 0;
tmax = 1;
s = double(PathIntegral(field, curve, tmin, tmax));

expected = -1;
assert(abs(s - expected) < absTol);

% Tests of usability
%% Wrong nargin
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t];
dcurve = @(t) [1, 1];
tmin = -2;
tmax = 1;

try
    PathIntegral(field); % Too few input arguments
    assert(false, 'Exception failed to be thrown');
catch me
    expectedError = 'PathIntegral:WrongNargin';
    assert(strcmp(me.identifier, expectedError));
end

try
    PathIntegral(field, curve, dcurve, tmin, tmax, 1, 2, 3); % Too many input arguments
    assert(false, 'Exception failed to be thrown');
catch me
    expectedError = 'PathIntegral:WrongNargin';
    assert(strcmp(me.identifier, expectedError));
end

%% Wrong symbolic
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t.^2];
tmin = 0;
tmax = 1;

try
    PathIntegral(field, curve, tmin, tmax); % 4 input arguments requires symbolic input
    assert(false, 'Exception failed to be thrown');
catch me
    expectedError = 'PathIntegralSymbolic:SymArginRequired';
    assert(strcmp(me.identifier, expectedError));
end