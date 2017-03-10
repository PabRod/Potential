function testPathIntegral

%% Tests parameters
absTol = 1e-12;

%% Parametric cases
%% Case: straight path
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t];
dcurve = @(t) [1, 1];
tmin = -2;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 3;
assertTrue(abs(s - expected) < absTol);

%% Case: Repeat with column input
field = @(x) [-x(2); -x(1)];
curve = @(t) [t; t];
dcurve = @(t) [1; 1];
tmin = -2;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 3;
assertTrue(abs(s - expected) < absTol);

%% Case: curvilinear path
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t.^2];
dcurve = @(t) [1, 2.*t];
tmin = 0;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = -1;
assertTrue(abs(s - expected) < absTol);

%% Case: closed cycle on a conservative system
field = @(x) [-x(2).^2, -2.*x(1).*x(2)];
curve = @(t) [cos(t), sin(t)];
dcurve = @(t) [-sin(t), cos(t)];
tmin = -pi;
tmax = pi;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = 0;
assertTrue(abs(s - expected) < absTol);

%% Case: 3-dimensional problem
field = @(x) [-x(2).*x(3), -x(1).*x(3), -x(2).*x(3)];
curve = @(t) [t, t, t];
dcurve = @(t)[1, 1, 1];
tmin = 0;
tmax = 1;

s = PathIntegral(field, curve, dcurve, tmin, tmax);

expected = -1;
assertTrue(abs(s - expected) < absTol);

%% Integration over conservative field
%% Case: cartesian
field = @(x) [-2.*x(1).*x(2), -x(1).^2];
x0 = [1 2];
x = [3 4];
s = PathIntegral(field, x0, x);

expected = -34;
assertTrue(abs(s - expected) < absTol);

%% Case: origin as default
field = @(x) [-x(2), -x(1)];
x = [1 1];
s = PathIntegral(field, x);

expected = -1;
assertTrue(abs(s - expected) < absTol);

%% Case: one dimensional integral
field = @(x) x;
x0 = 0;
x = 1;
s = PathIntegral(field, x0, x);

expected = 0.5;
assertTrue(abs(s - expected) < absTol);

%% Tests of usability
%% Case: wrong number of input arguments
field = @(x) [-x(2), -x(1)];
curve = @(t) [t, t];
dcurve = @(t) [1, 1];
tmin = -2;
tmax = 1;

assertExceptionThrown(@() PathIntegral(field), 'PathIntegral:WrongNargin');
assertExceptionThrown(@() PathIntegral(field, curve, dcurve, tmin, tmax, 1, 2, 3), 'PathIntegral:WrongNargin');

end