% Tests parameters
absTol = 1e-12;

%% Case: 1D correct potential
field = @(x) -cos(x);
Vexpected = @(x) sin(x);
x = -6:0.05:6;

x0 = -3;
V0 = Vexpected(x0);

V = Potential(field, V0, x0, x);

% plot(x, V);
% hold on; 
% plot(x, Vexpected(x)); 
% legend('Numerical', 'Exact');

errors = abs(V - Vexpected(x));
assert(all(errors(:) < absTol));

%% Case: 1D incorrect potential
field = @(x) -cos(x);
Vexpected = @(x) sin(x) + x;
x = -6:0.05:6;

x0 = -3;
V0 = Vexpected(x0);

V = Potential(field, V0, x0, x);

% plot(x, V);
% hold on; 
% plot(x, Vexpected(x)); 
% legend('Numerical', 'Exact');

errors = abs(V - Vexpected(x));
assert(~all(errors(:) < absTol));

%% Case: 2D correct potential
field = @(x) [ -x(2), -x(1)];
Vexpected = @(x,y) x.*y;
x = -1:0.025:1;
y = -1:0.05:1.2;

[xm, ym] = meshgrid(x,y);
x0 = zeros(1, 2);
V0 = Vexpected(x0(1), x0(2));

V = Potential(field, V0, x0, xm, ym);

% figure; surf(xm, ym, V);

errors = abs(V - Vexpected(xm, ym));
assert(all(errors(:) < absTol));

%% Case: 2D incorrect potential
field = @(x) [ -x(2), -x(1)];
Vexpected = @(x,y) x.*y.^2;
x = -1:0.025:1;
y = -1:0.05:1.2;

[xm, ym] = meshgrid(x,y);
x0 = zeros(1, 2);
V0 = Vexpected(x0(1), x0(2));

V = Potential(field, V0, x0, xm, ym);

errors = abs(V - Vexpected(xm, ym));
assert(~all(errors(:) < absTol));

%{
%% Case: 3D correct potential
field = @(x) [ -x(2).*x(3), -x(1).*x(3), -x(1).*x(2)];
Vexpected = @(x,y,z) x.*y.*z;
x = 0:0.025:1;
y = 0:0.05:1.2;
z = 0:0.05:1;

[xm, ym, zm] = meshgrid(x,y,z);
x0 = zeros(1, 3);
V0 = Vexpected(x0(1), x0(2), x0(3));

V = Potential(field, V0, x0, xm, ym, zm);

errors = abs(V - Vexpected(xm, ym, zm));
assert(all(errors(:) < absTol));

%}

% Tests of usability
%{
%% Error in input
% field = @(x) [-x(2), -x(1)];
% 
% assertError(matlab.unittest.TestCase.forInteractiveUse, @() Potential(field), 'Potential:WrongNargin');
%}
