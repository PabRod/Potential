function testPotential

%% Tests parameters
absTol = 1e-12;

%% Parametric cases
%% Case: one-dimensional
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

assertVectorsAlmostEqual(V, Vexpected(x), 'absolute', absTol);


%% Case: two-dimensional
field = @(x) [ -x(2), -x(1)];
Vexpected = @(x,y) x.*y;
x = -1:0.025:1;
y = -1:0.05:1.2;

[xm, ym] = meshgrid(x,y);
x0 = zeros(1, 2);
V0 = Vexpected(x0(1), x0(2));

V = Potential(field, V0, x0, xm, ym);

% figure; surf(V);

assertVectorsAlmostEqual(V, Vexpected(xm,ym), 'absolute', absTol);

%% Case: three-dimensional
%{
field = @(x) [ -x(2).*x(3), -x(1).*x(3), -x(1).*x(2)];
Vexpected = @(x,y,z) x.*y.*z;
x = 0:0.025:1;
y = 0:0.05:1.2;
z = 0:0.05:1;

[xm, ym, zm] = meshgrid(x,y,z);
x0 = zeros(1, 3);
V0 = Vexpected(x0(1), x0(2), x0(3));

V = Potential(field, V0, x0, xm, ym, zm);

assertVectorsAlmostEqual(V, Vexpected(xm,ym,zm), 'absolute', absTol);

%}

%% Tests of usability
%% Case: wrong number of input arguments
field = @(x) [-x(2), -x(1)];

assertExceptionThrown(@() Potential(field), 'Potential:WrongNargin');

end
