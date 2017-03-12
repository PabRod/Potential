%% 1D
field = @(x) cos(x).*sin(x);
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D gradient
field = @(x) [x(2), x(1)];
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D non gradient
field = @(x) [x(2), x(1).*x(2)];
isGradient = IsGradientField(field);

assert(~isGradient);

%% 3D gradient
field = @(x) [x(2).*x(3), x(1).*x(3), x(2).*x(3)];
isGradient = IsGradientField(field);

assert(isGradient);

%% 3D non gradient
field = @(x) [x(2).*x(3), x(3), x(2).*x(3)];
isGradient = IsGradientField(field);

assert(~isGradient);
