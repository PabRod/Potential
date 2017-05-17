%% 1D function_handle
field = @(x) cos(x).*sin(x);
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D gradient function_handle
field = @(x) [x(2), x(1)];
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D non gradient function_handle
field = @(x) [x(2), x(1).*x(2)];
isGradient = IsGradientField(field);

assert(~isGradient);

%% 3D gradient function_handle
field = @(x) [x(2).*x(3), x(1).*x(3), x(2).*x(3)];
isGradient = IsGradientField(field);

assert(isGradient);

%% 3D non gradient function_handle
field = @(x) [x(2).*x(3), x(3), x(2).*x(3)];
isGradient = IsGradientField(field);

assert(~isGradient);

%% 1D symfun
syms x;
field(x) = cos(x).*sin(x);
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D gradient symfun
syms x y;
field(x,y) = [y, x];
isGradient = IsGradientField(field);

assert(isGradient);

%% 2D non gradient symfun
syms x y;
field(x,y) = [y, x.*y];
isGradient = IsGradientField(field);

assert(~isGradient);

%% 3D gradient symfun

syms x y z;
field(x, y, z) = [y.*z, x.*z, x.*y];
isGradient = IsGradientField(field);

assert(isGradient);

%% 3D non gradient symfun
syms x y z;
field(x, y, z) = [y.*z, z, y.*z];
isGradient = IsGradientField(field);

assert(~isGradient);

% Usability tests

%% Too many dimensions Function_handle
field = @(x) [x(2).*x(3), x(1).*x(3), x(2).*x(3), x(4)];

try
    isGradient = IsGradientField(field); % Too many dimensions
    assert(false, 'Exception failed to be thrown');
catch me
    expectedError = 'IsGradientField:WrongDimension';
    assert(strcmp(me.identifier, expectedError));
end

%% Too many dimensions SymFun
syms x y z u;
field(x, y, z, u) = [y.*z, z, y.*z, u];

try
    isGradient = IsGradientField(field); % Too many dimensions
    assert(false, 'Exception failed to be thrown');
catch me
    expectedError = 'IsGradientField:WrongDimension';
    assert(strcmp(me.identifier, expectedError));
end