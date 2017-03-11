function V = Potential(field, V0, x0, varargin)
%POTENTIAL Computes numerically the potential of a gradient field
%   Gradient fields can be computed as the gradient of a scalar field. This
%   scalar field, reversed in sign, is know as potential.
%
% Syntax:
%
% V = Potential(field, V0, x0, x, y, ...);
%
% where:
%
% FIELD is the underlying field, as an anonymous function of the vector x
% V0 is the value of the potential at the point x0
% X, Y, ... are the coordinates of the points for computing the potential,
% in the form of a meshgrid
%
% Example: one-dimensional potential
%
% field = @(x) -4*x.^3 + 3*x.^2 + 10*x - 1;
% x0 = -2;
% V0 = 0
% x = -2:0.05:3;%
% V = Potential(field, V0, x0, x);
%
% Example: two-dimensional potential
%
% field = @(x) [ -x(2), -x(1)];
% Vexpected = @(x,y) x.*y;
% x = -1:0.025:1;
% y = -1:0.05:1.2;
% 
% [xm, ym] = meshgrid(x,y);
% x0 = zeros(1, 2);
% V0 = 0
% 
% V = Potential(field, V0, x0, xm, ym);
%
% Pablo Rodr�guez-S�nchez. March 2017

%% Input parse and control
% Error control
if nargin <= 3
    msgId = 'Potential:WrongNargin';
    errMsg = 'Wrong number of input arguments';
    error(msgId, errMsg);
end

% Inputs measurements
dims = nargin - 3; % Dimensions of the problem space = number of varargin
gridSize = size(varargin{1});

% Display warning for higher dimensional potentials
if dims > 1
    warning('The current implementation of Potential.m cannot guarantee that the field is gradient for dimensions higher than one');
end

% Input vectorization
Xmesh = cell(1, dims);
Xvect = zeros(dims, numel(varargin{1}));
for i = 1:dims
    Xmesh{i} = varargin{i};
    temp = Xmesh{i};
    Xvect(i, :) = temp(:); % Each row represents a coordinate X, Y, ...
end

nPoints = size(Xvect, 2); % Number of grid points

%% Computation of the potential
% Loop
Vvect = zeros(1, nPoints); % Potential, vectorized
for i = 1:nPoints
    xCurr = Xvect(:,i)';
    Vvect(i) = V0 - PathIntegral(field, x0, xCurr);
end

% Output
V = reshape(Vvect, gridSize); % Potential, reshaped as a matrix

end