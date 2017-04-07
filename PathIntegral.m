function s = PathIntegral(field, varargin)
%PATHINTEGRAL Computes the integral of a vector field along the
%given path
%
% PathIntegral(field, curve, dcurve, tmin, tmax) uses a parametric curve as
% an integration path. The derivative dcurve respective to the parameter
% should be provided, as well as the parameter limits tmin and tmax
%
% PathIntegral(field, curve, tmin, tmax), where field and curve are
% symfun objects perform the integral without the requirement of providing
% the derivative
%
% PathIntegral(field, x0, x) uses a straigth line between x0 and x as
% integration path. This is specially useful for working with conservative
% fields, where the path integral doesn't depend on the particular shape of
% the path, but only on the starting and ending points
%
% PathIntegral(field, x) is a shorthand for PathIntegral(field, x0, x) with
% x0 being an all-zeros vector
%
% Example 1:
%
% field = @(x) [-x(2), -x(1)];
% curve = @(t) [t, t];
% dcurve = @(t) [1, 1];
% tmin = 0;
% tmax = 1;
%
% s = PathIntegral(field, curve, dcurve, tmin, tmax);
%
% Example 2:
%
% field = @(x) [-x(2), -x(1)];
% x0 = [0 0];
% x = [1 2];
%
% s = PathIntegral(field, x0, x);
%
% Pablo Rodríguez-Sánchez. March 2017

%% Parse input and method choice
switch nargin
    
    case 2 % Integrate from point 0 to point x
        % Useful only for conservative fields where the integral depends only
        % in the initial and final points of the path
        
        x = varargin{1};
        x0 = zeros(size(x)); % Use origin as default initial point
        
        s = PathIntegralCartesian(field, x0, x);
        
    case 3 % Integrate from point x0 to point x
        % Useful only for conservative fields where the integral depends only
        % in the initial and final points of the path
        
        x0 = varargin{1};
        x = varargin{2};
        
        s = PathIntegralCartesian(field, x0, x);
        
    case 4 % Integrate along parametric curve using symbolic functions
        curve = varargin{1};
        tmin = varargin{2};
        tmax = varargin{3};
        
        s = PathIntegralSymbolic(field, curve, tmin, tmax);
        
    case 5 % Integrate along parametric curve
        % Proper and general definition of a path integral
        curve = varargin{1};
        dcurve = varargin{2};
        tmin = varargin{3};
        tmax = varargin{4};
        
        s = PathIntegralParametric(field, curve, dcurve, tmin, tmax);
        
    otherwise
        msgId = 'PathIntegral:WrongNargin';
        errMsg = 'Wrong number of input arguments';
        error(msgId, errMsg);
        
end

end

function s = PathIntegralParametric(field, curve, dcurve, tmin, tmax)
%PATHINTEGRALPARAMETRIC Computes the integral of a vector field along the
%given path
% 
% Example:
% field = @(x) [-x(2), -x(1)];
% curve = @(t) [t, t];
% dcurve = @(t) [1, 1];
% tmin = 0;
% tmax = 1;
% PathIntegral(field, curve, dcurve, tmin, tmax)

%% Compute
fieldt = @(t) field(curve(t)); % Evaluate the field on the path
integrand = @(t) sum(fieldt(t).*dcurve(t)); % Scalar integrand
s = integral(integrand, tmin, tmax, 'ArrayValue', 1);

%% TODO: Input control

%% TODO: compute dcurve numerically
% epsilon = 1e-14;
% dcurvenum = @(t,s) (curve(t+s) - curve(t))/s;
% dcurve = @(t) dcurvenum(t, epsilon);

end

function s = PathIntegralCartesian(field, x0, x)
%PATHINTEGRALCARTESIAN Computes the integral of a vector field along the
%straigth line joining two points. Specially useful for conservative
%fields, where the path integral doesn't depend on the particular path
%chosen
%
% Example:
% field = @(x) [-x(2), -x(1)];
% x0 = [0 0];
% x = [1 2];
%
% s = PathIntegralCartesian(field, x0, x);

%% Create parametric integration path
% The simplest curve joining x0 and x is a segment of a straight line
curve = @(t) x0 + (x - x0).*t;
dcurve = @(t) x - x0;
tmin = 0;
tmax = 1;

%% Call parametric method
s = PathIntegralParametric(field, curve, dcurve, tmin, tmax);

%% TODO: throw warning if field is non-gradient

end

function s = PathIntegralSymbolic(field, curve, tmin, tmax)
%PATHINTEGRALSYMBOLIC Computes the integral of a vector field along the
%provided curve. Both the field and the curve are provided as symbolic
%objects
%
% syms x y t;
% field(x,y) = [-y, -x];
% curve(t) = [t, t.^2];
% tmin = 0;
% tmax = 1;
% s = PathIntegral(field, curve, tmin, tmax));

%% Input control
arginOk = isa(field, 'symfun') && isa(curve, 'symfun');
if ~arginOk
    msgId = 'PathIntegral:SymArginRequired';
    errMsg = 'Arguments field and curve are expected to by of type symfun';
    error(msgId, errMsg);
end

%% Algorithm
syms t;

% This trick is required for adding curve(t) as the input to f(x, y, ...)
dims = numel(curve(t));
input = cell(1, dims);
for i = 1:dims
   versor = zeros(1, dims);
   versor(i) = 1;
   input{i} = curve(t)*versor'; % Example: input = {curve(t)*[1 0]', curve(t)*[0 1]'};
end

fieldt(t) = field(input{:});
dcurve(t) = diff(curve(t), t);
integrand(t) = fieldt(t)*dcurve(t)';

s = int(integrand(t), tmin, tmax);

end