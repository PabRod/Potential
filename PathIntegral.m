function s = PathIntegral(field, varargin)
%PATHINTEGRAL Computes the integral of a vector field along the
%given path
%
% PathIntegral(field, curve, dcurve, tmin, tmax) uses a parametric curve as
% an integration path. The derivative dcurve respective to the parameter
% should be provided, as well as the parameter limits tmin and tmax
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
        
    %TODO: add case 4 with numerical computation of dcurve
        
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

%% TODO: Input control

%% TODO: compute dcurve numerically
% epsilon = 1e-14;
% dcurvenum = @(t,s) (curve(t+s) - curve(t))/s;
% dcurve = @(t) dcurvenum(t, epsilon);

%% Compute
fieldt = @(t) field(curve(t)); % Evaluate the field on the path
integrand = @(t) sum(fieldt(t).*dcurve(t)); % Scalar integrand
s = integral(integrand, tmin, tmax, 'ArrayValue', 1);

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

end

