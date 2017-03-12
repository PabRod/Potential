function isGradient = IsGradientField(field)
%ISGRADIENTFIELD Returns true if the introduced field is gradient, false if
%not.
%
%Warning: may throw false positives, but all negatives are robust. That
%is, if the algorithm says the field is not gradient, then surely it's
%not.
%
% Pablo Rodríguez-Sánchez, March 2017

%% Determine the input size
% Surely this can be done in a cleaner way
try
    x = 0;
    field(x);
    dims = 1;
catch
    try
        x = [0, 0];
        field(x);
        dims = 2;
    catch
        try
            x = [0, 0, 0];
            field(x);
            dims = 3;
        catch
            error('Only systems up to 3D are supported');
        end
    end
    
end

%% Check using an integration along a random closed ellipse as a test
% Any integral of a gradient field along any closed curve should equal 0.
% Here, we only test a small subset of all possible curves. That is why the
% method can throw false positives
absTol = 1e-9;
switch dims
    case 1 % One-dimensional fields are always gradient
        isGradient = true;
        return;
    case 2
        % Create a parametric ellipse
        center = 10*rand(1, dims);
        rads = 5*rand(1, dims);
        
        curve = @(t) center + rads.*[cos(t), sin(t)];
        dcurve = @(t) rads.*[-sin(t), cos(t)];
        tmin = -pi;
        tmax = pi;
        
        % Integrate along it
        closedIntegral = PathIntegral(field, curve, dcurve, tmin, tmax);

        % Check if it adds up to zero
        isGradient = (abs(closedIntegral) < absTol);
        return;
    case 3
        % Create a parametric ellipse
        center = 10*rand(1, dims);
        rads = 5*rand(1, dims);
        
        curve = @(t) center + rads.*[cos(t), sin(t), 0];
        dcurve = @(t) rads.*[-sin(t), cos(t), 0];
        tmin = -pi;
        tmax = pi;
        
        % Integrate along it
        closedIntegral = PathIntegral(field, curve, dcurve, tmin, tmax);

        % Check if it adds up to zero
        isGradient = (abs(closedIntegral) < absTol);
        return;
end