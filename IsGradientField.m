function isGradient = IsGradientField(field)
%ISGRADIENTFIELD Returns true if the introduced field is gradient, false if
%not.
%
%Warning: may throw false positives, but all negatives are robust. That
%is, if the algorithm says the field is not gradient, then surely it's
%not.
%
% Pablo Rodríguez-Sánchez, March 2017

switch class(field)
    case {'symfun', 'sym'}
        % The algorithm is dimension-dependent...
        names = argnames(field);
        dimension = numel(names(:));
        
        % ... It is based in the fact that the skew part of the jacobian is 
        % zero everywhere for gradient fields
        
        switch dimension
            case 1
                % The skew part of a 1-D jacobian is always zero
                Jskew = 0;
                
            case 2
                J = jacobian(field);
                testPoint = rand(1, 2);
                Jskew = SkewSymmDecomposition(eval(J(testPoint(1), testPoint(2))));
                
            case 3
                J = jacobian(field);
                testPoint = rand(1, 3);
                Jskew = SkewSymmDecomposition(eval(J(testPoint(1), testPoint(2), testPoint(3))));

            otherwise
                ThrowWrongDimensionException();
        end
        
        isGradient = (sum(abs(Jskew(:))) == 0); % True iff all elements are zero
        return;
                 
    case 'function_handle'
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
                    ThrowWrongDimensionException();
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
            otherwise
                
        end
        
end

end

function ThrowWrongDimensionException()
%THROWWRONGDIMENSIONEXCEPTION Throws an exception when systems with more
%than 3 dimensions are introduced

msgId = 'IsGradientField:WrongDimension';
errMsg = 'Only systems up to 3D are supported';
error(msgId, errMsg);

end