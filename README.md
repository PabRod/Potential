# Potential
A couple of functions for computing path integrals and scalar potentials in an arbitrary number of dimensions

## Background
### Path integrals
A path integral, also known as line or curve integral, is an integral where the integrand is evaluated along a curve.

![GeneralPathIntegral](https://github.com/PabRod/Potential/blob/master/figs/path_general.png "General path integral")

The usual way of evaluating it requires the specification of a parameterization of the curve:

![GeneralParametric](https://github.com/PabRod/Potential/blob/master/figs/parameter_curve.png "General parametric curve")

So our path integral turns into the classical integral:

![PathIntegralEvaluation](https://github.com/PabRod/Potential/blob/master/figs/path_parametric.png "Path integral evaluation")

## Examples of usage
### Path integral along a parametric curve
```[Matlab]
# Underlying field
field = @(x) [-x(2), -x(1)];

# Parametric curve specification
curve = @(t) [t, t];
dcurve = @(t) [1, 1];
tmin = -2;
tmax = 1;

% Path integral
s = PathIntegral(field, curve, dcurve, tmin, tmax);
```
### Path integral over a conservative field
In the case of conservative fields the integral only depends in the initial and final points of the integration path. Thus, we can specify only those points regardless of the integration curve
```[Matlab]
# Underlying field
field = @(x) [-2.*x(1).*x(2), -x(1).^2];

# Initial and final points
x0 = [1 2];
x = [3 4];

# Path integral
s = PathIntegral(field, x0, x);
```

### Potential
Given a conservative field, we can compute the associated potential
```[Matlab]
# Underlying field
field = @(x) [ -x(2), -x(1)];

# Grid characteristics
x = -1:0.025:1;
y = -1:0.05:1.2;
[xm, ym] = meshgrid(x,y);

# Fixing of the potential at origin
x0 = zeros(1, 2);
V0 = Vexpected(x0(1), x0(2));

# Compute potential
V = Potential(field, V0, x0, xm, ym);
```

By Pablo Rodríguez-Sánchez, March 2017.
