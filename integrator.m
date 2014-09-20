function xnext = integrator(x0, u0, dt)

%This function integrates carousel_lagrange numerically with initial value
%x0, u0 over time dt.

options = odeset('RelTol', 5e-14);
[~, Y] = ode15s(@(t, x)carousel_lagrange(x, u0), [0, dt], x0, options);
xnext = Y(end, :)';

end