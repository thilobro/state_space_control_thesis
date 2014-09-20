function [T,X,U] = carousel_dynamics_SE(carouselspeed_sp, timestamp)

%This function takes the experimental data for the ddelta_motor setpoint as
%controls and simulates the system.

%set time step and number of steps
ts = 0.1;
n = 1000;

%set initial conditions
ddelta0 = 1.44;
xss0 = [0, 0, -57*pi/180, 0, ddelta0, ddelta0, 0, 0]';

%initialize matrices
nx = length(xss0);
X = zeros(nx, n + 1);
X(:, 1) = xss0;
T = 0:ts:(n*ts);

%resample data
carouselspeedsetpointresamp = interp1(timestamp,carouselspeed_sp,T,'spline');

%set U
U = carouselspeedsetpointresamp(200:end - 1);

for k = 1:(n - 200)
    x = X(:, k);
    u = U(k);
    x = integrator(x, u, ts);
    X(:, k + 1) = x;
end

end