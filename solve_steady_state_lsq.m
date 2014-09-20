function xopt = solve_steady_state_lsq(alpha_ref)

%This function solves for steady states by solving a least squares problem.

options = optimoptions('lsqnonlin', 'TolFun', 1e-25,'TolX', 1e-25);
xopt = lsqnonlin(@(x)g(x, alpha_ref), [0, 0, -0.5*pi/180, 1.6, 1.6, 1.6], [], [], options);

end


function ret = g(dv, alpha_ref)

delta_motor = dv(1);
delta_arm = dv(2);
alpha = alpha_ref;
beta = dv(3);
ddelta_motor = dv(4);
ddelta_arm = dv(5);
ddelta_sp = dv(6);
u = 0;

x = [delta_motor; delta_arm; alpha; beta; ddelta_arm; ddelta_arm; 0; 0; ddelta_sp; 0];
xdot = carousel_lagrange(x,u) - [ddelta_arm; ddelta_arm; 0; 0; 0; 0; 0; 0; 0; 0];
ret = xdot;

end