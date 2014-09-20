function [A, B, C] = disc_syscreator(x0, u0, ts)

%This function returns the matrices A, B and C of a discrete linearization
%of carousel_lagrange at x0, u0.

A = fingrad(@(x)integrator(x,u0,ts), x0, 1e-6);
B = fingrad(@(u)integrator(x0,u,ts), u0, 1e-6);
C = fingrad(@sensor, x0, 1e-6);
Contr = ctrb(A, B);

if rank(Contr) < rank(A)
    disp('Syscreator: System not controllable.')
end

end