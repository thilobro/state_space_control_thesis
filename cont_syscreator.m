function [A, B, C] = cont_syscreator(x0, u0)

%This function returns the matrices A, B and C of a continous linearization
%of carousel_lagrange at x0, u0.

A = fingrad(@(x)carousel_lagrange(x,u0), x0, 1e-6);
B = fingrad(@(u)carousel_lagrange(x0,u), u0, 1e-6);
C = fingrad(@sensor, x0, 1e-6);
Contr = ctrb(A ,B);

if rank(Contr) < rank(A)
    disp('Syscreator: System not controllable.')
end

end