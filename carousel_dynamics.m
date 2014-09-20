%This script simulates and plots the system behavior for a repeated step 
%from alpha0 to alpha1. It is possible to choose between a Discrete Linear
%Quadratic Regulator and a PID controller. In order to change the control
%for open loop simulations or simulate without the pseudo-force s, one has
%to change the matrix dimensions of the steady states and LQG matrices
%accordingly.

%constants
k_p = 400; %internal PID of model

%specify alpha reference values
alpha0 = -50*pi/180;
alpha1 = -55*pi/180;
alphass = (alpha0 + alpha1)/2;

%compute steady states
xopt0 = solve_steady_state_lsq(alpha0);
xopt1 = solve_steady_state_lsq(alpha1);
xoptss = solve_steady_state_lsq(alphass);
xss = [xoptss(1), xoptss(2), alphass, xoptss(3), xoptss(5), xoptss(5), 0, 0,xoptss(6), 0].';
xss0 = [xopt0(1), xopt0(2), alpha0, xopt0(3), xopt0(5), xopt0(5), 0, 0,xopt0(6), 0].';
xss1 = [xopt1(1), xopt1(2), alpha1, xopt1(3), xopt1(5), xopt1(5), 0, 0,xopt1(6), 0].';

%check if steady state were computed correctly
xdot_zero0 = carousel_lagrange(xss0, 0)
xdot_zero1 = carousel_lagrange(xss1, 0)

%set size of time step ts and number of time steps
n=700;
ts=0.1;

%get system matrices A, B and C by linerization
[A,B,C] = disc_syscreator(xss, 0, ts);

%check controllability and observability
rank_ctrb = rank(ctrb(A, B))
rank_obsv = rank(obsv(A, C)) 

%set PID gains
%Kp_pid = 0.3;
%Ki_pid = 0.001;
%Kd_pid = 0.08;

%dimension of state vector nx and output vector ny
nx = length(xss0);
ny = size(C, 1);

%set matrices Q and R for DLQR
penalty_delta_motor = 0;
penalty_delta_arm = 0;
penalty_alpha = 8;
penalty_beta = 0;
penalty_ddelta_motor = 0;
penalty_ddelta_arm = 0;
penalty_dalpha = 5;
penalty_dbeta = 0;
penalty_ddelta_motor_sp = 0;
penalty_s = 0;
Q = diag([penalty_delta_motor, penalty_delta_arm, penalty_alpha, penalty_beta, penalty_ddelta_motor, penalty_ddelta_arm, penalty_dalpha, penalty_dbeta, penalty_ddelta_motor_sp, penalty_s].^2);
R = (2).^2;

%set matrices G, QE and RE for DLQE (Kalman filter)
w_n = 0.001;
G = eye(nx);
variance_delta_motor = 0;
variance_delta_arm = 0;
variance_alpha = 0;
variance_beta = 0;
variance_ddelta_motor = w_n;
variance_ddelta_arm = w_n;
variance_dalpha = w_n;
variance_dbeta = w_n;
variance_ddelta_motor_sp = 0;
variance_s = w_n;
QE = diag([variance_delta_motor, variance_delta_arm, variance_alpha, variance_beta, variance_ddelta_motor, variance_ddelta_arm, variance_dalpha, variance_dbeta, variance_ddelta_motor_sp, variance_s].^2);
RE = diag([0.01,0.01].^2);
QE_sim = diag([variance_delta_motor, variance_delta_arm, 0, variance_beta, variance_ddelta_motor, variance_ddelta_arm, variance_dalpha, variance_dbeta, variance_ddelta_motor_sp, 0].^2);
RE_sim = diag([0.01,0.01].^2);

%get matrix K for controls and M for Kalman filter
[K, ~, ~] = dlqr(A, B, Q, R);
[M, ~, ~, ~] = dlqe(A, G, C, QE, RE);

%Create matrices X, Xref, Xest, U and T to log simulation results
X = zeros(nx, n + 1);
Xest = X;
U = zeros(1, n);
Ysens = zeros(ny, n);
T = 0:ts:(n*ts);

%Create matrices for integral control
ERROR = zeros(ny,n);

%Initialize logging matrices
Xest(:,1) = xss0;
X(:,1) = xss0;
Xref = zeros(nx,n);

for k = 1:n
    
    %get simulated sensor
    Ysens(:, k) = sensor(X(:, k)) + mvnrnd(zeros([2, 1]), RE_sim)';
    
    %run estimator
    yest = sensor(Xest(:, k));
    Xest(:, k) = Xest(:, k) + M*(Ysens(:, k) - yest);
    
    %reference value changing periodically from xss0 to xss1
    if mod(T(k),40) < 20
        xref = xss0;
    else
        xref = xss1;
    end
    
    %update reference value for delta_motor and delta_arm
    if k > 2
        xref(1) = Xest(1, k - 1) + Xest(5, k)*ts;
        xref(2) = Xest(2 ,k - 1) + Xest(6, k)*ts;
        
    end
    
    %log reference value
    Xref(:, k) = xref;
    
    %apply full state feedback controller
    U(k) = -K*(Xest(:, k)-Xref(:, k));
    
    %apply PID controller
%     ERROR(:, k + 1) = ERROR(:, k) + [X(3, k) - xref(3); X(4, k) - xref(4)];
%     bar_ak = (X(3, k) - Xref(3, k));
%     if k > 1
%         bar_akminus1 = X(3, k - 1) - Xref(3, k - 1);
%         U(k) = -(Kp_pid*(bar_ak) + Ki_pid*ts*ERROR(1, k) + Kd_pid*(bar_ak - bar_akminus1)/ts);
%     else
%         U(k) = -(Kp_pid*bar_ak + Ki_pid*ts*ERROR(1, k));
%     end
    
    %simulate system
    X(:, k + 1) = integrator(X(:,k),U(k),ts) + mvnrnd(zeros(nx, 1), QE_sim)';    
    
    %run estimator
    Xest(:, k + 1) = integrator(Xest(:, k), U(k), ts); %non-linear
    %Xest(:, k + 1) = xss + A*(Xest(:, k) - xss) + B*U(k); %linear
         
end

%plots
figure(1);
clf;

%plot elevation
ax(1) = subplot(3, 2, 1);
hold on;
plot(T, 180/pi*X(3, :), 'b');
plot(T, 180/pi*Xest(3, :), 'g');
axis([0 55 -60 -45])
xlabel('t [s]')
ylabel('Elevation [deg]')
plot(T(1:end-1), 180/pi*Xref(3, :), 'r')
legend('Simulation Data', 'Estimator', 'Reference Value')
grid on;

%plot arm speed
ax(2) = subplot(3, 2, 2);
hold on;
plot(T, X(6, :));
plot(T(1:end - 1), Xref(6, :), 'r')
plot(T, Xest(6, :), 'g');
grid on;
axis([0 55 1.45 1.65])
xlabel('t [s]')
ylabel('Arm Speed [rad/s]')
legend('Simulation Data', 'Estimator', 'Reference Value')

%plot motor speed
ax(3) = subplot(3, 2, 3);
hold on;
grid on;
axis([0 55 1.45 1.65])
plot(T, X(5, :));
plot(T(1:end - 1), Xref(5, :), 'r')
plot(T, Xest(5, :), 'g');
xlabel('t [s]')
ylabel('Motor Speed [rad/s]')
legend('Simulation Data', 'Estimator', 'Reference Value')

%plot u
ax(4) = subplot(3, 2, 4);
hold on;
grid on;
plot(T(1:end - 1), U);
xlabel('t [s]')
ylabel('Change of Motor Speed SP [rad/s^2]')

%plot azimuth
ax(5) = subplot(3, 2, 5);
hold on;
grid on;
axis([0 55 -1.5 1.5])
plot(T(1:end - 1), Xref(4, :)*180/pi, 'r')
plot(T, X(4, :)*180/pi, 'b')
plot(T, 180/pi*Xest(4, :), 'g')
xlabel('t [s]')
ylabel('Azimuth [deg]')
legend('Simulation Data', 'Estimator', 'Reference Value')

%plot motor torque
ax(6) = subplot(3, 2, 6);
hold on;
grid on;
plot(T, -k_p*(X(5, :) - X(9, :)), 'b')
axis([0 55 -20 20])
xlabel('t [s]')
ylabel('Motor Torque [Nm]')

linkaxes(ax, 'x')