%This script derives xdot2 = [dddelta_motor, dddelta_arm, ddalpha, ddbeta]'
%for the carousel model using the Lagrange formalism.

%ground frame: NED centered at the carousel axis at arm level

syms M_motor %torque of motor [Nm]
syms p_ball real %position of the ball in ground frame, [m]
syms r_arm real %length of arm [m]
syms l_tether real %length of tether [m]
syms m_ball real %mass of the ball [kg]
syms I_arm real %moment of inertia of carousel [kg m^2]
syms I_motor real %moment of inertia of motor [kg m^2]
syms I_tether real %moment of intertia of the tether [kg m^2]
syms p_sensor real %position of line angle sensor in ground frame [m]
syms delta_arm ddelta_arm dddelta_arm real %angle of arm rotating around carousel axis with 0 aligned with east, positive for counter-clockwise rotation [rad]
syms delta_motor ddelta_motor dddelta_motor real %angle of motor rotating around carousel axis with 0 aligned with east, positive for counter-clockwise rotation [rad]
syms alpha dalpha ddalpha real %elevation angle of rope, positive above horizontal [rad]
syms beta dbeta ddbeta real %azimuth angle of rope, positive if ball is ahead of the arm [rad]
syms k_beltspring real %spring constant of the belt [Nm/rad]
syms c_beltdampening real %dampening coefficient of the belt [Nm/rad^2]
syms my_shaft real %friction constant of motor shaft [Nm/rad]
syms my_beta_LA real %friction constant of LA-sensor in beta direction [Nm/rad]
syms my_alpha_LA real %friction constant of LA-sensor in alpha direction [Nm/rad]
syms g real %gravitational constant
syms roh_air real %density of the air
syms A_ball real %area of the ball
syms c_w %air friction constant
syms k_p %proportional constant of the motor control [-]
syms ddelta_motor_sp %setpoint for the motor controller [rad/s]

%defining generalized coordinates
q = [delta_motor; delta_arm; alpha; beta];
dq = [ddelta_motor; ddelta_arm; dalpha; dbeta];
ddq = [dddelta_motor; dddelta_arm; ddalpha; ddbeta];

%defining geometric subexpressions
p_sensor = r_arm*[sin(delta_arm);cos(delta_arm);0];
p_ball = p_sensor + l_tether*[cos(alpha)*sin(delta_arm + beta); cos(alpha)*cos(delta_arm + beta); -sin(alpha)];

%velocity of the ball
v_ball = jacobian(p_ball, q)*dq;

%kinetic energy
KE = 0;
KE = KE + 1/2*m_ball*(v_ball.'*v_ball);
KE = KE + 1/2*I_arm*ddelta_arm^2;
KE = KE + 1/2*I_motor*ddelta_motor^2;
KE = KE + 1/2*I_tether*(dalpha^2+dbeta^2);

%potential energy
PE = m_ball*g*(l_tether*sin(alpha));
PE = PE + 1/2*k_beltspring*(delta_motor-delta_arm)^2;

%lagrangian
L = KE - PE;

%define motor control
M_motor = -k_p*(ddelta_motor - ddelta_motor_sp);

%generalized forces
J = jacobian(p_ball,q);
f_airfriction = -1/2*roh_air*A_ball*c_w*v_ball*norm(v_ball);
M_beltfriction = -c_beltdampening*(ddelta_motor - ddelta_arm);
M_shaftfriction = -my_shaft;
gen_airfriction = J.'*f_airfriction;
M_LA_friction_alpha = -my_alpha_LA*dalpha;
M_LA_friction_beta = -my_beta_LA*dbeta;
gen_forces = [14.86*(3/2)*M_motor+M_beltfriction; -M_beltfriction + M_shaftfriction; M_LA_friction_alpha; M_LA_friction_beta] + gen_airfriction;

%implicit ODE
Lq = jacobian(L, q).';
Ldq = jacobian(L, dq);
Ldqt = jacobian(Ldq, q)*dq + jacobian(Ldq, dq)*ddq;
impl_ode = Ldqt - Lq - gen_forces;

%explicit ODE
sol = solve(impl_ode(1), impl_ode(2), impl_ode(3), impl_ode(4), dddelta_motor, dddelta_arm, ddalpha, ddbeta);
xdot2 = simplify(expand([sol.dddelta_motor; sol.dddelta_arm; sol.ddalpha; sol.ddbeta]));