
%% LINEAR MODEL PREDICTIVE CONTROL (LINEAR MPC)
% Cart-Pendulum System
%
% Features:
% 1. Correct linearization around upright equilibrium
% 2. Discrete-time state-space model
% 3. Quadratic programming (quadprog)
% 4. Input constraints
% 5. Nonlinear plant simulation using ode45
% 6. Performance evaluation
% 7. State plots and animation
%
% Note:
% This controller uses a linear prediction model and does not
% include a terminal cost. The actual plant is nonlinear.

clc;
clear;
close all;

%% =========================================================
% 1. SYSTEM PARAMETERS
% ==========================================================

M  = 1.0;       % Cart mass (kg)
m  = 0.2;       % Pendulum mass (kg)
lc = 0.3;       % Distance to center of mass (m)
Ic = 0.006;     % Pendulum moment of inertia (kg.m^2)

I = Ic + m*lc^2;

bc = 0.1;       % Cart damping
bp = 0.05;      % Pendulum damping

g = 9.81;       % Gravity (m/s^2)

%% =========================================================
% 2. MPC PARAMETERS
% ==========================================================

Ts = 0.05;      % Sampling time (s)
N = 15;         % Prediction horizon

Tsim = 10;      % Simulation duration (s)
Nsim = round(Tsim/Ts);

Fmax = 15;      % Maximum control force (N)

Q = diag([500, 1, 400, 20]);
R = 0.2;

nx = 4;
nu = 1;

x_ref = [0; 0; 0; 0];

%% =========================================================
% 3. LINEARIZATION AROUND UPRIGHT EQUILIBRIUM
% ==========================================================

D0 = (M+m)*I - (m*lc)^2;

A = [0, 1, 0, 0;
     0, -(bc*I)/D0, -(m^2*g*lc^2)/D0,  (bp*m*lc)/D0;
     0, 0, 0, 1;
     0, (bc*m*lc)/D0, (g*(M+m)*m*lc)/D0, -(bp*(M+m))/D0];

B = [0;
     I/D0;
     0;
     -m*lc/D0];

%% =========================================================
% 4. DISCRETIZATION
% ==========================================================

sysc = ss(A,B,eye(nx),zeros(nx,nu));

sysd = c2d(sysc,Ts);

Ad = sysd.A;
Bd = sysd.B;

%% =========================================================
% 5. CONTROLLABILITY CHECK
% ==========================================================

Ctrb = ctrb(Ad,Bd);

rankCtrb = rank(Ctrb);

fprintf('Controllability rank = %d / %d\n',rankCtrb,nx);

if rankCtrb < nx
    error('The linearized system is not fully controllable.');
end

%% =========================================================
% 6. BUILD MPC PREDICTION MATRICES
% ==========================================================

Qbar = kron(eye(N),Q);
Rbar = kron(eye(N),R);

Abar = zeros(nx*N,nx);
Bbar = zeros(nx*N,nu*N);

for i = 1:N

    rows = (i-1)*nx+1:i*nx;

    Abar(rows,:) = Ad^i;

    for j = 1:i

        cols = (j-1)*nu+1:j*nu;

        Bbar(rows,cols) = Ad^(i-j)*Bd;

    end

end

%% =========================================================
% 7. QUADRATIC PROGRAMMING MATRICES
% ==========================================================

H = 2*(Bbar'*Qbar*Bbar + Rbar);

% Ensure numerical symmetry
H = (H+H')/2;

% Input constraints
lb = -Fmax*ones(nu*N,1);
ub =  Fmax*ones(nu*N,1);

% QP solver options
qp_options = optimoptions('quadprog',...
    'Display','off',...
    'Algorithm','interior-point-convex');

%% =========================================================
% 8. INITIAL CONDITIONS
% ==========================================================

% State vector:
% x(1) = Cart position (m)
% x(2) = Cart velocity (m/s)
% x(3) = Pendulum angle (rad)
% x(4) = Angular velocity (rad/s)

x0 = [0;
      0;
      20*pi/180;
      0];

x = x0;

%% =========================================================
% 9. DATA STORAGE
% ==========================================================

X = zeros(nx,Nsim);
U = zeros(1,Nsim);

solve_time = zeros(1,Nsim);
exitflag_log = zeros(1,Nsim);

%% =========================================================
% 10. NONLINEAR PLANT MODEL
% ==========================================================

nonlinear_dyn = @(x,u) full_nonlinear(...
    x,u,M,m,lc,I,bc,bp,g);

%% =========================================================
% 11. MAIN LINEAR MPC CONTROL LOOP
% ==========================================================

fprintf('\nStarting Linear MPC simulation...\n');

u_previous = 0;

for k = 1:Nsim

    % Current state deviation from reference
    x_error = x-x_ref;

    % Linear MPC objective function
    f = 2*Bbar'*Qbar*Abar*x_error;

    % Solve quadratic programming problem
    tic;

    [Uopt,~,exitflag] = quadprog(...
        H,f,[],[],[],[],lb,ub,[],qp_options);

    solve_time(k) = toc;
    exitflag_log(k) = exitflag;

    % Select first optimal control input
    if exitflag > 0 && ~isempty(Uopt) && all(isfinite(Uopt))

        u = Uopt(1);

    else

        warning('QP optimization failed at step %d. Exitflag = %d.',...
                 k,exitflag);

        % Use previous control input as fallback
        u = u_previous;

    end

    % Input saturation
    u = max(min(u,Fmax),-Fmax);

    u_previous = u;

    % Apply control to actual nonlinear plant
    [~,Xstep] = ode45(...
        @(t,xx) nonlinear_dyn(xx,u),...
        [0 Ts],x);

    x = Xstep(end,:)';

    % Store results
    X(:,k) = x;
    U(k) = u;

end

fprintf('Simulation completed.\n');

%% =========================================================
% 12. SIMULATION TIME
% ==========================================================

t = (1:Nsim)*Ts;

%% =========================================================
% 13. PERFORMANCE EVALUATION
% ==========================================================

IAE_angle = sum(abs(X(3,:)))*Ts;

IAE_position = sum(abs(X(1,:)))*Ts;

Control_energy = sum(U.^2)*Ts;

Maximum_force = max(abs(U));

Mean_solve_time = mean(solve_time);

% Settling conditions:
% Pendulum angle within +/- 1 degree
% Cart position within +/- 0.02 m

settling_time_value = calculate_settling_time(...
    t,X,1*pi/180,0.02);

fprintf('\n========== LINEAR MPC PERFORMANCE ==========\n');

fprintf('Angular IAE: %.6f rad.s\n',IAE_angle);

fprintf('Position IAE: %.6f m.s\n',IAE_position);

fprintf('Control energy: %.6f N^2.s\n',Control_energy);

fprintf('Maximum force: %.4f N\n',Maximum_force);

fprintf('Mean solve time: %.6f seconds\n',Mean_solve_time);

fprintf('Settling time: %.4f seconds\n',settling_time_value);

fprintf('Successful QP steps: %d / %d\n',...
        sum(exitflag_log > 0),Nsim);

%% =========================================================
% 14. STATE RESPONSE PLOTS
% ==========================================================

figure('Name','Linear MPC State Responses',...
       'NumberTitle','off',...
       'Position',[200 100 850 850]);

subplot(5,1,1);

plot(t,X(1,:),'LineWidth',1.5);
ylabel('x (m)');
title('Linear MPC - Cart-Pendulum System');
grid on;

subplot(5,1,2);

plot(t,X(2,:),'LineWidth',1.5);
ylabel('xdot (m/s)');
grid on;

subplot(5,1,3);

plot(t,X(3,:)*180/pi,'LineWidth',1.5);
ylabel('\theta (deg)');
grid on;

subplot(5,1,4);

plot(t,X(4,:)*180/pi,'LineWidth',1.5);
ylabel('theta dot (deg/s)');
grid on;

subplot(5,1,5);

stairs(t,U,'LineWidth',1.5);
ylabel('Force (N)');
xlabel('Time (s)');
grid on;

%% =========================================================
% 15. PHASE PORTRAIT
% ==========================================================

figure('Name','Linear MPC Phase Portrait',...
       'NumberTitle','off');

plot(X(3,:)*180/pi,...
     X(4,:)*180/pi,...
     'LineWidth',1.5);

xlabel('Pendulum angle (deg)');
ylabel('Angular velocity (deg/s)');

title('Pendulum Phase Portrait - Linear MPC');

grid on;

%% =========================================================
% 16. PENDULUM ANIMATION
% ==========================================================

figure('Name','Linear MPC Animation',...
       'NumberTitle','off');

track_width = 3;

cart_w = 0.3;
cart_h = 0.15;

L = 2*lc;

for k = 1:2:Nsim

    xc = X(1,k);
    yc = 0;

    cart_x = [xc-cart_w/2,...
              xc+cart_w/2,...
              xc+cart_w/2,...
              xc-cart_w/2];

    cart_y = [yc,...
              yc,...
              yc+cart_h,...
              yc+cart_h];

    px = xc;
    py = yc+cart_h;

    bob_x = px+L*sin(X(3,k));
    bob_y = py+L*cos(X(3,k));

    cla;

    hold on;
    grid on;
    axis equal;

    axis([-track_width track_width -0.5 1.3]);

    % Track
    plot([-track_width track_width],[0 0],...
         'k','LineWidth',1.5);

    % Cart
    fill(cart_x,cart_y,[0.2 0.6 1],...
         'LineWidth',1.2);

    % Pendulum
    plot([px bob_x],[py bob_y],...
         'r','LineWidth',2.5);

    plot(bob_x,bob_y,'ko',...
         'MarkerFaceColor','k',...
         'MarkerSize',10);

    title(sprintf('Linear MPC - Time = %.2f s',t(k)));

    xlabel('Cart position (m)');
    ylabel('Height (m)');

    drawnow;

end

%% =========================================================
% 17. SAVE RESULTS
% ==========================================================

save('Linear_MPC_Results.mat',...
     'X','U','t','Q','R','Ts','N',...
     'solve_time','exitflag_log',...
     'IAE_angle','IAE_position',...
     'Control_energy','Maximum_force',...
     'Mean_solve_time','settling_time_value');

fprintf('\nLinear MPC results saved successfully.\n');

%% =========================================================
% 18. NONLINEAR CART-PENDULUM MODEL
% ==========================================================

function dx = full_nonlinear(x,u,M,m,lc,I,bc,bp,g)

    theta = x(3);
    th_dot = x(4);
    xd = x(2);

    Mmat = [M+m, m*lc*cos(theta);
            m*lc*cos(theta), I];

    f = [u-bc*xd-m*lc*th_dot^2*sin(theta);
         m*g*lc*sin(theta)-bp*th_dot];

    acc = Mmat\f;

    dx = [xd;
          acc(1);
          th_dot;
          acc(2)];

end

%% =========================================================
% 19. SETTLING TIME CALCULATION
% ==========================================================

function ts = calculate_settling_time(t,X,angle_limit,...
                                      position_limit)

    angle_ok = abs(X(3,:)) <= angle_limit;

    position_ok = abs(X(1,:)) <= position_limit;

    settled = angle_ok & position_ok;

    ts = NaN;

    for k = 1:length(t)

        if all(settled(k:end))

            ts = t(k);
            return;

        end

    end

end
