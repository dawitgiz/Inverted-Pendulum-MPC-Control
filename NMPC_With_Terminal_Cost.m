
%% NONLINEAR MODEL PREDICTIVE CONTROL (NMPC)
% Cart-Pendulum System with LQR-Based Terminal Cost
%
% Features:
% 1. Correct nonlinear model around upright equilibrium
% 2. Linearization for LQR terminal cost
% 3. Discrete-time LQR terminal penalty
% 4. RK4 nonlinear prediction
% 5. fmincon SQP optimization
% 6. Warm-start optimization
% 7. Input constraints
% 8. ode45 nonlinear plant simulation
% 9. State plots and animation

clc;
clear;
close all;

%% =========================================================
% 1. SYSTEM PARAMETERS
% ==========================================================

M  = 1.0;       % Cart mass (kg)
m  = 0.2;       % Pendulum mass (kg)
lc = 0.3;       % Distance to pendulum center of mass (m)
Ic = 0.006;     % Pendulum moment of inertia (kg.m^2)

I = Ic + m*lc^2;

bc = 0.1;       % Cart damping coefficient
bp = 0.05;      % Pendulum damping coefficient

g = 9.81;       % Gravitational acceleration (m/s^2)

%% =========================================================
% 2. MPC PARAMETERS
% ==========================================================

Ts = 0.05;       % Sampling time (s)
N = 15;          % Prediction horizon
Tsim = 10;       % Simulation duration (s)

Nsim = round(Tsim/Ts);

Fmax = 15;       % Maximum control force (N)

Q = diag([500, 1, 400, 20]);
R = 0.2;

nx = 4;

x_ref = [0; 0; 0; 0];

%% =========================================================
% 3. LINEARIZATION AROUND UPRIGHT EQUILIBRIUM
% ==========================================================

D0 = (M+m)*I - (m*lc)^2;

A = [0, 1, 0, 0;
     0, -(bc*I)/D0, -(m^2*g*lc^2)/D0, (bp*m*lc)/D0;
     0, 0, 0, 1;
     0, (bc*m*lc)/D0, (g*(M+m)*m*lc)/D0, -(bp*(M+m))/D0];

B = [0;
     I/D0;
     0;
     -m*lc/D0];

%% =========================================================
% 4. DISCRETIZATION
% ==========================================================

sysc = ss(A,B,eye(nx),zeros(nx,1));

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
    error('Linearized system is not fully controllable.');
end

%% =========================================================
% 6. COMPUTE LQR TERMINAL COST MATRIX
% ==========================================================

try

    [K_lqr,P,~] = dlqr(Ad,Bd,Q,R);

catch ME

    warning('dlqr failed: %s',ME.message);

    try
        [P,~,~] = dare(Ad,Bd,Q,R);
    catch ME2
        error('Unable to calculate terminal matrix P: %s', ...
              ME2.message);
    end

end

% Ensure symmetry
P = (P+P')/2;

fprintf('\nTerminal cost matrix P:\n');
disp(P);

%% =========================================================
% 7. INITIAL CONDITIONS
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
% 8. OPTIMIZATION SETTINGS
% ==========================================================

fmin_opts = optimoptions('fmincon',...
    'Display','off',...
    'Algorithm','sqp',...
    'MaxIterations',200,...
    'OptimalityTolerance',1e-6,...
    'ConstraintTolerance',1e-6);

% Input constraints
lb = -Fmax*ones(N,1);
ub =  Fmax*ones(N,1);

% Initial optimization guess
U_guess = zeros(N,1);

%% =========================================================
% 9. DATA STORAGE
% ==========================================================

X = zeros(nx,Nsim);
U = zeros(1,Nsim);

exitflag_log = zeros(1,Nsim);
solve_time = zeros(1,Nsim);

%% =========================================================
% 10. MAIN NMPC CONTROL LOOP
% ==========================================================

fprintf('\nStarting NMPC simulation...\n');

for k = 1:Nsim

    % Define NMPC objective function
    costfun = @(U) nmpc_cost_with_terminal( ...
        U,x,N,Ts,Q,R,P,...
        M,m,lc,I,bc,bp,g);

    % Solve nonlinear optimization problem
    tic;

    [Uopt,~,exitflag] = fmincon( ...
        costfun,U_guess,...
        [],[],[],[],lb,ub,[],fmin_opts);

    solve_time(k) = toc;
    exitflag_log(k) = exitflag;

    % Select control action
    if exitflag > 0 && all(isfinite(Uopt))

        u = Uopt(1);

        % Shift solution for next optimization
        U_guess = [Uopt(2:end); Uopt(end)];

    else

        warning('Optimization failed at step %d. Exitflag = %d.',...
                 k,exitflag);

        % Use previous feasible input sequence
        u = U_guess(1);

    end

    % Enforce input saturation
    u = max(min(u,Fmax),-Fmax);

    % Apply control to actual nonlinear plant
    [~,Xstep] = ode45( ...
        @(t,xx) full_nonlinear( ...
        xx,u,M,m,lc,I,bc,bp,g),...
        [0 Ts],x);

    x = Xstep(end,:)';

    % Store results
    X(:,k) = x;
    U(k) = u;

end

fprintf('Simulation completed.\n');

%% =========================================================
% 11. SIMULATION TIME
% ==========================================================

t = (1:Nsim)*Ts;

%% =========================================================
% 12. PERFORMANCE INFORMATION
% ==========================================================

fprintf('\n========== NMPC PERFORMANCE ==========\n');

fprintf('Mean optimization time: %.6f seconds\n',...
        mean(solve_time));

fprintf('Maximum control force: %.4f N\n',...
        max(abs(U)));

fprintf('Angular IAE: %.6f rad.s\n',...
        sum(abs(X(3,:))*Ts));

fprintf('Position IAE: %.6f m.s\n',...
        sum(abs(X(1,:))*Ts));

fprintf('Control energy: %.6f N^2.s\n',...
        sum(U.^2)*Ts);

fprintf('Successful optimization steps: %d / %d\n',...
        sum(exitflag_log > 0),Nsim);

%% =========================================================
% 13. STATE RESPONSE PLOTS
% ==========================================================

figure('Name','NMPC State Responses',...
       'NumberTitle','off',...
       'Position',[200 100 850 850]);

subplot(5,1,1);

plot(t,X(1,:),'LineWidth',1.5);
ylabel('x (m)');
title('NMPC with Terminal Cost');
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
% 14. PHASE PORTRAIT
% ==========================================================

figure('Name','Pendulum Phase Portrait',...
       'NumberTitle','off');

plot(X(3,:)*180/pi,...
     X(4,:)*180/pi,...
     'LineWidth',1.5);

xlabel('Pendulum angle (deg)');
ylabel('Angular velocity (deg/s)');

title('Pendulum Phase Portrait');

grid on;

%% =========================================================
% 15. PENDULUM ANIMATION
% ==========================================================

figure('Name','NMPC Animation',...
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

    plot([-track_width track_width],[0 0],...
         'k','LineWidth',1.2);

    fill(cart_x,cart_y,[0.2 0.6 1],...
         'LineWidth',1.2);

    plot([px bob_x],[py bob_y],...
         'r','LineWidth',2.5);

    plot(bob_x,bob_y,'ko',...
         'MarkerFaceColor','k',...
         'MarkerSize',8);

    title(sprintf('NMPC - Time = %.2f s',t(k)));

    xlabel('Cart position (m)');
    ylabel('Height (m)');

    drawnow;

end

%% =========================================================
% 16. SAVE RESULTS
% ==========================================================

save('NMPC_Terminal_Cost_Results.mat',...
     'X','U','t','P','Q','R','Ts','N',...
     'solve_time','exitflag_log');

fprintf('\nResults saved successfully.\n');

%% =========================================================
% 17. NONLINEAR CART-PENDULUM MODEL
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
% 18. RK4 INTEGRATOR
% ==========================================================

function x_next = rk4_step(x,u,Ts,M,m,lc,I,bc,bp,g)

    k1 = full_nonlinear(x,u,M,m,lc,I,bc,bp,g);

    k2 = full_nonlinear( ...
        x+0.5*Ts*k1,u,M,m,lc,I,bc,bp,g);

    k3 = full_nonlinear( ...
        x+0.5*Ts*k2,u,M,m,lc,I,bc,bp,g);

    k4 = full_nonlinear( ...
        x+Ts*k3,u,M,m,lc,I,bc,bp,g);

    x_next = x+(Ts/6)*(k1+2*k2+2*k3+k4);

end

%% =========================================================
% 19. NMPC COST FUNCTION WITH TERMINAL PENALTY
% ==========================================================

function J = nmpc_cost_with_terminal( ...
    U,x0,N,Ts,Q,R,P,M,m,lc,I,bc,bp,g)

    x = x0;
    J = 0;

    x_ref = [0;0;0;0];

    for i = 1:N

        u = U(i);

        % Nonlinear prediction using RK4
        x = rk4_step(x,u,Ts,M,m,lc,I,bc,bp,g);

        err = x-x_ref;

        % Stage cost
        J = J+err'*Q*err+u'*R*u;

    end

    % Terminal penalty
    errN = x-x_ref;

    J = J+errN'*P*errN;

end
