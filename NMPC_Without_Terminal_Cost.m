
%% NONLINEAR MODEL PREDICTIVE CONTROL (NMPC)
% Cart-Pendulum System WITHOUT Terminal Cost
%
% Refined version of the original implementation
%
% - Correct upright nonlinear model
% - RK4 prediction
% - fmincon SQP optimization
% - Warm-start optimization
% - Input saturation
% - ode45 plant simulation
% - Performance evaluation

clc;
clear;
close all;

%% =========================================================
% 1. PARAMETERS
% ==========================================================

M  = 1.0;
m  = 0.2;
lc = 0.3;
Ic = 0.006;

I = Ic + m*lc^2;

bc = 0.1;
bp = 0.05;

g = 9.81;

%% =========================================================
% 2. MPC PARAMETERS
% ==========================================================

Ts = 0.05;
N = 15;

Fmax = 15;

Q = diag([500,1,400,20]);
R = 0.2;

nx = 4;

Tsim = 10;
Nsim = round(Tsim/Ts);

x_ref = [0;0;0;0];

%% =========================================================
% 3. INITIAL CONDITIONS
% ==========================================================

x0 = [0;
      0;
      20*pi/180;
      0];

x = x0;

%% =========================================================
% 4. OPTIMIZATION SETTINGS
% ==========================================================

fmin_opts = optimoptions('fmincon',...
    'Display','off',...
    'Algorithm','sqp',...
    'MaxIterations',200,...
    'OptimalityTolerance',1e-6,...
    'ConstraintTolerance',1e-6);

lb = -Fmax*ones(N,1);
ub =  Fmax*ones(N,1);

% Warm-start initial guess
U_guess = zeros(N,1);

%% =========================================================
% 5. DATA STORAGE
% ==========================================================

X = zeros(nx,Nsim);
U = zeros(1,Nsim);

solve_time = zeros(1,Nsim);
exitflag_log = zeros(1,Nsim);

%% =========================================================
% 6. MAIN NMPC SIMULATION LOOP
% ==========================================================

fprintf('\nStarting conventional NMPC simulation...\n');

for k = 1:Nsim

    % NMPC objective function
    costfun = @(U) nmpc_cost(...
        U,x,N,Ts,Q,R,M,m,lc,I,bc,bp,g);

    % Solve optimization problem
    tic;

    [Uopt,~,exitflag] = fmincon(...
        costfun,U_guess,...
        [],[],[],[],lb,ub,[],fmin_opts);

    solve_time(k) = toc;
    exitflag_log(k) = exitflag;

    % Check optimization result
    if exitflag > 0 && ~isempty(Uopt) && all(isfinite(Uopt))

        u = Uopt(1);

        % Shift solution for next iteration
        U_guess = [Uopt(2:end);Uopt(end)];

    else

        warning('Optimization failed at step %d. Exitflag = %d.',...
                 k,exitflag);

        % Previous feasible sequence as fallback
        u = U_guess(1);

    end

    % Input saturation
    u = max(min(u,Fmax),-Fmax);

    % Simulate actual nonlinear plant
    [~,Xstep] = ode45(...
        @(t,xx) full_nonlinear(...
        xx,u,M,m,lc,I,bc,bp,g),...
        [0 Ts],x);

    x = Xstep(end,:)';

    % Store results
    X(:,k) = x;
    U(k) = u;

end

fprintf('Simulation completed.\n');

%% =========================================================
% 7. TIME VECTOR
% ==========================================================

t = (1:Nsim)*Ts;

%% =========================================================
% 8. PERFORMANCE EVALUATION
% ==========================================================

IAE_angle = sum(abs(X(3,:)))*Ts;

IAE_position = sum(abs(X(1,:)))*Ts;

Control_energy = sum(U.^2)*Ts;

Maximum_force = max(abs(U));

Mean_solve_time = mean(solve_time);

settling_time_value = calculate_settling_time(...
    t,X,1*pi/180,0.02);

fprintf('\n========== NMPC PERFORMANCE ==========\n');

fprintf('Angular IAE: %.6f rad.s\n',IAE_angle);

fprintf('Position IAE: %.6f m.s\n',IAE_position);

fprintf('Control energy: %.6f N^2.s\n',Control_energy);

fprintf('Maximum force: %.4f N\n',Maximum_force);

fprintf('Mean solve time: %.6f seconds\n',Mean_solve_time);

fprintf('Settling time: %.4f seconds\n',settling_time_value);

fprintf('Successful optimization steps: %d / %d\n',...
        sum(exitflag_log > 0),Nsim);

%% =========================================================
% 9. STATE RESPONSE PLOTS
% ==========================================================

figure('Name','NMPC State Responses',...
       'NumberTitle','off',...
       'Position',[200 100 850 850]);

subplot(5,1,1);

plot(t,X(1,:),'LineWidth',1.5);
ylabel('x (m)');
title('NMPC Without Terminal Cost');
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
% 10. PHASE PORTRAIT
% ==========================================================

figure('Name','NMPC Phase Portrait',...
       'NumberTitle','off');

plot(X(3,:)*180/pi,...
     X(4,:)*180/pi,...
     'LineWidth',1.5);

xlabel('Pendulum angle (deg)');
ylabel('Angular velocity (deg/s)');

title('NMPC Phase Portrait');

grid on;

%% =========================================================
% 11. ANIMATION
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
         'k','LineWidth',1.5);

    fill(cart_x,cart_y,[0.2 0.6 1],...
         'LineWidth',1.2);

    plot([px bob_x],[py bob_y],...
         'r','LineWidth',2.5);

    plot(bob_x,bob_y,'ko',...
         'MarkerFaceColor','k',...
         'MarkerSize',10);

    title(sprintf('NMPC - Time = %.2f s',t(k)));

    xlabel('Cart position (m)');
    ylabel('Height (m)');

    drawnow;

end

%% =========================================================
% 12. SAVE RESULTS
% ==========================================================

save('NMPC_Without_Terminal_Cost_Results.mat',...
     'X','U','t','Q','R','Ts','N',...
     'solve_time','exitflag_log',...
     'IAE_angle','IAE_position',...
     'Control_energy','Maximum_force',...
     'Mean_solve_time','settling_time_value');

fprintf('\nResults saved successfully.\n');

%% =========================================================
% 13. NONLINEAR CART-PENDULUM MODEL
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
% 14. RK4 INTEGRATOR
% ==========================================================

function x_next = rk4_step(x,u,Ts,M,m,lc,I,bc,bp,g)

    k1 = full_nonlinear(x,u,M,m,lc,I,bc,bp,g);

    k2 = full_nonlinear(...
        x+0.5*Ts*k1,u,M,m,lc,I,bc,bp,g);

    k3 = full_nonlinear(...
        x+0.5*Ts*k2,u,M,m,lc,I,bc,bp,g);

    k4 = full_nonlinear(...
        x+Ts*k3,u,M,m,lc,I,bc,bp,g);

    x_next = x+(Ts/6)*(k1+2*k2+2*k3+k4);

end

%% =========================================================
% 15. NMPC COST FUNCTION WITHOUT TERMINAL COST
% ==========================================================

function J = nmpc_cost(U,x0,N,Ts,Q,R,M,m,lc,I,bc,bp,g)

    x = x0;
    J = 0;

    x_ref = [0;0;0;0];

    for i = 1:N

        u = U(i);

        % RK4 nonlinear prediction
        x = rk4_step(x,u,Ts,M,m,lc,I,bc,bp,g);

        err = x-x_ref;

        % Stage cost
        J = J+err'*Q*err+u'*R*u;

    end

end

%% =========================================================
% 16. SETTLING TIME FUNCTION
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
