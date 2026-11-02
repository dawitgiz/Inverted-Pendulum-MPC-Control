
clc;
clear;
close all;

%% =====================================================
% THREE-CONTROLLER COMPARISON
% 1. Linear MPC
% 2. Nonlinear MPC without terminal cost
% 3. Nonlinear MPC with LQR terminal cost
%
% All controllers use the SAME nonlinear plant.
% Upright equilibrium: theta = 0.
%% =====================================================

%% 1. SYSTEM PARAMETERS

M  = 1.0;
m  = 0.2;
lc = 0.3;
Ic = 0.006;
I  = Ic + m*lc^2;

bc = 0.1;
bp = 0.05;
g  = 9.81;

Ts = 0.05;
N = 15;
Tsim = 10;
Nsim = round(Tsim/Ts);

Q = diag([500,1,400,20]);
R = 0.2;

Fmax = 15;

nx = 4;

% Initial condition
x0 = [0; 0; 20*pi/180; 0];

%% 2. CORRECT LINEARIZATION ABOUT UPRIGHT

D0 = (M+m)*I - (m*lc)^2;

A = [0, 1, 0, 0;
     0, -(bc*I)/D0, -(m^2*g*lc^2)/D0, (bp*m*lc)/D0;
     0, 0, 0, 1;
     0, (bc*m*lc)/D0, ...
     (g*(M+m)*m*lc)/D0, -(bp*(M+m))/D0];

B = [0;
     I/D0;
     0;
     -m*lc/D0];

sysc = ss(A,B,eye(nx),zeros(nx,1));
sysd = c2d(sysc,Ts);

Ad = sysd.A;
Bd = sysd.B;

%% 3. LQR TERMINAL COST

[P_lqr,~,~] = dare(Ad,Bd,Q,R);
P_lqr = (P_lqr + P_lqr')/2;

fprintf('LQR terminal cost calculated successfully.\n');

%% 4. BUILD LINEAR MPC MATRICES

Qbar = kron(eye(N),Q);
Rbar = kron(eye(N),R);

Abar = zeros(nx*N,nx);
Bbar = zeros(nx*N,N);

for i = 1:N

    rows = (i-1)*nx+1:i*nx;

    Abar(rows,:) = Ad^i;

    for j = 1:i

        cols = j;

        Bbar(rows,cols) = Ad^(i-j)*Bd;

    end
end

H = 2*(Bbar'*Qbar*Bbar + Rbar);
H = (H+H')/2;

lb = -Fmax*ones(N,1);
ub =  Fmax*ones(N,1);

qp_options = optimoptions('quadprog',...
    'Display','off');

%% 5. INITIALIZE STORAGE

X_L = zeros(nx,Nsim);
X_N = zeros(nx,Nsim);
X_T = zeros(nx,Nsim);

U_L = zeros(1,Nsim);
U_N = zeros(1,Nsim);
U_T = zeros(1,Nsim);

time_L = zeros(1,Nsim);
time_N = zeros(1,Nsim);
time_T = zeros(1,Nsim);

exit_N = zeros(1,Nsim);
exit_T = zeros(1,Nsim);

xL = x0;
xN = x0;
xT = x0;

% NMPC warm-start sequences
Uguess_N = zeros(N,1);
Uguess_T = zeros(N,1);

%% 6. NONLINEAR OPTIMIZATION OPTIONS

fmin_opts = optimoptions('fmincon',...
    'Display','off',...
    'Algorithm','sqp',...
    'MaxIterations',200,...
    'OptimalityTolerance',1e-6);

%% =====================================================
% 7. MAIN SIMULATION
%% =====================================================

fprintf('\nStarting three-controller simulation...\n');

for k = 1:Nsim

    %% -----------------------------------------------
    % CONTROLLER 1: LINEAR MPC
    %% -----------------------------------------------

    tic;

    f = 2*Bbar'*Qbar*Abar*xL;

    Uopt_L = quadprog(H,f,[],[],[],[],...
        lb,ub,[],qp_options);

    time_L(k) = toc;

    if isempty(Uopt_L)
        error('Linear MPC optimization failed at step %d.',k);
    end

    uL = Uopt_L(1);

    U_L(k) = uL;

    [~,Xstep] = ode45(@(t,x) ...
        full_nonlinear(x,uL,M,m,lc,I,bc,bp,g),...
        [0 Ts],xL);

    xL = Xstep(end,:)';

    X_L(:,k) = xL;


    %% -----------------------------------------------
    % CONTROLLER 2: NMPC WITHOUT TERMINAL COST
    %% -----------------------------------------------

    tic;

    [Uopt_N,~,exit_N(k)] = fmincon(...
        @(U) nmpc_cost(U,xN,N,Ts,Q,R,...
        M,m,lc,I,bc,bp,g,false,[]),...
        Uguess_N,[],[],[],[],lb,ub,[],fmin_opts);

    time_N(k) = toc;

    if exit_N(k) <= 0
        warning('NMPC without terminal cost: solver issue at step %d.',k);
    end

    if isempty(Uopt_N) || any(~isfinite(Uopt_N))
        error('NMPC without terminal cost returned invalid solution.');
    end

    uN = Uopt_N(1);

    U_N(k) = uN;

    % Shift solution for next optimization
    Uguess_N = [Uopt_N(2:end); Uopt_N(end)];

    [~,Xstep] = ode45(@(t,x) ...
        full_nonlinear(x,uN,M,m,lc,I,bc,bp,g),...
        [0 Ts],xN);

    xN = Xstep(end,:)';

    X_N(:,k) = xN;


    %% -----------------------------------------------
    % CONTROLLER 3: NMPC WITH TERMINAL COST
    %% -----------------------------------------------

    tic;

    [Uopt_T,~,exit_T(k)] = fmincon(...
        @(U) nmpc_cost(U,xT,N,Ts,Q,R,...
        M,m,lc,I,bc,bp,g,true,P_lqr),...
        Uguess_T,[],[],[],[],lb,ub,[],fmin_opts);

    time_T(k) = toc;

    if exit_T(k) <= 0
        warning('Terminal NMPC: solver issue at step %d.',k);
    end

    if isempty(Uopt_T) || any(~isfinite(Uopt_T))
        error('Terminal NMPC returned invalid solution.');
    end

    uT = Uopt_T(1);

    U_T(k) = uT;

    % Shift solution for next optimization
    Uguess_T = [Uopt_T(2:end); Uopt_T(end)];

    [~,Xstep] = ode45(@(t,x) ...
        full_nonlinear(x,uT,M,m,lc,I,bc,bp,g),...
        [0 Ts],xT);

    xT = Xstep(end,:)';

    X_T(:,k) = xT;

end

t = (0:Nsim-1)*Ts;

fprintf('\nSimulation completed.\n');

%% =====================================================
% 8. STATE COMPARISON PLOTS
%% =====================================================

figure('Name','Three MPC Controllers',...
    'Color','w','Position',[100 50 1100 900]);

labels = {'Cart Position [m]',...
          'Cart Velocity [m/s]',...
          'Pendulum Angle [deg]',...
          'Angular Velocity [deg/s]',...
          'Control Force [N]'};

YL = {X_L(1,:),X_L(2,:),...
      X_L(3,:)*180/pi,X_L(4,:)*180/pi,U_L};

YN = {X_N(1,:),X_N(2,:),...
      X_N(3,:)*180/pi,X_N(4,:)*180/pi,U_N};

YT = {X_T(1,:),X_T(2,:),...
      X_T(3,:)*180/pi,X_T(4,:)*180/pi,U_T};

for i = 1:5

    subplot(5,1,i);

    plot(t,YL{i},'b','LineWidth',1.5);
    hold on;

    plot(t,YN{i},'r--','LineWidth',1.5);

    plot(t,YT{i},'k-.','LineWidth',1.5);

    ylabel(labels{i});
    grid on;

    if i == 1
        title('Linear MPC vs NMPC vs NMPC with Terminal Cost');

        legend('Linear MPC',...
               'NMPC without terminal cost',...
               'NMPC with terminal cost',...
               'Location','best');
    end

    if i == 5
        xlabel('Time [s]');
    end

end

%% =====================================================
% 9. PHASE PORTRAIT
%% =====================================================

figure('Color','w');

plot(X_L(3,:)*180/pi,...
     X_L(4,:)*180/pi,'b','LineWidth',1.5);

hold on;

plot(X_N(3,:)*180/pi,...
     X_N(4,:)*180/pi,'r--','LineWidth',1.5);

plot(X_T(3,:)*180/pi,...
     X_T(4,:)*180/pi,'k-.','LineWidth',1.5);

xlabel('Pendulum Angle [deg]');
ylabel('Angular Velocity [deg/s]');

title('Pendulum Phase Portrait');

legend('Linear MPC','NMPC','NMPC with terminal cost');
grid on;

%% =====================================================
% 10. COMPUTATION TIME COMPARISON
%% =====================================================

figure('Color','w');

bar([mean(time_L),mean(time_N),mean(time_T)]);

set(gca,'XTickLabel',...
    {'Linear MPC','NMPC','NMPC + Terminal'});

ylabel('Average Optimization Time [s]');

title('Computational Time Comparison');
grid on;

%% =====================================================
% 11. PERFORMANCE METRICS
%% =====================================================

% Angle IAE (rad.s)
IAE_angle = Ts*[sum(abs(X_L(3,:))),...
                sum(abs(X_N(3,:))),...
                sum(abs(X_T(3,:)))];

% Position IAE (m.s)
IAE_position = Ts*[sum(abs(X_L(1,:))),...
                   sum(abs(X_N(1,:))),...
                   sum(abs(X_T(1,:)))];

% Control energy (N^2.s)
Energy = Ts*[sum(U_L.^2),...
             sum(U_N.^2),...
             sum(U_T.^2)];

% Maximum force
MaxForce = [max(abs(U_L)),...
            max(abs(U_N)),...
            max(abs(U_T))];

% Average computation time
MeanTime = [mean(time_L),...
            mean(time_N),...
            mean(time_T)];

% Settling time:
% Both |theta| <= 1 degree and |cart position| <= 0.02 m
Settling = [settling_time(X_L,t),...
            settling_time(X_N,t),...
            settling_time(X_T,t)];

Results = table(IAE_angle',IAE_position',Energy',...
    MaxForce',MeanTime',Settling',...
    'VariableNames',{'IAE_Angle','IAE_Position',...
    'Control_Energy','Maximum_Force',...
    'Mean_Solve_Time','Settling_Time'},...
    'RowNames',{'Linear_MPC','NMPC','NMPC_Terminal'});

fprintf('\n========== PERFORMANCE RESULTS ==========\n');
disp(Results);

fprintf('\nNMPC solver non-positive exit flags: %d\n',...
    sum(exit_N<=0));

fprintf('Terminal NMPC non-positive exit flags: %d\n',...
    sum(exit_T<=0));

%% =====================================================
% 12. SAVE RESULTS
%% =====================================================

save('Three_MPC_Comparison_Results.mat',...
    't','X_L','X_N','X_T','U_L','U_N','U_T',...
    'time_L','time_N','time_T','Results',...
    'exit_N','exit_T');

writetable(Results,...
    'Three_MPC_Performance_Comparison.csv',...
    'WriteRowNames',true);

fprintf('\nResults saved successfully.\n');

%% =====================================================
% LOCAL FUNCTIONS
%% =====================================================

function dx = full_nonlinear(x,u,M,m,lc,I,bc,bp,g)

    theta = x(3);
    th_dot = x(4);
    xd = x(2);

    Mmat = [M+m, m*lc*cos(theta);
            m*lc*cos(theta), I];

    % Upright equilibrium convention
    f = [u-bc*xd-m*lc*th_dot^2*sin(theta);
         m*g*lc*sin(theta)-bp*th_dot];

    acc = Mmat\f;

    dx = [xd;acc(1);th_dot;acc(2)];

end


function J = nmpc_cost(U,x0,N,Ts,Q,R,...
    M,m,lc,I,bc,bp,g,use_terminal,P)

    x = x0;
    J = 0;

    for i = 1:N

        u = U(i);

        % RK4 prediction
        k1 = full_nonlinear(x,u,M,m,lc,I,bc,bp,g);

        k2 = full_nonlinear(x+0.5*Ts*k1,...
            u,M,m,lc,I,bc,bp,g);

        k3 = full_nonlinear(x+0.5*Ts*k2,...
            u,M,m,lc,I,bc,bp,g);

        k4 = full_nonlinear(x+Ts*k3,...
            u,M,m,lc,I,bc,bp,g);

        x = x+(Ts/6)*(k1+2*k2+2*k3+k4);

        % Stage cost
        J = J+x'*Q*x+R*u^2;

    end

    % Optional terminal cost
    if use_terminal

        J = J+x'*P*x;

    end

end


function Ts_settle = settling_time(X,t)

    angle_limit = 1*pi/180;
    position_limit = 0.02;

    settled = abs(X(3,:)) <= angle_limit & ...
              abs(X(1,:)) <= position_limit;

    Ts_settle = NaN;

    for i = 1:length(t)

        if all(settled(i:end))

            Ts_settle = t(i);
            return;

        end

    end

end
