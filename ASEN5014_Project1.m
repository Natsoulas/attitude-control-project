%% ASEN 5014 Project 1
clc
clear all
close all

%given inertia matrix
I = [10,0,0;0,100,0;0,0,50]; %kg*m^2

%ss matrices
A = [zeros(3,3),(eye(3,3)*0.5);zeros(3,6)];
B = [zeros(3,3);inv(I)];
C = eye(6,6);
D = zeros(6,3);

%% Problem 1
%eigvectors and eigvalues
[Eig_vec,Eig_val] = eigs(A);

% Augment the eigvectors with the additional generalized eigenvectors
% (Matlab can't solve for them so we have to do it manually)
Eig_vec = [Eig_vec(:,1:3),[0;0;0;2;0;0],[0;0;0;0;2;0],[0;0;0;0;0;2]];

%setup initial conditions and system response
dt = 0.01;
ss_system = ss(A,B,C,D);

% Full System Stability Analysis with Jordan Form
% Check eigenvalues
[V, Diag] = eig(A);
disp('Eigenvalues:');
disp(diag(Diag));

% Analyze Jordan form
[Vj, J] = jordan(A);
disp('Jordan Form:');
disp(J);

disp('Vj:');
disp(Vj);

% Exponential solution (if needed)
syms t
Phi = expm(A*t); % Matrix exponential
disp('State Transition Matrix (Phi):');
disp(Phi);

%problem 1 plot parameter
plot_1 = 0;

if plot_1 == 1
    % Create figures directory if it doesn't exist
    if ~exist('../figures', 'dir')
        mkdir('../figures');
    end

    % Create two figures - one for quaternions, one for angular velocities
    figure('Position', [100 100 1200 800])
    sgtitle('Quaternion Components for All Eigenvector Initial Conditions')
    figure('Position', [100 100 1200 800])
    sgtitle('Angular Velocity Components for All Eigenvector Initial Conditions')
end

% Create descriptive labels for each eigenvector
eig_labels = {'IC: [1 0 0 0 0 0]', 'IC: [0 1 0 0 0 0]', 'IC: [0 0 1 0 0 0]', ...
              'IC: [0 0 0 2 0 0]', 'IC: [0 0 0 0 2 0]', 'IC: [0 0 0 0 0 2]'};

for i = 1:6
    x0{i} = Eig_vec(:,i);
    t = [0:dt:10]';
    u = ones(length(t),1).*[0,0,0];
    [yout{i},tout{i}] = lsim(ss_system,u,t,x0{i});

    if plot_1 == 1
        %quaternion component plots
        figure(1)
        subplot(6,3,(i-1)*3+1)
        plot(tout{i},yout{i}(:,1),'r', 'LineWidth', 1.5)
        ylabel('q1')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end
        
        subplot(6,3,(i-1)*3+2)
        plot(tout{i},yout{i}(:,2),'k', 'LineWidth', 1.5)
        ylabel('q2')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end
        
        subplot(6,3,(i-1)*3+3)
        plot(tout{i},yout{i}(:,3),'g', 'LineWidth', 1.5)
        ylabel('q3')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end

        %angular velocity component plots
        figure(2)
        subplot(6,3,(i-1)*3+1)
        plot(tout{i},yout{i}(:,4),'r', 'LineWidth', 1.5)
        ylabel('ω1')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end
        
        subplot(6,3,(i-1)*3+2)
        plot(tout{i},yout{i}(:,5),'k', 'LineWidth', 1.5)
        ylabel('ω2')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end
        
        subplot(6,3,(i-1)*3+3)
        plot(tout{i},yout{i}(:,6),'g', 'LineWidth', 1.5)
        ylabel('ω3')
        title(eig_labels{i})
        if i == 6, xlabel('time'), end

        % Save quaternion figure
        saveas(gcf, '../figures/problem1_quaternions.png')
        
        % Save angular velocity figure
        saveas(gcf, '../figures/problem1_angular_velocities.png')
    end
end

%% Problem 2
%reachability
U = [B, A*B, (A^2)*B, (A^3)*B, (A^4)*B, (A^5)*B];
U_rank = rank(U);
%reachable subspace of the state space is R6

%orthonormal basis Q for R6
[Q,R] = qr(Eig_vec);

% print out the orthonormal basis Q
disp('Orthonormal Basis Q:');
disp(Q);

%{
Determine the energy required by the minimum-energy control to restore 
the state to the equilibrium at zero from unit-perturbations in each of 
the orthonormal basis vector directions for the reachable subspace.
%}
x_final = [0;0;0;0;0;0];
W = zeros(6,6);
for i =1:length(t)
    W = W + dt*(expm(A*(-t(i)))*B*B'*expm(A'*(-t(i))));
end
for i = 1:length(Eig_val)
    zeta(:,i) = expm(-A*(t(end)-t(1)))*x_final - Q(:,i);
    E(i) = zeta(:,i)'*inv(W)*zeta(:,i);
    v(:,i) = W\zeta(:,i);
    for j = 1:length(t)
        u_t{i}(j,:) = [B'*expm(-A'*(t(j)-t(1)))*v(:,i)]';
    end
end

%problem 2 plot parameter
plot_2 = 0;

for i = 1:6
    x0{i} = Q(:,i);
    t = [0:0.01:10]';
    u = u_t{i};
    [yout{i},tout{i}] = lsim(ss_system,u,t,x0{i});

    if plot_2 == 1
        %quaternion component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,1),'r', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('q1')
        title('Quaternion Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,2),'k', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('q2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,3),'g', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('q3')

        %angular velocity component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,4),'r', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('w1')
        title('Angular Velocity Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,5),'k', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('w2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,6),'g', 'LineWidth', 1.5)
        xlabel('time')
        ylabel('w3')
    end
end

%% Problem 3

% First, let's properly check controllability of each mode using the PBH test
% PBH Test: For each eigenvalue λ, check rank([λI - A, B]) = n

% Get system eigenvalues
[V, Diag] = eig(A);
eigenvalues = diag(Diag);
n = size(A,1);  % system order

fprintf('PBH Test Results:\n');
for i = 1:length(eigenvalues)
    lambda = eigenvalues(i);
    PBH_matrix = [lambda*eye(n) - A, B];
    PBH_rank = rank(PBH_matrix);
    fprintf('Mode %d (λ = %.3f): rank = %d (of %d possible)\n', ...
            i, lambda, PBH_rank, n);
    if PBH_rank == n
        fprintf('  → Mode is controllable\n');
    else
        fprintf('  → Mode is NOT controllable\n');
    end
end

% Define desired closed-loop eigenvalues
desired_eigs = [
    -0.1;  % Slower mode for q1 (roll angle)
    -0.1;  % Slower mode for q2 (pitch angle)
    -0.1;   % Slower mode for q3 (yaw angle)
    -0.2;  % Fast mode for w1 (roll rate)
    -0.2;  % Fast mode for w2 (pitch rate)
    -0.2   % Fast mode for w3 (yaw rate)
];

% Display time constants
fprintf('\nTime Constants:\n');
fprintf('Quaternions (q1,q2,q3): %.1f seconds\n', -1/desired_eigs(1));
fprintf('Angular Rates (w1,w2,w3): %.1f seconds\n', -1/desired_eigs(4));


%% Problem 4
% (30 pts)
% Design a state-feedback controller to place Eigenvalues (that can be placed) at the
% desired locations. Find the closed loop Eigenvectors and closed loop Modal spaces. Simulate the
% closed loop system response to unit initial state perturbations in each of the orthonormal basis
% directions for the reachable subspace. Compare this to the expected behavior from the choice
% of Eigenvalues. Compare the control signals for state feedback control to the corresponding
% ones (i.e. in corresponding orthogonal basis vector directions) for minimum-energy open loop
% control.

% Plot parameter for problem 4
plot_4 = 1;

% Calculate feedback gain matrix K using place command
K = place(A, B, desired_eigs);

% Create closed-loop system
Acl = A - B*K;
sys_cl = ss(Acl, B, C, D);

% Find closed-loop eigenstructure
[V_cl, D_cl] = eig(Acl);
eigenvalues_cl = diag(D_cl);

% Normalize eigenvectors to get modal matrix
M_cl = zeros(size(V_cl));
for i = 1:size(V_cl,2)
    M_cl(:,i) = V_cl(:,i) / norm(V_cl(:,i));
end

% Display closed-loop eigenvalues and modal matrix
disp('Closed-loop eigenvalues:');
disp(eigenvalues_cl);
disp('Closed-loop eigenvectors:');
disp(V_cl);
disp('Closed-loop modal matrix:');
disp(M_cl);

% Time vector for simulation
t = 0:0.01:100;

% Initialize plots
if plot_4 == 1
    % Create three figures - quaternions, angular rates, and controls
    figure('Position', [100 100 1200 800])
    sgtitle('Quaternion Components for All Initial Conditions')
    
    figure('Position', [100 100 1200 800])
    sgtitle('Angular Velocity Components for All Initial Conditions')
    
    figure('Position', [100 100 1200 800])
    sgtitle('Control Signals for All Initial Conditions')
    
    % Create color map
    colors = lines(3);  % Get 3 colors for 3 components :)
end

% Initialize cell arrays to store all data
y_cl_all = cell(6,1);
u_cl_all = cell(6,1);
u_ol_all = cell(6,1);

% Simulate responses for each basis vector
for i = 1:6
    x0 = Q(:,i);  % Initial condition from orthonormal basis
    
    % Closed-loop simulation
    [y_cl_all{i}, t_cl, x_cl] = initial(sys_cl, x0, t);
    u_cl_all{i} = -K * x_cl';  % Feedback control signal
    
    % Open-loop minimum energy control
    W = zeros(6,6);
    for j = 1:length(t)
        W = W + 0.01*(expm(A*(-t(j)))*B*B'*expm(A'*(-t(j))));
    end
    zeta = expm(-A*(t(end)-t(1)))*zeros(6,1) - x0;
    v = W\zeta;
    u_ol = zeros(3, length(t));
    for j = 1:length(t)
        u_ol(:,j) = B'*expm(-A'*(t(j)-t(1)))*v;
    end
    u_ol_all{i} = u_ol;
    
    if plot_4 == 1
        % Plot quaternions
        figure(1)
        subplot(3,2,i)
        hold on
        for j = 1:3
            plot(t, y_cl_all{i}(:,j), 'Color', colors(j,:), 'LineWidth', 1.5)
        end
        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        ylabel('Quaternions')
        legend('q1', 'q2', 'q3')
        grid on
        
        % Plot angular rates
        figure(2)
        subplot(3,2,i)
        hold on
        for j = 1:3
            plot(t, y_cl_all{i}(:,j+3), 'Color', colors(j,:), 'LineWidth', 1.5)
        end
        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        ylabel('Angular Rates')
        legend('ω1', 'ω2', 'ω3')
        grid on
        
        % Plot control signals
        figure(3)
        subplot(3,2,i)
        hold on
        % Plot feedback control
        for j = 1:3
            plot(t, u_cl_all{i}(j,:), 'Color', colors(j,:), 'LineWidth', 1.5, 'LineStyle', '-')
        end
        % Plot open-loop control
        for j = 1:3
            plot(t, u_ol_all{i}(j,:), 'Color', colors(j,:), 'LineWidth', 1.5, 'LineStyle', '--')
        end
        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        ylabel('Control Torque (N⋅m)')
        legend('FB τ1', 'FB τ2', 'FB τ3', 'OL τ1', 'OL τ2', 'OL τ3')
        grid on
    end
end

if plot_4 == 1
    % Save figures
    figure(1)
    saveas(gcf, '../figures/problem4_quaternions.png')
    figure(2)
    saveas(gcf, '../figures/problem4_angular_velocities.png')
    figure(3)
    saveas(gcf, '../figures/problem4_controls.png')
end
