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
    t = [0:dt:75]';
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
P = [B, A*B, (A^2)*B, (A^3)*B, (A^4)*B, (A^5)*B];
P_rank = rank(P);
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
G = zeros(6,6);
for i =1:length(t)
    G = G + dt*(expm(A*(-t(i)))*B*B'*expm(A'*(-t(i))));
end
for i = 1:length(Eig_val)
    zeta(:,i) = expm(-A*(t(end)-t(1)))*x_final - Q(:,i);
    E(i) = zeta(:,i)'*inv(G)*zeta(:,i);
    v(:,i) = G\zeta(:,i);
    for j = 1:length(t)
        u_t{i}(j,:) = [B'*expm(-A'*(t(j)-t(1)))*v(:,i)]';
    end
end

% Print minimum energy for each mode (with units of N⋅m²⋅s):
fprintf('\nMinimum Energy for Each Mode:\n');
for i = 1:6
    fprintf('Mode %d: %.2e\n', i, E(i));
end

%problem 2 plot parameter
plot_2 = 0;

% Simulate the system response
for i = 1:6
    x0{i} = Q(:,i);
    t = [0:0.01:75]';
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
desired_eigs = [-.1 -.1 -.1 -.2 -.2 -.2]';  % First 3 for quaternions (τ=1s), last 3 for rates (τ=0.5s)

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
plot_4 = 0;

% Define desired eigenvalues with slight offsets to avoid multiplicity
% Order them to have angular velocities first since they are directly controlled, and will naturally appear first.
% This doesn't change how they relate to states, just represents the order of the modes.
desired_eigs = [-0.2 -0.201 -0.2001 -0.1 -0.1001 -0.1004]';  % Rates first, then quaternions

% Compute controllability matrix and verify rank
Co = ctrb(A, B);
rank_Co = rank(Co);
fprintf('Controllability matrix rank: %d\n', rank_Co);

% Try pole placement with ordered poles
K = place(A, B, desired_eigs);

% Create closed-loop system
Acl = A - B*K;
sys_cl = ss(Acl, B, C, D);

% Verify eigenvalues and eigenvectors
[V_cl, D_cl] = eig(Acl);
[eigvals, idx] = sort(diag(D_cl), 'ComparisonMethod', 'real');  % Sort by real part
V_cl = V_cl(:,idx);

fprintf('\nFinal closed-loop eigenvalues:\n')
disp(eigvals)

% Time vector for simulation
t = 0:0.001:75;

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
    G = zeros(6,6);
    for j = 1:length(t)
        G = G + 0.01*(expm(A*(-t(j)))*B*B'*expm(A'*(-t(j))));
    end
    zeta = expm(-A*(t(end)-t(1)))*zeros(6,1) - x0;
    v = G\zeta;
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

% First verify pole placement worked
fprintf('Desired eigenvalues:\n')
disp(desired_eigs)

% Now look at modal responses
if plot_4 == 1
    figure('Position', [100 100 1200 600])
    sgtitle('Modal Response Verification')
    
    % Get modal transformation and sort by eigenvalue magnitude
    [V_cl, D_cl] = eig(Acl);
    
    % Initialize storage for modal responses
    modal_resp = zeros(6, length(t));
    
    % For each initial condition
    for i = 1:6
        x0 = Q(:,i);
        z0 = inv(V_cl)*x0;
        
        % For each mode, compute response
        for j = 1:6
            z = z0(j)*exp(D_cl(j,j)*t);
            modal_resp(j,:) = max(abs(modal_resp(j,:)), abs(z));
        end
    end
    
    % Normalize responses
    for i = 1:6
        if modal_resp(i,1) > 0
            modal_resp(i,:) = modal_resp(i,:) / modal_resp(i,1);
        end
    end
    
    % Plot with distinct colors
    hold on
    % Plot quaternion modes
    for i = 1:3
        plot(t, modal_resp(i,:), 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
    end
    % Plot angular rate modes
    for i = 4:6
        plot(t, modal_resp(i,:), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
    end
    
    % Plot expected decays
    plot(t, exp(-0.1*t), 'r--', 'LineWidth', 2)  % τ=10s for quaternions
    plot(t, exp(-0.2*t), 'k--', 'LineWidth', 2)  % τ=5s for angular rates
    
    xlabel('Time (s)')
    ylabel('Normalized Modal Response Magnitude')
    legend({'q1 mode', 'q2 mode', 'q3 mode', 'ω1 mode', 'ω2 mode', 'ω3 mode', ...
            'Expected τ=10s', 'Expected τ=5s'}, 'Location', 'southwest')
    set(gca, 'YScale', 'log')
    grid on
    ylim([0.01 1])
    xlim([0 max(t)])

    % save figure
    saveas(gcf, '../figures/problem4_modal_response.png')
end

% Calculate energies for each mode
mode_names = {'ω1', 'ω2', 'ω3', 'q1', 'q2', 'q3'};
energies = struct();

% Calculate controllability Gramian W
G = zeros(6,6);
for i = 1:length(t)
    G = G + dt*(expm(A*(-t(i)))*B*B'*expm(A'*(-t(i))));
end

% Calculate energies for each mode
for i = 1:6
    x0 = Q(:,i);  % Initial condition from orthonormal basis
    
    % Theoretical open loop minimum energy
    zeta = expm(-A*(t(end)-t(1)))*zeros(6,1) - x0;
    energies.theoretical(i) = zeta'*inv(G)*zeta;
    
    % Actual feedback energy used
    energies.feedback(i) = trapz(t, sum(u_cl_all{i}.^2, 1));
end

% Plot parameter for energy comparison
plot_energy = 0;

% Create simple bar plot comparing energies
if plot_energy == 1
    figure('Name', 'Control Energy Comparison', 'Position', [100 100 800 400]);
    bar([energies.theoretical; energies.feedback]');
    set(gca, 'YScale', 'log');
    set(gca, 'XTickLabel', mode_names);
    xlabel('Mode');
    ylabel('Control Energy (N⋅m)²⋅s');
    legend('Open Loop Minimum', 'Feedback');
    title('Control Energy Comparison by Mode');
    grid on;

    % save figure
    saveas(gcf, '../figures/problem4_energy_comparison.png')
end

%% Problem 5
% (20 pts)
% Construct an input gain matrix F such that reference inputs are accurately tracked by the corresponding outputs at low frequency.
% Simulate the response of the closed loop system to a unit step input (one input at a time if there are multiple inputs).
% Discuss the behavior of the step response relative to the closed loop state matrix eigenvalues.

% Plot parameter for problem 5
plot_5 = 0;

% For reference tracking, we need to ensure that the steady-state output
% matches the reference input for the quaternion components (q1, q2, q3)

% Select quaternion outputs for tracking (first 3 states)
C_track = C(1:3,:);  % Take first 3 rows (q1, q2, q3)
D_track = D(1:3,:);  % Take first 3 rows

% Calculate F using steady-state tracking condition
% For step reference tracking: y(∞) = -C(A-BK)^(-1)BF*r = r
% Therefore: F = -(C(A-BK)^(-1)B)^(-1)
Acl = A - B*K;  % Closed-loop A matrix from Problem 4
F = inv(C_track*inv(-Acl)*B);  % Input gain matrix for tracking

% Create augmented closed-loop system with reference input
sys_tracking = ss(Acl, B*F, C, D*F);

% Time vector for simulation
t = 0:0.01:75;

% Initialize plots
if plot_5 == 1
    figure('Position', [100 100 1200 1000])  % Made taller to accommodate extra subplot
    sgtitle('Step Response Analysis')

    % First simulate individual step responses
    for i = 1:3  
        % Create step input (one channel at a time)
        r = zeros(3, length(t));
        r(i,:) = ones(1,length(t));
        
        % Simulate system response
        [y, t, x] = lsim(sys_tracking, r', t, zeros(6,1));
        u = -K*x' + F*r;  % Total control input (feedback + feedforward)
        
        % Plot states
        subplot(4,2,2*i-1)
        plot(t, y(:,1:3), 'LineWidth', 1.5)  % Plot quaternions
        hold on
        plot(t, r(i,:), 'k--', 'LineWidth', 1.5)  % Plot reference
        title(sprintf('Quaternion Response - Step in q_%d', i))
        xlabel('Time (s)')
        ylabel('Quaternion Value')
        legend('q1', 'q2', 'q3', 'Reference')
        grid on
        
        % Plot control signals
        subplot(4,2,2*i)
        plot(t, u', 'LineWidth', 1.5)
        title(sprintf('Control Signals - Step in q_%d', i))
        xlabel('Time (s)')
        ylabel('Control Torque (N⋅m)')
        legend('\tau_1', '\tau_2', '\tau_3')
        grid on
    end

    % Create sequential step inputs
    r = zeros(3, length(t));
    r(1, t >= 0) = ones(1, sum(t >= 0));  % Step in q1 at t=0
    r(2, t >= 25) = 0.5*ones(1, sum(t >= 25));  % Step in q2 at t=25
    
    % Simulate system response
    [y, t, x] = lsim(sys_tracking, r', t, zeros(6,1));
    u = -K*x' + F*r;  % Total control input (feedback + feedforward)
    
    % Plot states
    subplot(4,2,7)
    plot(t, y(:,1:3), 'LineWidth', 1.5)  % Plot quaternions
    hold on
    plot(t, r(1,:), 'k--', 'LineWidth', 1.5)  % Plot q1 reference
    plot(t, r(2,:), 'r--', 'LineWidth', 1.5)  % Plot q2 reference
    title('Quaternion Response - Combined Step')
    xlabel('Time (s)')
    ylabel('Quaternion Value')
    legend('q1', 'q2', 'q3', 'q1 ref', 'q2 ref')
    grid on
    
    % Plot control signals
    subplot(4,2,8)
    plot(t, u', 'LineWidth', 1.5)
    title('Control Signals - Combined Step')
    xlabel('Time (s)')
    ylabel('Control Torque (N⋅m)')
    legend('\tau_1', '\tau_2', '\tau_3')
    grid on

    % Save figure
    %saveas(gcf, '../figures/problem5_step_response.png')
end

% Calculate and display steady-state errors
fprintf('\nSteady-State Analysis:\n')
for i = 1:3
    r = zeros(3,1);
    r(i) = 1;
    
    % Calculate theoretical steady state
    x_ss = -inv(Acl)*B*F*r;
    y_ss = C_track*x_ss;
    
    fprintf('\nStep in q%d:\n', i)
    fprintf('Steady-state error in q%d: %.2e\n', i, abs(y_ss(i) - 1))
    for j = 1:3
        if j ~= i
            fprintf('Steady-state coupling to q%d: %.2e\n', j, abs(y_ss(j)))
        end
    end
end

% Analysis of step response vs eigenvalues
fprintf('\nStep Response Analysis:\n')
fprintf('Rise time corresponds to faster eigenvalues (%.2f, %.2f, %.2f)\n', ...
    desired_eigs(1), desired_eigs(2), desired_eigs(3))
fprintf('Settling time dominated by slower eigenvalues (%.2f, %.2f, %.2f)\n', ...
    desired_eigs(4), desired_eigs(5), desired_eigs(6))

%% Project 2 P1

%ANALYSIS PROBABLY NEEDS TO BE ON A NOT Acl, update
O = [];
for i = 1:size(Acl,1)
    O = [O;C*(Acl^(i-1))];
end
O_rank = rank(O);
O_RN = null(Acl);
%The rank of the observability matrix O is 6=n, therefore all x within R6
%are observable. The RN of Acl is trivial, containing only the zero vector.

%observability grammian
OG = zeros(6,6);
for i =1:length(t)
    OG = OG + dt*(expm(Acl*(t(i)))*C*C'*expm(Acl'*(t(i))));
end
[OG_eigvec,OG_eigval] = eigs(OG);
%small eigvals of OG, less observable states

%energy in output signal from t0 to t1, each modal space unit ICs
for i = 1:length(Eig_val)
    E_y(i) = Q(:,i)'*OG*Q(:,i);
end
%All modes can be changed in an observer designed around the plant dynamics
%as all modes in the plant system are observable? (need to verify)

%% Project 2 P2
%Luenberger observer slow
desired_eigs_L_slow = -desired_eigs/5;
L_slow = [place(A',C',desired_eigs_L_slow)]';
%Luenberger observer equal
desired_eigs_L_eq = -desired_eigs;
%[-0.1001;-0.1002;-0.1003;-0.2001;-0.2002;-0.2003];
L_eq = [place(A',C',desired_eigs_L_eq)]';
%Luenberger observer fast 
desired_eigs_L_fast = -5*desired_eigs;
L_fast = [place(A',C',desired_eigs_L_fast)]';

%% Project 2 P3
%Luenberger ss matrices
AL_slow = [Acl,B*K;zeros(6,6),L_slow];
AL_equal = [Acl,B*K;zeros(6,6),L_eq];
AL_fast = [Acl,B*K;zeros(6,6),L_fast];
BL = [B*F;zeros(6,3)];

%L_slow: set integration time, ICs, ODE call, and plots
slowplot = 1;
if slowplot == 1
t = 0:0.01:600;

x0_test = zeros(6,1);
x0hat_test = zeros(6,1);
for i = 1:3
    r = zeros(3,1);
    r(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_slow,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end

%L_eq: set integration time, ICs, ODE call, and plots
eqplot = 0;
if eqplot == 1
t = 0:0.01:400;

r = zeros(3,1);
x0hat_test = zeros(6,1);
for i = 1:6
    x0_test = zeros(6,1);
    x0_test(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_eq,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end

%L_fast: set integration time, ICs, ODE call, and plots
fastplot = 1;
if fastplot == 1
t = 0:0.01:200;

r = zeros(3,1);
x0hat_test = zeros(6,1);
for i = 1:6
    x0_test = zeros(6,1);
    x0_test(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_fast,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end

%% Project 2 P4
%L_slow: set integration time, ICs, ODE call, and plots
slowplot = 1;
if slowplot == 1
t = 0:0.01:600;

r = zeros(3,1);
x0hat_test = zeros(6,1);
for i = 1:6
    x0_test = zeros(6,1);
    x0_test(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_slow,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end

%L_eq: set integration time, ICs, ODE call, and plots
eqplot = 0;
if eqplot == 1
t = 0:0.01:400;

r = zeros(3,1);
x0hat_test = zeros(6,1);
for i = 1:6
    x0_test = zeros(6,1);
    x0_test(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_eq,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end

%L_fast: set integration time, ICs, ODE call, and plots
fastplot = 0;
if fastplot == 1
t = 0:0.01:200;

r = zeros(3,1);
x0hat_test = zeros(6,1);
for i = 1:6
    x0_test = zeros(6,1);
    x0_test(i) = 1;
    e0_test = x0_test - x0hat_test;
    [ttest{i},xtest{i}] = ode45(@(t,x) CL_Luenberger(t,x,r,L_fast,Acl,B,K,F),t,[x0_test;e0_test]);
    xhattest{i} = xtest{i}(:,1:6)-xtest{i}(:,7:12);

    figure ()
    plot(ttest{i},xtest{i}(:,1),'k')
    hold on
    plot(ttest{i},xtest{i}(:,2),'k')
    plot(ttest{i},xtest{i}(:,3),'k')
    plot(ttest{i},xhattest{i}(:,1),'r')
    plot(ttest{i},xhattest{i}(:,2),'r')
    plot(ttest{i},xhattest{i}(:,3),'r')
    hold off

    figure ()
    plot(ttest{i},xtest{i}(:,4),'k')
    hold on
    plot(ttest{i},xtest{i}(:,5),'k')
    plot(ttest{i},xtest{i}(:,6),'k')
    plot(ttest{i},xhattest{i}(:,4),'r')
    plot(ttest{i},xhattest{i}(:,5),'r')
    plot(ttest{i},xhattest{i}(:,6),'r')
    hold off
end
end


