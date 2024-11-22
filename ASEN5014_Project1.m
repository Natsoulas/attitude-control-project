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
%eigvecs and eigvals
[Eig_vec,Eig_val] = eigs(A);
Eig_vec = [Eig_vec(:,1:3),[0;0;0;2;0;0],[0;0;0;0;2;0],[0;0;0;0;0;2]];

%problem 1 plot parameter
plot_1 = 0;

%setup initial conditions and system response
dt = 0.01;
ss_system = ss(A,B,C,D);
for i = 1:6
    x0{i} = Eig_vec(:,i);
    t = [0:dt:10]';
    u = ones(length(t),1).*[0,0,0];
    [yout{i},tout{i}] = lsim(ss_system,u,t,x0{i});

    if plot_1 == 1
        %quaternion component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,1),'r')
        xlabel('time')
        ylabel('q1')
        title('Quaternion Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,2),'k')
        xlabel('time')
        ylabel('q2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,3),'g')
        xlabel('time')
        ylabel('q3')

        %angular velocity component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,4),'r')
        xlabel('time')
        ylabel('w1')
        title('Angular Velocity Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,5),'k')
        xlabel('time')
        ylabel('w2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,6),'g')
        xlabel('time')
        ylabel('w3')
    end
end

%% Problem 2
%reachability
U = [B, A*B, (A^2)*B, (A^3)*B, (A^4)*B, (A^5)*B];
U_rank = rank(U);
%reachable subspace of the state space is R6

%orthonormal basis Q for R6
[Q,R] = qr(Eig_vec);

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

%problem 1 plot parameter
plot_2 = 1;

for i = 1:6
    x0{i} = Q(:,i);
    t = [0:0.01:10]';
    u = u_t{i};
    [yout{i},tout{i}] = lsim(ss_system,u,t,x0{i});

    if plot_2 == 1
        %quaternion component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,1),'r')
        xlabel('time')
        ylabel('q1')
        title('Quaternion Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,2),'k')
        xlabel('time')
        ylabel('q2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,3),'g')
        xlabel('time')
        ylabel('q3')

        %angular velocity component plots
        figure()
        subplot(3,1,1)
        plot(tout{i},yout{i}(:,4),'r')
        xlabel('time')
        ylabel('w1')
        title('Angular Velocity Components')
        subplot(3,1,2)
        plot(tout{i},yout{i}(:,5),'k')
        xlabel('time')
        ylabel('w2')
        subplot(3,1,3)
        plot(tout{i},yout{i}(:,6),'g')
        xlabel('time')
        ylabel('w3')
    end
end