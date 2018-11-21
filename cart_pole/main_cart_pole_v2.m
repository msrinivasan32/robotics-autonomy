%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE4803  Fall  2018                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Evangelos Theodorou                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Uses the following fuctions:
% fnState_and_Control_Transition_Matrices_2
% fnsimulate_2
% fnCost_2
% fnCostComputation
% fnLS

clear all;
close all;
clc;

%% SET UP GLOBAL VARIABLES

global mp;
global mc;
global g;
global l;

global mphat;
global mchat;
global lhat;

% masses in Kgr
 mp = .01;
 mc = 1;

% length parameter in meters
l = 0.25;

%gravity term in meters/seconds sqaured 
g=9.8;

%initial guesses
mphat = .01;
mchat = 1;
lhat = 0.25;

paramchange(1,1)=mphat;
paramchange(2,1)=mchat;
paramchange(3,1)=lhat;

%% TIME HORIZON AND INITIAL & TARGET STATE

% Time horizon parameters
horizon = 100; % 1.5sec
dt = 0.01; % for discretization

% Initial Configuration:
xo = zeros(4,1); % state: [x; x_dot; theta; theta_dot]
xo(3) = xo(3) + 0.8; % initial perturbation to solve for non-singularity w
u_k = zeros(1,horizon-1); % control
du_k = zeros(1,horizon-1); % control
x_traj = zeros(4,horizon); % trajectory

% Target Configuration: 
p_target(1,1) = 10; 
p_target(2,1) = 0;
p_target(3,1) = pi;
p_target(4,1) = 0;


%% PARAMETERS TO TUNE

% Weight in Final State:
Q_f(1,1) = 1;
Q_f(2,2) = 150;
Q_f(3,3) = 400;
Q_f(4,4) = 50;

% Weight in the Control:
R = 2 * eye(1,1);

% Learning Rate:
gamma = 0.5;

% Number of iterations
num_iter = 50;
plotlength=1:num_iter;

%% DIFFERENTIAL DYNAMIC PROGRAMMING
for k = 1:num_iter
    
    % Computing the final value function values
    Vxx(:,:,horizon)= Q_f;
    Vx(:,horizon) = Q_f * (x_traj(:,horizon) - p_target); 
    V(horizon) = 0.5 * (x_traj(:,horizon) - p_target)' * Q_f * (x_traj(:,horizon) - p_target); 

    % Loop for the quad approx of cost function and linearization of dynamics
    for  j = 1:(horizon-1)

        % Quadratic approximations of the cost function
        [l0,l_x,l_xx,l_u,l_uu,l_ux, l_xu] = fnCost_2(x_traj(:,j), u_k(:,j), j,R,dt);
        % Breaking up the pieces of the appromix for ease of access later
        q0(j) = dt * l0;
        q_k(:,j) = dt * l_x;
        Q_k(:,:,j) = dt * l_xx;
        r_k(:,j) = dt * l_u;
        R_k(:,:,j) = dt * l_uu;
        P_k(:,:,j) = dt * l_ux;
        S_k(:,:,j) = dt * l_xu;

        % Linearization of the dynamics
        [dfx,dfu] = fnState_And_Control_Transition_Matrices_2(x_traj(:,j),u_k(:,j));
        
        % Phi and Beta from the dx(t_k+1) equation
        phi(:,:,j) = eye(4,4) + dfx * dt;
        beta(:,:,j) = dfu * dt; 

    end
    
    % Loop to backpropagate the value function
    for j = (horizon-1):-1:1

        % setting up Q matrices to simplify the DDP code
        Qo(j) = q0(j) + V(j+1);
        Qx(:,j) = q_k(:,j) + phi(:,:,j)' * Vx(:,j+1);
        Qu(:,j) = r_k(:,j) + beta(:,:,j)' * Vx(:,j+1);
        Qxx(:,:,j) = Q_k(:,:,j) + phi(:,:,j)' * Vxx (:,:,j+1) * phi(:,:,j);
        Qxu(:,:,j) = S_k(:,:,j) + phi(:,:,j)' * Vxx(:,:,j+1) * beta(:,:,j);
        Qux(:,:,j) = P_k(:,:,j) + beta(:,:,j)' * Vxx(:,:,j+1) * phi(:,:,j);
        Quu(:,:,j) = R_k(:,:,j) + beta(:,:,j)' * Vxx(:,:,j+1) * beta(:,:,j);

        % Backpropagating the value function
        inv_Quu = inv(Quu(:,:,j));
        Vxx(:,:,j) = Qxx(:,:,j) - Qxu(:,:,j) * inv_Quu * Qux(:,:,j);
        Vx(:,j) = (Qx(:,j)' - Qu(:,j)' * inv_Quu * Qux(:,:,j))';
        V(j) = Qo(j) - 0.5 * Qu(:,j)' * inv_Quu * Qu(:,j);

    end 

    % Determining optimal du and computing dx
    dx = zeros(4,1);
    for i=1:(horizon-1)    
       du = -inv(Quu(:,:,i)) * (Qux(:,:,i)*dx+Qu(:,i));
       dx = phi(:,:,i) * dx + beta(:,:,i) * du;  
       u_new(:,i) = u_k(:,i) + gamma * du;
    end

    % Updating u to the new value
    u_k = u_new;

    % Simulation of the Nonlinear System (Forward pass)
	[x_traj,xd] = fnsimulate_2(xo,u_new,horizon,dt,1,0);
	[cost(:,k)] = fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    
    % Printing to the console for debugging
    fprintf('iLQG Iteration %d,  Current Cost = %e \n',k,cost(1,k));

    % Linear Regression to get new parameters
    M=10;
    x_LS=zeros(4,M);
    xd_LS=zeros(4,M);
    u_LS=zeros(1,M)
    for m=1:M
        u_input=u_new(1)+randn
        [x,xd]=fnsimulate_2(xo,u_input,2,dt,0.01,8)
        x_LS(:,m)=x(:,2);
        xd_LS(:,m)=xd(:,2);
        u_LS(:,m)=u_input
    end
    
    [mphat,mchat,lhat]=fnLS(x_LS,u_LS,xd_LS);
    paramchange(1,k+1)=mphat;
    paramchange(2,k+1)=mchat;
    paramchange(3,k+1)=lhat;
 
end

%% PLOTTING

% Setting up a time vector to span the horizon to plot against
time(1)=0;
for i= 2:horizon
time(i) =time(i-1) + dt;  
end

% TODO: Add ylabels?

figure(1)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

subplot(4,2,1)
hold on
plot(time,x_traj(1,:),'b', 'linewidth', 1.5) 
plot(time,p_target(1,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('X')
xlabel('Time in sec')
hold off
grid

subplot(4,2,2)
hold on
plot(time,x_traj(2,:),'b', 'linewidth', 1.5);
plot(time,p_target(2,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('X dot')
xlabel('Time in sec')
hold off
grid

subplot(4,2,3)
hold on
plot(time,x_traj(3,:),'b', 'linewidth', 1.5);
plot(time,p_target(3,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta')
xlabel('Time in sec')
hold off
grid

subplot(4,2,4)
hold on
plot(time,x_traj(4,:), 'b', 'linewidth', 1.5)
plot(time,p_target(4,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta dot')
xlabel('Time in sec')
hold off
grid

subplot(4,2,5);
hold on
plot(cost, 'g', 'linewidth', 1.5)
xlabel('Iterations')
title('Cost')
hold off
grid

subplot(4,2,6);
hold on
plot(plotlength,mp*ones(1,length(plotlength)),'r','linewidth',1.5)
plot(plotlength,paramchange(1,2:end),'b','linewidth',1.5)
title('Mp')
xlabel('Number of Iterations')
hold off
grid

subplot(4,2,7);
hold on
plot(plotlength,mc*ones(1,length(plotlength)),'r','linewidth',1.5)
plot(plotlength,paramchange(2,2:end),'b','linewidth',1.5)
title('Mc')
xlabel('Number of Iterations')
hold off
grid

subplot(4,2,8)
hold on
plot(plotlength,l*ones(1,length(plotlength)),'r','linewidth',1.5)
plot(plotlength,paramchange(3,2:end),'b','linewidth',1.5)
title('Length')
xlabel('Number of Iterations')
hold off 
grid
