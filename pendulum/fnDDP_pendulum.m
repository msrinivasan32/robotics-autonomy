%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Differential Dynamic Programming               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Course: Robotics and Autonomy                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%  AE4803  Fall  2018                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Author: Evangelos Theodorou                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Uses the following fuctions:
% fnState_and_Control_Transition_Matrices_1
% fnsimulate_1
% fnCost_1
% fnCostComputation

function [u_k] = fnDDP_pendulum(xo, p_target, Q_f, R, T, dt, gamma, num_iter)

global I
global b
global m
global g
global l

b = 1;
m = 1;
g = 9.81;
l = 1;
I = m.*l.^2;


% Time horizon parameters
horizon = T;

% Initial Configuration:
u_k = zeros(1,horizon-1); % control
du_k = zeros(1,horizon-1); % control
x_traj = zeros(2,horizon); % trajectory

% Target Configuration: 
% Given by p_target input

%% DIFFERENTIAL DYNAMIC PROGRAMMING
for k = 1:num_iter
    
    % Computing the final value function values
    Vxx(:,:,horizon)= Q_f;
    Vx(:,horizon) = Q_f * (x_traj(:,horizon) - p_target); 
    V(horizon) = 0.5 * (x_traj(:,horizon) - p_target)' * Q_f * (x_traj(:,horizon) - p_target); 

    % Loop for the quad approx of cost function and linearization of dynamics
    for  j = 1:(horizon-1)

        % Quadratic approximations of the cost function
        [l0,l_x,l_xx,l_u,l_uu,l_ux, l_xu] = fnCost_1(x_traj(:,j), u_k(:,j), j,R,dt);
        % Breaking up the pieces of the appromix for ease of access later
        q0(j) = dt * l0;
        q_k(:,j) = dt * l_x;
        Q_k(:,:,j) = dt * l_xx;
        r_k(:,j) = dt * l_u;
        R_k(:,:,j) = dt * l_uu;
        P_k(:,:,j) = dt * l_ux;
        S_k(:,:,j) = dt * l_xu;

        % Linearization of the dynamics
        [dfx,dfu] = fnState_And_Control_Transition_Matrices_1(x_traj(:,j),u_k(:,j));
        
        % Phi and Beta from the dx(t_k+1) equation
        phi(:,:,j) = eye(2) + dfx * dt;
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
    dx = zeros(2,1);
    for i=1:(horizon-1)    
       du = -inv(Quu(:,:,i)) * (Qux(:,:,i)*dx+Qu(:,i));
       dx = phi(:,:,i) * dx + beta(:,:,i) * du;  
       u_new(:,i) = u_k(:,i) + gamma * du;
    end

    % Updating u to the new value
    u_k = u_new;

    % Simulation of the Nonlinear System (Forward pass)
	[x_traj] = fnsimulate_1(xo,u_new,horizon,dt,0);
	[cost(:,k)] =  fnCostComputation(x_traj,u_k,p_target,dt,Q_f,R);
    
    % Printing to the console for debugging
    %fprintf('iLQG Iteration %d,  Current Cost = %e, Control input = %d \n',k,cost(1,k));

 
end

end