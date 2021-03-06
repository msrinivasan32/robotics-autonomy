 clear;
 clc;
 close all;
 
 % Time initialization
 horizon = 300;
 dt = 0.01;
 T = 10; % shorter DDP time horizon
 
 % Trajectory initialization
 xo = zeros(2,1); % state: [theta; theta_dot]
 x_traj = zeros(2,horizon); % trajectory
 x_traj_different=zeros(2,horizon);
 xo_different=zeros(2,1);
 % Final configuration
 p_target(1,1) = pi;
 p_target(2,1) = 0;
 
 % Weight in Final State:
 Q_f(1,1) = 400;
 Q_f(2,2) = 100;

% Weight in the Control:
 R = 1;
 
 % Learning Rate:
gamma = 1;

% Number of iterations
num_iter = 100;

% Changing the nominal model from the real one
uncert = 1;
otheruncert= 1.0001;
 for i = 1:horizon-1
    
     [u] = fnDDP_pendulum(xo,p_target,Q_f,R,T,dt,gamma,num_iter,1);
     [x_traj_next] = fnsimulate_1(xo,u(1),2,dt,0,uncert);
     x_traj(:,i+1) = x_traj_next(:,end);
     xo = x_traj_next(:,end);
     
     
     [u_different]=fnDDP_pendulum(xo,p_target,Q_f,R,T,dt,gamma,num_iter,1);
     [x_traj_next_different] = fnsimulate_1(xo_different,u(1),2,dt,0,otheruncert);
      x_traj_different(:,i+1) = x_traj_next_different(:,end);
      xo_different = x_traj_next_different(:,end);
      
     fprintf('iLQG Horizon %d,  Current State %d %d,   Control applied %d \n',i, xo(1), xo(2), u(1));
 end  
 
 %% PLOTTING

% Setting up a time vector to span the horizon to plot against
time(1)=0;
for i= 2:horizon
time(i) =time(i-1) + dt;  
end

% Plots
figure(1)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

subplot(2,2,1)
hold on
plot(time,x_traj(1,:),'b', 'linewidth', 1.5) 
plot(time,p_target(1,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta')
xlabel('Time in sec')
hold off
grid

subplot(2,2,2)
hold on
plot(time,x_traj(2,:),'b', 'linewidth', 1.5);
plot(time,p_target(2,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta dot')
xlabel('Time in sec')
hold off
grid

subplot(2,2,3)
hold on
plot(time,x_traj_different(1,:),'b', 'linewidth', 1.5) 
plot(time,p_target(1,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta')
xlabel('Time in sec')
hold off
grid

subplot(2,2,4)
hold on
hold on
plot(time,x_traj_different(2,:),'b', 'linewidth', 1.5);
plot(time,p_target(2,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta dot')
xlabel('Time in sec')
hold off
grid

