function[mphat,mchat,lhat]=fnLS(x_traj,u_k,dt)
global g

theta=x_traj(3,:);
thetadot=x_traj(4,:);
xdd= x_traj(4,:)/dt;
z= u_k;
phi= [xdd.*(sin(theta)).^2-g.*sin(theta).*cos(theta); xdd; -sin(theta).*(thetadot).^2];
w = inv(transpose(phi)*phi)*transpose(phi)*z;

mphat=w(1);
mchat=w(2);
lhat =w(3);


end
