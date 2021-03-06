function [x] = fnsimulate_1(xo,u_new,Horizon,dt,sigma,uncert)

global b
global m
global g
global l
global I

x = xo;

for k = 1:(Horizon-1)
    
    Fx(1,1) = x(2, k);
    Fx(2,1) = -m*uncert*g*l*uncert*sin(x(1,k))/(I*uncert) + -b*uncert/(I*uncert) .* x(2, k);
    
    Gx(1,1) = 0;
    Gx(2,1) = 1/(I*uncert);
    
    x(:,k+1) = x(:,k) + Fx .* dt + Gx .* u_new(:,k) .* dt + sqrt(dt) * sigma * randn;

end


