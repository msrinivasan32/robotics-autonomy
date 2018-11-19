function [x,xd] = fnsimulate_2(xo,u_new,Horizon,dt,uncert,sigma)

global mp;
global mc;
global g;
global l;

x = xo;
xd= xo;
for k = 1:(Horizon-1)
    
    coeff1 = 1/(mc*uncert + mp*uncert.*sin(x(3,k)).^2);
    coeff2 = coeff1./l;

    Fx(1,1) = x(2,k);
	Fx(2,1) = coeff1.*(mp*uncert.*sin(x(3,k)).*(l.*x(4,k).^2 + g.*cos(x(3,k))));
	Fx(3,1) = x(4,k);
	Fx(4,1) = coeff2.*(-mp*uncert.*l*uncert.*x(4,k).^2.*cos(x(3,k)).*sin(x(3,k)) - (mc*uncert+mp*uncert).*g.*sin(x(3,k)));
    
    G_x(1,1) = 0;
    G_x(2,1) = coeff1;
    G_x(3,1) = 0;
    G_x(4,1) = -cos(x(3,k)).*coeff2;

    x(:,k+1) = x(:,k) + Fx * dt + G_x * u_new(:,k) * dt  + G_x * u_new(:,k) * sqrt(dt) * randn * sigma;
    xd(:,k+1) = xd(:,k) + Fx + G_x * u_new(:,k) + G_x * u_new(:,k) * sqrt(dt) * randn * sigma;
    
end