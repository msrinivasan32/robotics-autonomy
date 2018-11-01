function [dfx,dfu] = fnState_And_Control_Transition_Matrices_1(x,u)

b = 1;
m = 1;
g = 9.81;
l = 1;
I = m.*l.^2;

    x1 = x(1,1);
    x2 = x(2,1);

    dfx = zeros(2);
    dfu = zeros(2,1);

    dfx(1,1) = 0;
    dfx(1,2)  = 1;
    dfx(2,1) = -(m.*g.*l./I.*(cos(x1)));
    dfx(2,2) = -b./I;

    dfu(1,1) = 0;
    dfu(2,1) = 1./I;

end