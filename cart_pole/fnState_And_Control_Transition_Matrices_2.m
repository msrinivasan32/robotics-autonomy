function [A,B] = fnState_And_Control_Transition_Matrices_2(x,u)

global mphat
global mchat
global lhat
global g


    x1 = x(1,1);
    x2 = x(2,1);
    x3 = x(3,1);
    x4 = x(4,1);

    u1 = u(1,1);

    A = zeros(4,4);
    B = zeros(4,1);

    coeff1 = 1/(mchat + mphat.*sin(x3).^2);
    coeff2 = coeff1./lhat;

    A(1,2) = 1;
    A(2,3)=(x4^2*mphat*cos(x3)*(mchat-mphat*(sin(x3)^2))-2*u1*mphat*sin(x3)*cos(x3)-mphat*g*(sin(x3)^2)*(mphat*(sin(x3)^2)+mchat)+(cos(x3)^2)*(mphat*(sin(x3)^2)-mchat))...
        /((mphat*(sin(x3)^2)+mchat)^2);
    A(2,4)=2*mphat*lhat*sin(x3)*x4/(mchat+mphat*(sin(x3)^2));
    A(3,4) = 1;
    A(4,3)=((u1*sin(x3)*(mphat*(sin(x3)^2)+mchat)-g*cos(x3)*(mchat-mphat*(sin(x3)^2)))/(lhat*(mchat+mphat*(sin(x3)^2))))...
        +((-mphat*x4^2*(sin(x3)^2)*(mchat+mphat*(sin(x3)^2)))+(cos(x3)^2)*(mphat*(sin(x3)^2)-mchat))/((mchat+mphat*(sin(x3)^2))^2);
    A(4,4)=2*mphat*x4*cos(x3)*sin(x3)/(mchat+mphat*(sin(x3)^2));

    B(2,1)=coeff1;
    B(4,1)=-cos(x3).*coeff2;

end





