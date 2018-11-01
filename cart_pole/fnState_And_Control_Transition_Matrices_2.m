function [A,B] = fnState_And_Control_Transition_Matrices_2(x,u)

mp = .01;
mc = 1;
l = 0.25;
g=9.8;

    x1 = x(1,1);
    x2 = x(2,1);
    x3 = x(3,1);
    x4 = x(4,1);

    u1 = u(1,1);

    A = zeros(4,4);
    B = zeros(4,1);

    coeff1 = 1/(mc + mp.*sin(x3).^2);
    coeff2 = coeff1./l;

    A(1,2) = 1;
    A(2,3)=(x4^2*mp*cos(x3)*(mc-mp*(sin(x3)^2))-2*u1*mp*sin(x3)*cos(x3)-mp*g*(sin(x3)^2)*(mp*(sin(x3)^2)+mc)+(cos(x3)^2)*(mp*(sin(x3)^2)-mc))...
        /((mp*(sin(x3)^2)+mc)^2);
    A(2,4)=2*mp*l*sin(x3)*x4/(mc+mp*(sin(x3)^2));
    A(3,4) = 1;
    A(4,3)=((u1*sin(x3)*(mp*(sin(x3)^2)+mc)-g*cos(x3)*(mc-mp*(sin(x3)^2)))/(l*(mc+mp*(sin(x3)^2))))...
        +((-mp*x4^2*(sin(x3)^2)*(mc+mp*(sin(x3)^2)))+(cos(x3)^2)*(mp*(sin(x3)^2)-mc))/((mc+mp*(sin(x3)^2))^2);
    A(4,4)=2*mp*x4*cos(x3)*sin(x3)/(mc+mp*(sin(x3)^2));

    B(2,1)=coeff1;
    B(4,1)=-cos(x3).*coeff2;

end





