function dx = stateRK_c(x)
global b r k delT xi dN qT dT qE a  g mu c

nsv = size(x);
dx = zeros(nsv(1),1);

dx(1) = -b*r*x(2)*x(1);   %S
dx(2) = b*r*x(2)*x(1) - k*x(2)*x(3) - delT*x(2)*x(4);              %I
dx(3) = xi*x(2) - dN*x(3);                %N
dx(4) = qT*x(2)*x(4) - dT*x(4);%T
dx(5) = qE*(1+a*x(4))*x(2) - c(1)*x(5);
dx(6) = g - c(2)*x(5)*x(6) - mu*x(6);
