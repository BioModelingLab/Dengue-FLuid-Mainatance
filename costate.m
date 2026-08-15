function dld = costate(x,ld,u,A1,A2)
global b r k delT xi dN qT dT qE a  dE dP g mu eta w c K h P0
nsv = size(x);
dld = zeros(nsv(1),1);
dld(1) = b*r*x(2)*(ld(1)-ld(2));    %lambda1
dld(2) = b*r*x(1)*(ld(1)-ld(2)) + ld(2)*(k*x(3)+delT*x(4)) - ld(3)*xi - ld(4)*qT*x(4) - ld(5)*qE*(1+a*x(4));  %lambda3
dld(3) = ld(2)*k*x(2)+ld(3)*dN;   %lambda3
dld(4) = ld(2)*delT*x(2) - ld(4)*(qT*x(2)-dT)-ld(5)*qE*a*x(2);   %lambda4
dld(5) = ld(5)*c(1)+c(2)*ld(6)*x(6);
% dld(6) = A2 + ld(6)*(eta*u*w+mu)+c(2)*ld(6)*x(5);   %lambda6
dld(6) = A2 + ld(6)*(eta*u+c(2)*x(5)+mu*(h+1)*(x(6)/P0)^h);   %lambda6
% dld(5) = ld(5)*dE+ld(6)*c(1)*x(6)*exp((c(2)-x(5))/c(3))/(c(3)*(1+exp((c(2)-x(5))/c(3)))^2);               %typeIV
% dld(6) = A2 + ld(6)*(eta*w*u+c(1)/(1+exp((c(2)-x(5))/c(3)))+mu);   %typeIV
% dld(6) = A2 + ld(6)*(c(1)*x(5)/(1+c(2)*x(5))+mu);   %lambda6