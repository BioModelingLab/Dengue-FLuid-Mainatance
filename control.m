function u = control(x,ld,nc,M,B)
global b r k delT xi dN qT dT qE a  dE dP g w K eta
u = zeros(nc,M+1); 
u = -ld(6,:).*w.*eta.*(K-x(6,:))./B;
% if u<0.7125
%     u=u;
% else
%     u=0;
% end
% u = -ld(6,:).*w.*eta./B;