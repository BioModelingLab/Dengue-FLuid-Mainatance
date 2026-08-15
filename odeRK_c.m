function [t,x] = odeRK_c(ix,T,M)%make change according to the objective function#######
delta = 0.0001;%tollerence
t = linspace(0,T,M+1);dt = T/M; dt2 = dt/2; dt6 = dt/6;
n = max(size(ix));
cache = zeros(n,4);
global b r k delT xi dN qT dT qE a  dE dP g mu eta w K
x = zeros(n,M+1);
x(:,1) = ix;

test = -1;
while(test<0)
ox = x;
for i = 1:M
    cache(:,1) = stateRK_c(x(:,i));%[m11 m12 m13 m14]'
    cache(:,2) = stateRK_c(x(:,i)+dt2.*cache(:,1));%[m21 m22 m23 m24]'
    cache(:,3) = stateRK_c(x(:,i)+dt2.*cache(:,2));%[m31 m32 m33 m34]'
    cache(:,4) = stateRK_c(x(:,i)+dt.*cache(:,3));%[m41 m42 m43 m44]'
    x(:,i+1) = x(:,i) + dt6.*cache*[1 2 2 1]';
end
test = min(delta.*sum(abs(x),2) - sum(abs(ox - x),2));
end