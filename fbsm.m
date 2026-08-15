function [t,x,u] = fbsm(ix,T,M,nc,A1,A2,B,ef)%make change according to the objective function#######
delta = 0.0001;%tollerence
t = linspace(0,T,M+1);dt = T/M; dt2 = dt/2; dt6 = dt/6;
n = max(size(ix));
cache = zeros(n,4);

x = zeros(n,M+1); ld = zeros(n,M+1); u = zeros(nc,M+1);
x(:,1) = ix;

test = -1;
while(test<0)
ou = u; ox = x; old = ld;
for i = 1:M
    cache(:,1) = state_c(x(:,i),u(:,i));%[m11 m12 m13 m14]'
    cache(:,2) = state_c(x(:,i)+dt2.*cache(:,1),(u(:,i)+u(:,i+1))./2);%[m21 m22 m23 m24]'
    cache(:,3) = state_c(x(:,i)+dt2.*cache(:,2),(u(:,i)+u(:,i+1))./2);%[m31 m32 m33 m34]'
    cache(:,4) = state_c(x(:,i)+dt.*cache(:,3),u(:,i+1));%[m41 m42 m43 m44]'
    x(:,i+1) = x(:,i) + dt6.*cache*[1 2 2 1]';
end

for i = M+1:-1:2
    cache(:,1) = costate(x(:,i),ld(:,i),u(:,i),A1,A2);%make change according to the objective function#######
    cache(:,2) = costate((x(:,i)+x(:,i-1))./2,ld(:,i)-dt2.*cache(:,1),(u(:,i)+u(:,i-1))./2,A1,A2);%make change according to the objective function#######
    cache(:,3) = costate((x(:,i)+x(:,i-1))./2,ld(:,i)-dt2.*cache(:,2),(u(:,i)+u(:,i-1))./2,A1,A2);%make change according to the objective function#######
    cache(:,4) = costate(x(:,i-1),ld(:,i)-dt.*cache(:,3),u(:,i-1),A1,A2);%make change according to the objective function#######
    ld(:,i-1) = ld(:,i) - dt6.*cache*[1 2 2 1]';
end

u = control(x,ld,nc,M,B);%make change according to the objective function#######
%Bddu = min(ef,max(0,u));
%u = (ou+Bddu)./2;
Bddu1 = min(ef,max(-10,u));%ef is efficiency of controls in fraction
u = (ou+Bddu1)./2;


t1 = min(delta.*sum(abs(u),2) - sum(abs(ou - u),2));
t2 = min(delta.*sum(abs(x),2) - sum(abs(ox - x),2));
t3 = min(delta.*sum(abs(ld),2) - sum(abs(old - ld),2));
test = min(min(t1,t2),t3);
end