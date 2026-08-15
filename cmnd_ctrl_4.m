clc
clear all
close all
T = 20; %time length
nc = 1; %number of controls#######
M = T*24; %number of grid
ef = 1;
%%%%%%%%%%%%%%%%%%%%%%%%% without control %%%%%%%%%%%%%%%%%%%%%%%%%
global b r k delT xi dN qT dT qE a  dE dP g mu eta w c K 
r = 1*1.1e5; k = 5.74e-4;  xi = 0.52; dN = 0.07; qE = 0.01; 
a = 2e-5;  dE = 5;  g = 2200; K = 3.5e3; mu = g/K; 
dP = 0.02;   
%b = 1.54e-10; %1.85e-10;%
b=4e-10;
delT = 3e-5;%0;% 
qT  = 1.2e-4;%0;
dT = 0.1;%0;


c = [0.823624685090885	0.00173698844810918];

ix = [3.5e5 8.62e-18 0 1000 0 K];

eta = 40; w = 1;

[t,x] = odeRK_c(ix,T,M);
%%%%%%%%%%%%%%%%%%%%%%%%% with control %%%%%%%%%%%%%%%%%%%%%%%%%
A1 = 0; A2 = 1; B = 1000;

[tc xc u] = fbsm(ix,T,M,nc,A1,A2,B,ef);%make change according to the objective function#######



%% figure
subplot(2,2,1)
     plot(t,x(2,:),'r-',tc,xc(2,:),'b--','LineWidth',2)
     title('(a)'); ylabel('I'); legend('without control','with control'); hold on;
   
subplot(2,2,2)
    plot(t,x(6,:),'r-',tc,xc(6,:),'b-.','LineWidth',2)
    title('(b)'); ylabel('P'); hold on; %legend('without control','with control'); 
%     ylim([2800 3500])
%     yticks([2800 2900 3000 3100 3200 3300 3400 3500])
subplot(2,2,3)  
    plot(tc,u(:),'b-','LineWidth',2);
    hold on 
   title('(c)'); ylabel('u^*'); xlabel('Days'); hold on;
ll=length(u);  
subplot(2,2,4)
   %HourlyFluidRate = eta*(K-xc(6,:)).*u*T/M;
   DailyFluidRate = eta*(K-xc(6,:)).*u;
   plot(tc,cumsum(eta*(K-xc(6,:)).*u)*T/M,'b-',tc,DailyFluidRate,'b--','LineWidth',2);
   title('(d)'); ylabel('Fluid(ml)'); xlabel('Days'); legend('Cumulative Volume','Daily Infusion Rate');hold on;

