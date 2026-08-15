clear all;
close all;
clc;
T = 20; %time length
%%%%%%%%%%%%%%%%%%%%%%%%% without control %%%%%%%%%%%%%%%%%%%%%%%%%
global b r k delT xi dN qT dT qE a  dE dP g mu eta w K c
r = 1*1.1e5; k = 5.74e-4;  xi = 0.52; dN = 0.07; qE = 0.01; 
a = 2e-5;  dE = 5;  g = 1200+1000; K = 3.5e3; mu = g/K; 
dP = 0.02; eta = 40; w = 1; 
b = 1.85e-10;%1.54e-10;
delT = 3e-5;%0;% 
qT  = 1.2e-4;%0;
dT = 0.1;%0;
ix = [3.5e5 8.62e-18 0 1000 0 K];
% T1 =  6.5:0.5:9.0;%8:-0.25:6.5;%[6.75 7 7.25 7.5 7.75 8];
T1=[6.5000  6.75  7.0000  7.25  7.5000    8.0000    8.5000    9.0000];
totalFluid = zeros(numel(T1),1);
initialDose = zeros(numel(T1),1);
c = [0.823624685090885	0.00173698844810918];

M = T*24; 
[t,x] = odeRK(ix,T,M);
subplot(1,2,1)
% set(gca,'fontsize',12);
plot(t,x(6,:),'r-','LineWidth',2); hold on;
%%%%%%%%%%%%%%%%%%%%%%%%% with control %%%%%%%%%%%%%%%%%%%%%%%%%
ef = 1;A1 = 0; A2 = 1; B = 1000;
nc = 1; %number of controls#######
linS = {'--k','-g','--g',':g','-b','--b',':b','-m','--m',':m','-c','--c',':c'};

for ii = 2:1:numel(T1)
    M = T1(ii)*24; %delay
    [t1,x1] = odeRK(ix,T1(ii),M);
    Tc = T-T1(ii);
   
    M = Tc*24; %number of grid
    
    ixc = x1(:,end);
    [tc xc u] = fbsm(ixc,Tc,M,nc,A1,A2,B,ef);%make change according to the objective function#######
    totalFluid(ii) = sum(eta*(K-xc(6,:)).*u)*Tc/M;
    sum(u>0.01)*Tc/M
    HourlyFluidRate = eta*(K-xc(6,:)).*u*Tc/M;
    initialDose(ii) = HourlyFluidRate(1);

    subplot(1,2,1)
        plot(tc+T1(ii),xc(6,:),linS{ii},'LineWidth',2)
        ylabel('P(ml)'); xlabel('Days'); hold on;
    subplot(1,2,2)
       DilyFluidRate = eta*(K-xc(6,:)).*u; HourlyFluidRate = eta*(K-xc(6,:)).*u*Tc/M;
       plot(24*(tc+T1(ii)),HourlyFluidRate,linS{ii},'LineWidth',2)
       %axis([0 T 0 max(ef)+0.02]);
       ylabel('Fluid(ml)'); xlabel('Hours'); hold on;
end
ii = 1
M = T1(ii)*24; %delay
    [t1,x1] = odeRK(ix,T1(ii),M);
    Tc = T-T1(ii);
   
    M = Tc*24; %number of grid
    
    ixc = x1(:,end);
    [tc xc u] = fbsm(ixc,Tc,M,nc,A1,A2,B,ef);%make change according to the objective function#######
    totalFluid(ii) = sum(eta*(K-xc(6,:)).*u)*Tc/M;
    sum(u>0.01)*Tc/M
    HourlyFluidRate = eta*(K-xc(6,:)).*u*Tc/M;
    initialDose(ii) = HourlyFluidRate(1);

    subplot(1,2,1)
        plot(tc+T1(ii),xc(6,:),linS{ii},'LineWidth',2)
        ylabel('P(ml)'); xlabel('Days'); hold on;
        xlim([6.0 9.5]);
        legend('Without Fluid Support','Started from 6.75th day','Started from 7.0th day','Started from 7.25th day','Started from 7.5th day',...
            'Started from 8.0th day', 'Started from 8.5th day', 'Started from 9.0th day',...
            'Started from 6.5th day','Location','southwest');
    
       subplot(1,2,2)
       DilyFluidRate = eta*(K-xc(6,:)).*u; HourlyFluidRate = eta*(K-xc(6,:)).*u*Tc/M;
       plot(24*(tc+T1(ii)),HourlyFluidRate,linS{ii},'LineWidth',2)
       %axis([0 T 0 max(ef)+0.02]);
       ylabel('Fluid(ml)'); xlabel('Hours'); hold on;
       xlim([160 260]);
%        set(gca,'fontsize',12);
       print(gcf, 'HourlyDelay', '-depsc', '-r300')
figure
   %\{subplot(1,2,1)
   %plot(T1,initialDose,'bo','LineWidth',2)
   %ylabel('Initial Dose(ml/hour)'); xlabel('Days'); hold on;
   %subplot(1,2,2)
   plot(T1,totalFluid,'bo','LineWidth',2)
   ylabel('Total Fluid(ml)'); xlabel('Days'); hold on;
   xlim([6.4 9.0]);
   print(gcf, 'TotalFluid', '-depsc', '-r300')

