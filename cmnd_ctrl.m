clc
clear all
close all
T = 20; %time length
nc = 1; %number of controls#######
M = T*24; %number of grid
ef = 1;
%%%%%%%%%%%%%%%%%%%%%%%%% without control %%%%%%%%%%%%%%%%%%%%%%%%%
global b r k delT xi dN qT dT qE a  dE dP g mu eta w c K h P0
r = 1*1.1e5; k = 5.74e-4;  xi = 0.52; dN = 0.07; qE = 0.01; 
a = 2e-5;  dE = 5;  g = 2200+500; h=5; 
% K = 3.5e3*(2200/g)^(1/(1+h)); mu = g/3.5e3; h=5; P0=3.5e3*(2200/g)^(1/(1+h));
K = 3.5e3; mu = 2200/3.5e3; P0=K;
dP = 0.02;   
%b = 1.54e-10; %1.85e-10;%
b=1.85e-10;
delT = 3e-5;%0;% 
qT  = 1.2e-4;%0;
dT = 0.1;%0;


c = [0.823624685090885	0.00173698844810918];

ix = [3.5e5 8.62e-18 0 1000 0 K];

eta = 40; w = 1;

[t,x] = odeRK(ix,T,M);
plot(t,x(6,:),'r-','LineWidth',2)
%%%%%%%%%%%%%%%%%%%%%%%%% with control %%%%%%%%%%%%%%%%%%%%%%%%%
A1 = 0; A2 = 1; B = 1000;

[tc xc u] = fbsm(ix,T,M,nc,A1,A2,B,ef);%make change according to the objective function#######
%[tc xc u] = Copy_of_fbsm(ix,T,M,nc,A1,A2,B,ef);%make change according to the objective function#######

% [peakVal, peakLoc] = findpeaks(u);
% totalFluid = sum((eta*(K-xc(6,:)).*u)*T/M)
% u_smooth = smooth(u, 5);
% du = diff(u);
% 
% % Define threshold for sharp rise
% threshold = max(u) * 0.003;  % Adjust as needed
% 
% % Find index where the sharp rise begins
% sharp_rise_idc=find(du > threshold, 1, 'first')
% 
% sharp_rise_time = tc(find(du > threshold, 1, 'first'))
% max_fluid_infusion_time=tc(peakLoc)
% 
% %Find point where it starts decreasing and reaches a minumum value 
% % start_decreasing_idx = find(du < 0, 1, 'first');
% start_decreasing_time = t(find(du < 0, 1, 'first'))
% 
% % Find the minimum value and its time
% [min_val, min_idx] = min(x(6,:));
% min_time = t(min_idx)
threshold = 0.2;

% First nonzero infusion
idx_start = find(u > 0.01, 1, 'first');

% Peak
[~, idx_max] = max(u);

% First point after the peak where u <= threshold
idx_cross = idx_max - 1 + find(u(idx_max:end) <= threshold, 1, 'first');

% Duration
duration= (idx_cross - idx_start + 1)*T/M;

% totalFluid_Loss = sum((K-x(6,:)))*T/M
% totalFluid = sum(eta*(K-xc(6,indices(1):indices(end))).*u(u>0.15))*T/M
% sum((u>0.15)*T/M)
%sum(u(idx_max:end) > threshold)*T/M

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
%     line([tc(sharp_rise_idx), tc(sharp_rise_idx)], [0 1], ...
%      'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5);
%  hold on
%  line([t(peakLoc(1)), t(peakLoc(1))], [0 1], ...
%      'Color', 'r', 'LineStyle', '-', 'LineWidth', 1.5);
% default_ticks = get(gca, 'XTick');
% custom_ticks = [t(sharp_rise_idx), t(peakLoc(1))];
% all_ticks = unique([default_ticks, custom_ticks]);
% set(gca, 'XTick', all_ticks);
%     hold on;
%     plot(tc(indices(1):indices(end)),u(indices(1):indices(end)),'r-','LineWidth',2)
%     hold on;
%     line([indices(1) indices(1)],[0 1.02],'r-','LineWidth',2)
%     hold on;
%     line([indices(end) indices(end)],[0 1.02],'r-','LineWidth',2)
%     hold on;
   %axis([0 T 0 max(ef)+0.02]);
   title('(c)'); ylabel('u^*'); xlabel('Days'); hold on;
   
subplot(2,2,4)
   %HourlyFluidRate = eta*(K-xc(6,:)).*u*T/M;
   DailyFluidRate = eta*(K-xc(6,:)).*u;
   plot(tc,cumsum(eta*(K-xc(6,:)).*u)*T/M,'b-',tc,DailyFluidRate,'b--','LineWidth',2);
   title('(d)'); ylabel('Fluid(ml)'); xlabel('Days'); legend('Cumulative Volume','Daily Infusion Rate');hold on;
%    ylim([0 5000])

set(gca,'FontSize',12)
print(gcf, 'optcontrol', '-depsc', '-r300')