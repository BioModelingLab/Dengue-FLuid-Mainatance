clc;
clear all;
close all;
T = 20; %time length
nc = 1; %number of controls#######
M = T*24; %number of grid
ef = 1;
%%%%%%%%%%%%%%%%%%%%%%%%% without control %%%%%%%%%%%%%%%%%%%%%%%%%
global b r k delT xi dN qT dT qE a  dE dP g mu eta w c K
r = 1*1.1e5; k = 5.74e-4;  xi = 0.52; dN = 0.07; qE = 0.01; 
a = 2e-5;  dE = 5;  g = 1200+1000; K = 3.5e3; mu = g/K; 
 w = 1; 
% b = 1.85e-10; 
b=4e-10;
delT = 3e-5; qT  = 1.2e-4; dT = 0.1;

c = [0.823624685090885	0.00173698844810918];

ix = [3.5e5 8.62e-18 0 1000 0 K];

A1 = 0; A2 = 1; 
eta_ = 5:2.5:60;
B_ = 1000:2500:50000;
totalFluid = zeros(numel(B_),numel(eta_)); maxRate = zeros(numel(B_),numel(eta_));
duration = zeros(numel(B_),numel(eta_)); maxDeficite = zeros(numel(B_),numel(eta_));
maxDeficite_c = zeros(numel(B_),numel(eta_));

%%
for i1 = 1:1:numel(B_)
    for i2 = 1:1:numel(eta_)
        B = B_(i1); eta = eta_(i2);
        [t,x] = odeRK_c(ix,T,M);
[tc xc u] = fbsm(ix,T,M,nc,A1,A2,B,ef);%make change according to the objective function#######
totalFluid(i1,i2) = sum(eta*(K-xc(6,:)).*u)*T/M;

maxRate(i1,i2) = max(eta*(K-xc(6,:)).*u)*M/T;
duration(i1,i2) = sum(u>0.05)*T/M;
maxDeficite(i1,i2) = max(K-x(6,:))*100/K;
maxDeficite_c(i1,i2) = max((K-xc(6,:)))*100/K;

    end
end

figure('Units','centimeters',...
       'Position',[2 2 18 6.5],...
       'Color','w');


ax1 = subplot(1,3,1);

surf(eta_,B_,maxDeficite,...
    'EdgeColor','none',...
    'FaceColor',[0.251 0.1843 0.01569],...
    'FaceAlpha',0.7);
hold on

surf(eta_,B_,maxDeficite_c,...
    'EdgeColor','none',...
    'FaceColor','interp');

[C, h]=contour(eta_,B_,maxDeficite_c,...
        [0 15],...
        'r',...
        'ShowText','on',...
        'LineWidth',2);
clabel(C,h,...
       'FontSize',9,...
       'Color','k');

[C,h]=contour3(eta_,B_,maxDeficite_c,...
               [5 5],...
               'k',...
               'ShowText','on',...
               'LineWidth',2);

clabel(C,h,...
       'FontSize',9,...
       'Color','k');

 t = title('(g) Maximum Fluid Deficit (%)',...
          'FontWeight','normal',...
          'FontSize',11);

set(t,'Units','normalized');
set(t,'Position',[0.5 1.07 0]);    % Default is around 1.02

lgd = legend('Without control','With control');

set(lgd,...
    'Units','normalized',...
    'Position',[0.12 0.80 0.16 0.08]);


xlabel('\eta','FontSize',10)
ylabel('B','FontSize',10)
zlabel('Fluid Deficit (%)','FontSize',10)

axis([eta_(1) eta_(end) B_(1) B_(end)])

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',9,...
    'LineWidth',1,...
    'TickDir','out',...
    'Box','off');

hold off


ax2 = subplot(1,3,2);

[C,h] = contour(eta_,B_,totalFluid,...
                10,...
                'LineWidth',2);

clabel(C,h,...
       'FontSize',9,...
       'Color','k',...
       'FontName','Times New Roman');

xlabel('Infusion coefficient (\eta)','FontSize',10)
ylabel('Overload constraint (B)','FontSize',10)

zlabel('Total Fluid','FontSize',10)

t2 = title('(h) Total Fluid (ml)',...
           'FontWeight','normal',...
           'FontSize',11);

set(t2,'Units','normalized');
set(t2,'Position',[0.5 1.06 0]);    % Increase 1.06 if needed

axis([eta_(1) eta_(end) B_(1) B_(end)])

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',9,...
    'LineWidth',1,...
    'TickDir','out',...
    'Box','on');


ax3 = subplot(1,3,3);

[C,h] =contour(eta_,B_,duration,10,...
        'ShowText','on',...
        'LineWidth',2);
clabel(C,h,...
       'FontSize',9,...
       'Color','k',...
       'FontName','Times New Roman');

xlabel('Infusion coefficient (\eta)','FontSize',10)
ylabel('Overload constraint (B)','FontSize',10)
% zlabel('Duration of Fluid Support','FontSize',10)

t3 = title('(i) Duration of Fluid Support (Days)',...
           'FontWeight','normal',...
           'FontSize',11);

set(t3,'Units','normalized');
set(t3,'Position',[0.5 1.06 0]);

axis([eta_(1) eta_(end) B_(1) 5e4])

axis tight;




