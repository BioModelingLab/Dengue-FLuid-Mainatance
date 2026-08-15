global b r k delT xi dN qT dT qE a g mu K T M 

r = 1*1.1e5; k = 5.74e-4;  xi = 0.52; dN = 0.07; qE = 0.01; 
a = 2e-5; g = 1200+1000; K = 3.5e3; mu = g/K; 
b = 1.85e-10; delT = 3e-5; qT  = 1.2e-4; dT = 0.1;

expTime=5:1:14;
V = 5500; hctB = 42; hctV = 5500*hctB/100; 
% K = 3.5e3;
DATA = hctB+[0 13.5 9 20 24 25 9 5 10.5 0]*hctB/100;
pdata=V*(ones(length(DATA),1)-DATA'/100);

T = 20; %time length
M = T*24; %number of grid
times=0:1:T;
%% Fitting related %%%%%%%%% Change
x0fcn = @(params) [3.5e5; 8.62e-18; 0; 1000; 0; K];
yfcn = @(x,params) (x(:,6));
params = [0.01, 0.001];
p=2;
paramnames = {'d_C', 'd_P'};
%% Parameter Estimation ML
[paramestsML, fval] = fminsearch(@(p) seiqrCost(times,p,pdata,x0fcn,yfcn),params,optimset('Display','iter','MaxFunEvals',5000,'MaxIter',5000));
ODE_Sol=odeRK(@stateRK,x0fcn,T,M,paramestsML);
simY = ODE_Sol.deval(times); 
yest = yfcn(simY',paramestsML);

%% Extract Solutions
test=ODE_Sol.t;
xest=ODE_Sol.y(6,:);
set(gca,'LineWidth',2,'FontSize',10,'FontName','Times')
plot(expTime,pdata,'o-','Linewidth',2)
hold on
plot(test,xest,'-','Linewidth',2)
legend('Data','model ','Location','nw'); 
ylabel('Plasma level');  
xlabel('time(days)'); box on
axis tight;
print(gcf, 'Fitting_BILINEAR_Final', '-depsc', '-r300');

%% Fisher Information Matrix ML
FIM = MiniFisher(times,paramestsML,x0fcn,@stateRK_wrapper,yfcn);
rank(FIM);
%% Generate Profile Likelihoods ML
figure(2)
profiles = [];
% Wrapper function for parameter estimation
costfun = @(p) seiqrCost(times,p,pdata,x0fcn,yfcn);
threshold = chi2inv(0.95,length(paramestsML))/2 + fval;
profrange = 0.1; %0.1

YLABEL={'PL_{d_C}', 'PL_{d_P}'};
for i=1:length(paramestsML)
    profiles(:,:,i) = ProfLike(paramestsML,i,costfun,profrange);
    subplot(1,2,i)
    plot(profiles(:,1,i),profiles(:,2,i),'LineWidth',2);
    hold on;
    plot(paramestsML(i),fval,'r*','LineWidth',2);
    plot([profiles(1,1,i) profiles(end,1,i)],[threshold threshold],'r--')
    hold on;
    xlabel(paramnames{i})
    ylabel(YLABEL{i})
    box on;
end
