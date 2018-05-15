%% INITIALIZE
clear all
%close all

% load pressure and flow data files
load pqdata.mat

% Make the data periodic by adding last samole equal to the first
% N=length(pdata)+1;
% pdata(N)=pdata(1);
% qdata(N)=qdata(1);

% % string data together to create multiple heartbeats.
nrOfHb = 10; % nr of heartbeats to use
pdata_s = pdata;
qdata_s = qdata;
for hb = 1:nrOfHb-1
    pdata = [pdata(:);pdata_s(:)];
    qdata = [qdata(:);qdata_s(:)];
end
N = length(pdata);
N_s = length(pdata_s);

%% DEFINE PARAMETERS
% pre defined parameters
f=1;              % cardiacal frequency [Hz]
L = 0.008;
Z = 0.008;

% calculate parameters
T=1/f;            % cardiacal period [s]

% create a time axis vector - transpose to create a row vector instead of a
% column vector
t=linspace(0,nrOfHb*T,N);%0:dt:Tt;
t=t';

% calculate peripheral resistance
SV = mean(qdata_s)*T; %stroke volume [ml] = mean flow [ml/s] * total time of one stroke [s]
HR = f; %heart rate [1/s] = frequency [Hz]
CO = SV*HR; %cardiac output [ml/min] = stroke volume [ml] * heart rate [1/s] 
R = mean(pdata_s)/CO;% peripheral resistance [kPa*s/ml]
%R = 0.1781;

% calculate C from RC time (decay time method)
C=(.5/R)/log(pdata_s(round(N_s/2))/pdata_s(N_s));
% C = 12.2738;

%% 2 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R C Z L]';

% execute model function with inital parameters
[sse1, p1, ssei1] = pqmodel(x0, 1, f, pdata, qdata, nrOfHb);

% optimize parameters with fminsearch --- model 1
x = fminsearch(@(x)pqmodel(x, 1, f, pdata, qdata, nrOfHb),x0);

% execute model function with optimized parameters
[sse1o, p1o, ssei1o] = pqmodel(x, 1, f, pdata, qdata, nrOfHb);

%output data
tData1 = {'--- Model 1 ---' 0 0;
    'Peripheral resistance [kPa*s/ml]' x0(1) x(1);
	'Total arterial compliance [ml/kPa]' x0(2) x(2);
	'tau (RC) [s]' x0(1)*x0(2) x(1)*x(2);
	'SSE' sse1 sse1o};



%% 3 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R-Z C Z L]';

% execute model function with inital parameters
[sse2, p2, ssei2] = pqmodel(x0, 2, f, pdata, qdata, nrOfHb);

% optimize parameters with fminsearch --- model 2
x = fminsearch(@(x)pqmodel(x, 2, f, pdata, qdata, nrOfHb),x0);

% execute model function with optimized parameters
[sse2o, p2o, ssei2o] = pqmodel(x, 2, f, pdata, qdata, nrOfHb);

%output data

tData2 = {'--- Model 2 ---' 0 0;
	'Peripheral resistance [kPa*s/ml]' x0(1) x(1);
	'Total arterial compliance [ml/kPa]' x0(2) x(2);
	'tau (RC) [s]' x0(1)*x0(2) x(1)*x(2);
	'Aortic characteristic impedance [kPa*s/ml]' x0(3) x(3);
	'SSE' sse2 sse2o};


%% 3 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R C Z L]';

% execute model function with inital parameters
[sse3, p3, ssei3] = pqmodel(x0, 3, f, pdata, qdata, nrOfHb);

% optimize parameters with fminsearch --- model 3
x = fminsearch(@(x)pqmodel(x, 3, f, pdata, qdata, nrOfHb),x0);

% execute model function with optimized parameters
[sse3o, p3o, ssei3o] = pqmodel(x, 3, f, pdata, qdata, nrOfHb);

% output data to uitable
tData3 = {'--- Model 3 ---' 0 0;
	'Peripheral resistance [kPa*s/ml]' x0(1) x(1);
	'Total arterial compliance [ml/kPa]' x0(2) x(2);
	'tau (RC) [s]' x0(1)*x0(2) x(1)*x(2);
	'Aortic characteristic impedance [kPa*s/ml]' x0(3) x(3);
	'Total arterial inertance ' x0(4) x(4);
	'SSE' sse3 sse3o};

%% PLOT FIGURES

% create figure with resizing uipanel
f = figure('Visible','off');
% set figure to fullscreen
set(f,'Units','normalized','outerposition',[0 0 1 1])

% create table with numeric data
table = uitable;
set(table, 'Units', 'normalized', 'Position', [0.05 0.55 .45 .4], 'ColumnName', {'Value', 'Estimated', 'Optimized'}, 'ColumnFormat',  {'char', 'numeric', 'numeric'}, 'ColumnWidth',{250 'auto' 'auto'});
tData = vertcat(tData1, tData2, tData3);
set(table,'Data',tData)


% plot estimated pressure data with standard deviation and measured data
subplot(2,2,2)
hold on
plot(t,pdata, 'bo-')
plot(t,p1,'r')
plot(t,p2,'g')
plot(t,p3,'y')
bg = [ bar(t,ssei1,'FaceColor',[1 .5 .5]); bar(t,ssei2,'FaceColor',[.5 1 .5]); bar(t,ssei3,'FaceColor',[1 1 .5])];


gcl = legend('Pressure data','Model 1 Pressure', 'Model 2 Pressure','Model 3 Pressure','SD model 1','SD model 2','SD model 3');
legend('boxoff')
set(gcl,'TextColor',[1 1 1])
title('Estimated pressure with 3 Windkessel models')
xlabel('time [s]')
ylabel('pressure [kPa]')
axis([0 nrOfHb*T 0 20])
grid on
set(gca,'Color',[.4 .4 .4]);
set(gca,'FontSize',8)


% plot optimized pressure data with standard deviation and measured data
subplot(2,2,4)
hold on
plot(t,pdata, 'bo-')
plot(t,p1o,'r')
plot(t,p2o,'g')
plot(t,p3o,'y')
bg = [bg; bar(t,ssei1o,'FaceColor',[1 .5 .5]); bar(t,ssei2o,'FaceColor',[.5 1 .5]); bar(t,ssei3o,'FaceColor',[1 1 .5])];



gcl = legend('Pressure data','Model 1 Pressure', 'Model 2 Pressure','Model 3 Pressure','SD model 1','SD model 2','SD model 3');
legend('boxoff')
set(gcl,'TextColor',[1 1 1])
title('Optimised fit for pressure with 3 Windkessel models')
xlabel('heartbeat [n]')
ylabel('pressure [kPa]')
axis([0 nrOfHb*T 0 20])
grid on
set(gca,'Color',[.4 .4 .4]);
set(gca,'FontSize',8)
for j = 1:length(bg)
    set(bg(j),'edgecolor','none')
end


set(f,'Visible','on');
