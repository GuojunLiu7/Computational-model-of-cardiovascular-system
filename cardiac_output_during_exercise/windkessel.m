%% INITIALIZE
clear all
close all


%% Load data
% load pressure and flow data files
load pqdata.mat

% Make the data periodic and an even amount of datapoints by adding last samole equal to the first
N=length(pdata)+1;
pdata(N)=pdata(1);
qdata(N)=qdata(1);

%create simplified flow data from actual flow data
t_points = [7 12 27]; % manually picked points from data
qdata_s = create_flow_data(t_points, qdata, 1);
qdata_s_double_sv = create_flow_data(t_points, qdata, 2);
t_points = [3 10 18]; % derived from literature and data
qdata_s_double_sv_and_hr = create_flow_data(t_points, qdata, 3);

pdata_s = pdata;
qdata = qdata_s;

% string data together to create multiple heartbeats.
nrOfHb = 50; % nr of heartbeats to use - has to be a even number
for hb = 1:nrOfHb-1
    pdata = [pdata(:);pdata_s(:)];
    if hb==20  % at the 20th heartbeat, increase stroke volume
        qdata_s = qdata_s_double_sv;
    end
    if hb==40 % at the 40th hearbeat, double bpm
        qdata_s = [qdata_s_double_sv_and_hr(:);qdata_s_double_sv_and_hr(:)];
    end
    qdata = [qdata(:);qdata_s(:)];
end
N = length(pdata);
N_s = length(pdata_s);

%% Define optimized parameters
R1=0.1785;
C1=11.4846;

R2=0.1695;
C2=12.8231;
Z2=0.0086;

R3=0.1777;
C3=12.4197;
Z3=0.0086;
L3=0.0243;


%% DEFINE  TIME PARAMETERS
% pre defined parameters
f=1;              % cardiacal frequency [Hz]

% calculate parameters
T=1/f;            % cardiacal period [s]

% create a time axis vector - transpose to create a row vector instead of a
% column vector
t=linspace(0,nrOfHb*T,N);%0:dt:Tt;
t=t';


%% 2 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R1 C1 0 0]';

% execute model function with inital parameters
p1 = pqmodel(x0, 1, f, pdata, qdata, nrOfHb);


%% 3 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R2-Z2 C2 Z2 0]';

% execute model function with inital parameters
p2 = pqmodel(x0, 2, f, pdata, qdata, nrOfHb);



%% 3 COMPONENT WINDKESSEL MODEL 
% set parameter vector
x0 = [R3 C3 Z3 L3]';

% execute model function with inital parameters
p3 = pqmodel(x0, 3, f, pdata, qdata, nrOfHb);


%% Extract beat data in different circumstances (5th, 15th, 25th, 35th, 45th)
mean_pressures = zeros(5,1);
pulse_pressures = zeros(5,1);
for bt = 1:5
    start_n = N_s*(bt*10 - 5);
    end_n = start_n+N_s;
    p_e = p2(start_n:end_n);
    mean_pressures(bt) = mean(p_e);
    pulse_pressures(bt) = max(p_e) - min(p_e);
end

%% PLOT FIGURES

% create figure with resizing uipanel
f = figure('Visible','off');
% set figure to fullscreen
set(f,'Units','normalized','outerposition',[0 0 1 1])

% plot estimated pressure data with standard deviation and measured data
subplot(2,1,1)
hold on
plot(t,qdata, 'b-')
subplot(2,1,2)
hold on
plot(t,p1,'r')
plot(t,p2,'g')
plot(t,p3,'y')

gcl = legend('Model 1 Pressure', 'Model 2 Pressure','Model 3 Pressure');
legend('boxoff')
set(gcl,'TextColor',[1 1 1])
xlabel('heartbeat [n]')
ylabel('pressure [kPa]')
axis([0 nrOfHb*T 0 20])
grid on
set(gca,'Color',[.4 .4 .4]);
set(gca,'FontSize',8)

set(f,'Visible','on');
