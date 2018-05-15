function [qdata] = create_flow_data(t_points, qdata, mode)
% t_points = [t1 t2 t3]; % important datapoints
    % t1 = start of systole
    % t2 = peak of systole flow
    % t3 = minimum of systole flow (negative flow due to closing valve)
% qdata = actual flow
% mode = 1, 2, 3:
    % 1: normal mode, 
    % 2: doubled amount of maximum flow, simulating increased strokevolume
    % 3: doubled flow and doubled heartrate, further increasing cardiac
    % output

t1 = t_points(1);
t2 = t_points(2);
t3 = t_points(3);

n = length(qdata); % n = amount of datapoints
max_q = max(qdata); % max = maximum flow (at t2)
min_q = min(qdata); % min = minimum flow (at t3)
std_q = mean(qdata(t3+1:end)); % std = standard flow during diastole (before t1 and after t3)

if mode > 1
    max_q = max_q*2;
end
if mode == 3
    n = round(n/2);
    std_q = mean(qdata(t3*2+1:end));
end
    

qdata = zeros(n,1);
for t = 1:n
    qdata(t) = std_q;
    
    if t>t1 && t<=t2
        delta_t = t2-t1;
        delta_q = max_q-std_q;
        qdata(t) = (t-t1)*delta_q/delta_t + std_q;
    end
    if t>t2 && t<=t3
        delta_t = t3-t2;
        delta_q = min_q-max_q;
        qdata(t) = (t-t2)*delta_q/delta_t + max_q;
    end    
end