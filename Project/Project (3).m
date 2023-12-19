%% Project
clear; clc; close all
load proj23.mat

rain = ElGeneina.rain;
rain_org = ElGeneina.rain_org;

plot(rain,'r')
figure
plot(rain_org,'b')

%% kolla om transform behövs
figure; 
lambda_max = bcNormPlot(rain_org,1);        % log-transform är lämplig
rain_org_trans = log(rain_org+1);
figure
plot(rain_org_trans)

%% Utöka rain_org 3ggr med tomma punkter 
close all
%data = rain_org_trans;      % MED TRANSFORM
data = rain_org;            % UTAN TRANSFORM

% Initialize a new array with space for additional points
rain_3 = nan(1, 3 * length(data));

% Copy existing data and insert empty points
rain_3(1:3:end) = data;
rain_3(2:3:end) = nan;
rain_3(3:3:end) = nan;

%% Kalman
close all
y = rain_3;
N = length(rain_3);
a1 = 0.1446;        % 0.0139 with Rv=0.0001
A = [a1 1 0;
      0 a1 1;
      0 0 a1];
Re = 2*eye(3);
Rv = 3;
Rx_t1 = 10*eye(3);                                      % Initial covariance matrix, R_{1|0}^{x,x}

x_hat = zeros(3,N);  % [x(t) x(t-1) x(t-2)]
yhat = zeros(N,1);
ehat = zeros(N,1);

for t = 3:N
    C = [1 1 1];
    % Prediction step
    x_t1 = A*x_hat(:,t-1);
    Ry = C*Rx_t1*C' + Rv;
    K = Rx_t1 * C' / Ry;
    yhat(t) = C*x_t1;
 if ~isnan(y(t))
    % Ry = C*Rx_t1*C' + Rv;
    % K = Rx_t1 * C' / Ry;
    % yhat(t) = C*x_t1;
    ehat(t) = y(t) - yhat(t);
    x_hat(:,t) = x_t1 + K*ehat(t);
    Rx_t = Rx_t1 - K*Ry*K';
 else
     x_hat(:,t) = x_t1;
     Rx_t = Rx_t1;
 end
    % Update step
    Rx_t1 = A*Rx_t*A' + Re;
end
rain_3_est = x_hat(1,:);
%plot(yhat)
%% Plot results WITHOUT transform
for k=1:length(rain_3_est)
    if rain_3_est(k) < 0
        rain_3_est(k) = 0;
    end
end
figure;
hold on
plot(rain_3_est)
plot(rain./3)
hold off
legend('rain_3_est - kalman','rain - interpolated / 3')
sum(rain_3_est)-sum(rain)/3
%% sum last three data points to check if rain_3_est is the same as rain_org
close all
rain_sum = zeros(1,length(rain_org));
for i = 1:3:N
    rain_sum((i+2)/3) = sum(rain_3_est(i:i+2));
end
figure
hold on
plot(rain_sum,'r')          % looks good
plot(rain_org,'b')
hold off
legend('rain sum - kalman','rain org')

