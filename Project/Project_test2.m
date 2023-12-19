% % %% Project
% % clear; clc; close all
% load proj23.mat
% % 
% rain = ElGeneina.rain;
% rain_org = ElGeneina.rain_org;
% % 
% % plot(rain,'r')
% % figure
% % plot(rain_org,'b')
%%
% Get monitor positions
monitorPositions = get(0, 'MonitorPositions');

% Assuming your second monitor is the second row
secondMonitor = monitorPositions(1, :);

% Define the size of your figure [width, height]
figSize = [1200, 550];

% Calculate position for the new figure
% The figure is placed at the center of the second monitor
figX = secondMonitor(1) + (secondMonitor(3) - figSize(1)) / 2;
figY = secondMonitor(2) + (secondMonitor(4) - figSize(2)) / 2;

% Create a new figure and set its position
figure;
set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);


%% SIMULATE TEST DATA
rng(1)                                          % Set the seed (just done for the lecture!)
N  = 999;
a1 = 0.9
plott = true
x_real = zeros(N,1);
e = randn( N, 1 );

%Create AR(1) process
for k=2:N                                       % Implement filter by hand to allow a1 to change.
    x_real(k) = e(k) - a1*x_real(k-1);
end


y = zeros(N/3,1)
for i = 3:3:N
    y(i/3)=x_real(i) + x_real(i-1) + x_real(i-2)
end



% Initialize a new array with space for additional points
datat = y
y = nan(1, 3 * length(y));
y= y'


% Copy existing data and insert empty points
y(3:3:end) = datat;
y(2:3:end) = nan;
y(1:3:end) = nan
figure;
set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
if plott == true
    figure(1)
    subplot(211)
    plot(datat)
    title("Summed test")

    subplot(212)
    plot(x_real)
    title("AR(1) - process (test)")

end
%% FEV

% AR(1) coefficient
phi = 0.9;  % Example value, replace with your coefficient

% Noise variances
Q = 1;  % Process noise variance
R = 0.01;   % Measurement noise variance

% Number of time steps
N = length(y);  % Assuming y is your input vector

% Initialize state estimate and covariance
x_est = zeros(N, 1);  % Initial state estimate
P = 1;                % Initial state covariance

% Kalman filter implementation
for k = 2:N
    % Prediction Step
    x_pred = phi * x_est(k-1);
    P_pred = phi^2 * P + Q;

    % Update Step
    if isnan(y(k))
        % Skip update step if y is NaN
        x_est(k) = x_pred;
        P = P_pred;
    else
        % Compute Kalman gain
        K = P_pred / (P_pred + R);

        % Update state estimate and covariance
        x_est(k) = x_pred + K * (y(k) - x_pred);
        P = (1 - K) * P_pred;
    end
end

% Output
disp('Estimated x values:');
disp(x_est);
%% FEV 2
% AR(1) coefficient
phi = 0.9; % Example value, replace with your coefficient

% Noise variances
Q = 1; % Process noise variance for each x
R = 0.1;  % Measurement noise variance

% State transition matrix (for three consecutive states)
A = [phi, 1, 0; 0, phi, 1; 0, 0, phi];

% Measurement matrix (to sum three states)
H = [1, 1, 1];

% Initialize state estimate and covariance
x_est = zeros(3, length(y)/3); % State estimates for each set of 3 x values
P = eye(3); % Initial state covariance, 3x3 identity matrix

% Process noise covariance matrix
Q_matrix = Q * eye(3);

% Iterate over y
for k = 1:3:length(y)-2
    idx = ceil(k/3);

    % Prediction Step
    if idx == 1
        x_pred = zeros(3,1); % Initial prediction
    else
        x_pred = A * x_est(:, idx-1);
    end
    P_pred = A * P * A' + Q_matrix;

    % Update Step (only if y(k) is not NaN)
    if ~isnan(y(k))
        % Compute Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);

        % Update state estimate and covariance
        x_est(:, idx) = x_pred + K * (y(k) - H * x_pred);
        P = (eye(3) - K * H) * P_pred;
    else
        % If y(k) is NaN, use the prediction as the estimate
        x_est(:, idx) = x_pred;
        P = P_pred;
    end
end

% Output
disp('Estimated x values (three at a time):');
disp(x_est);

%% Figure(1)
figure(1)
plot(x_est)
hold on
figure(2)
plot(x)
title("this is real")
hold off

sums = zeros(N/3,1)
for i = 3:3:N
    sums(i/3) = x_est(i) + x_est(i - 1) + x_est(i - 2)
end
%%
figure(1)
plot(sums)
hold on
figure(2)
plot(datat)
title("this is real")
hold off

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
rain_3 = nan(1, 3 * length(data)-2);

% Copy existing data and insert empty points
rain_3(1:3:end) = data;
rain_3(2:3:end) = nan;
rain_3(3:3:end) = nan;

%% Kalman
close all
y = y% rain_3%y;
N = length(y);
y1 = y
a1 = 0.5;
A = [a1 0 0;
      0 a1 0 ;
      0  0 a1 ];
scale=1
Re = eye(3)*scale;
Rv = var(datat);
Rx_t1 = 10*eye(3);      %P                                % Initial covariance matrix, R_{1|0}^{x,x}

x_hat = zeros(3,N);  % [x(t) x(t-1) x(t-2)]
yhat = zeros(N,1);
ehat = zeros(N,1);

for t = 3:N %N-2 originally
    
    % Prediction step
    x_t1 = A*x_hat(:,t-1); %x_t1 is x_(t-1)
    C = [1 1 1 ];
    %NEW STUFF
    Ry = C*Rx_t1*C' + Rv;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    K = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;
    
 if isnan(y(t)) %If we have a number. 
    
    xt(:,t) = x_t1;                         % x_{t|t} = x_{t|t-1} 
    Rx_t    = Rx_t1;                        % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} 
    y1(t)   = yhat(t);
 else
   h_et(t) = y(t)-yhat(t);                 % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
   xt(:,t) = x_t1 + Kt*( h_et(t) );        % x_{t|t} = x_{t|t-1} + K_t ( y_t -  \hat{y}_{t|t-1} ) 
   Rx_t    = Rx_t1 - Kt*Ry*Kt';  
 end
    
    Rx_t1 = A*Rx_t*A' + Re;
end

%% FROM BOOK
p0 = 1;                                         % Number of unknowns in the A polynomial.
q0 = 0;                                         % Number of unknowns in the C polynomial.
y1 = y
A     = eye(p0+q0);
Rw    = 1;                                      % Measurement noise covariance matrix, R_w. Note that Rw has the same dimension as Ry.
Re    = 1e-6*eye(p0+q0);                        % System noise covariance matrix, R_e. Note that Re has the same dimension as Rx_t1.
Rx_t1 = eye(p0+q0);                             % Initial covariance matrix, R_{1|0}^{x,x}
h_et  = zeros(N,1);                             % Estimated one-step prediction error.
xt    = zeros(p0+q0,N);                         % Estimated states. Intial state, x_{1|0} = 0.
yhat  = zeros(N,1);                             % Estimated output.
xStd  = zeros(p0+q0,N);                         % Stores one std for the one-step prediction.
for t=4:N                                       % We use t-3, so start at t=4.
    % Update the predicted state and the time-varying state vector.
    x_t1 = A*xt(:,t-1);                         % x_{t|t-1} = A x_{t-1|t-1}
    C    = [ -y1(t-1)];    % Use earlier prediction errors as estimate of e_t.
    
    % Update the parameter estimates.
    Ry = C*Rx_t1*C' + Rw;                       % R_{t|t-1}^{y,y} = C R_{t|t-1}^{x,x} + Rw
    Kt = Rx_t1*C'/Ry;                           % K_t = R^{x,x}_{t|t-1} C^T inv( R_{t|t-1}^{y,y} )
    yhat(t) = C*x_t1;                           % One-step prediction, \hat{y}_{t|t-1}.

    % If a sample is missing, just retain the earlier state.
    if isnan( y(t) )
        xt(:,t) = x_t1;                         % x_{t|t} = x_{t|t-1} 
        Rx_t    = Rx_t1;                        % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} 
        y1(t)   = yhat(t);                      % Replace the missing sample with the estimated value. 
    else
        h_et(t) = y(t)-yhat(t);                 % One-step prediction error, \hat{e}_t = y_t - \hat{y}_{t|t-1}
        xt(:,t) = x_t1 + Kt*( h_et(t) );        % x_{t|t} = x_{t|t-1} + K_t ( y_t -  \hat{y}_{t|t-1} ) 
        Rx_t    = Rx_t1 - Kt*Ry*Kt';            % R^{x,x}_{t|t} = R^{x,x}_{t|t-1} - K_t R_{t|t-1}^{y,y} K_t^T
    end
    Rx_t1 = A*Rx_t*A' + Re;                     % R^{x,x}_{t+1|t} = A R^{x,x}_{t|t} A^T + Re

    % Estimate a one std confidence interval of the estimated parameters.
    xStd(:,t) = sqrt( diag(Rx_t) );
end
tt = 1:N;
figure
plot(tt, [y yhat] )
xlabel('Days')
if sum(isnan(y))
    hold on
    plot( tt(noVal), y0(noVal), 'b*')
    hold off
    legend('Realisation', 'Kalman estimate', 'Missing sample', 'Location','SW')
    title('One-step prediction using the Kalman filter with missing samples')
else 
    legend('Realisation', 'Kalman estimate', 'Location','SW')
    title('One-step prediction using the Kalman filter')
end




%% calculate errors.
y_pred = sum(x_hat())
y_resid = y_pred - datat

figure(1)
plot(y_resid)
title("residuals of sum")

%% sum last three data points to check if rain_3_est is the same as rain_org
close all
rain_sum = zeros(1,length(rain_org));
for i = 1:3:N-2
    rain_sum((i+2)/3) = sum(rain_3_est(i:i+2));
end
figure
hold on
plot(rain_sum,'r')          % looks good
plot(rain_org,'b')
hold off
legend('rain sum - kalman','rain org')
%% Plot rain_sum WITH INVERSE TRANSFORM
figure
hold on
plot(exp(rain_sum)-1,'r')          % looks good!
plot(rain_org,'b')
hold off
%% AR(1) model av rain_sum - använd ej??? EJ transform
y = rain_3_est';
ar_model = arx(y,[1]); 
rar = resid(ar_model,y);      rar = rar(length(ar_model.A):end );
present(ar_model)

plot(rar.y)
plotter(rar.y, 360)
%% simulera AR(1) med funna a1
close all
N = length(y);
y1 = zeros(N,1);
e = randn(N,1);
for k=2:N                                       % Implement filter by hand to allow a1 to change.
    y1(k) = e(k) - ar_model.a(2)*y1(k-1);
end

figure
plot(y1)
hold on

plot(y)