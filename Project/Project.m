load proj23.mat

rain = ElGeneina.rain_org

%% Should you transform the data?
figure
hold on
lambda_max = bcNormPlot( rain );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
hold off
xlabel('Lambda')
ylabel('Log-likelihood')
%% Transform data using a offset and log 
offset = 0.1
rain_t = log(rain + offset)

%% Define the polynomials and initiate the model. 
A = [1 0]
B=[]
C = [ ]
model_init = idpoly( A, B ,C) ;
%%
model = pem( rain_t , model_init)
error = resid(model, rain_t); error(length(A):end)
[autoseas,pautoseas] = acf_pacf_norm(error)
whitenessTest(error)
%%
data = rain_t;      % MED TRANSFORM

% Initialize a new array with space for additional points
rain_3 = nan(1, 3 * length(data) - 2);

% Copy existing data and insert empty points
rain_3(1:3:end) = data;
rain_3(2:3:end) = nan;
rain_3(3:3:end) = nan;

%% Kalmar
close all
y = rain_3;
N = length(rain_3);
A = eye(3);
Re = eye(3);
Rv = 0.1;
Rx_t1 = 10*eye(3);                                      % Initial covariance matrix, R_{1|0}^{x,x}

x_hat = zeros(3,N);  % [x(t) x(t-1) x(t-2)]
yhat = zeros(N,1);
ehat = zeros(N,1);

for t = 3:N-2
    C = [1 1 1];
    % Prediction step
    x_t1 = A*x_hat(:,t-1);
    
 if ~isnan(y(t))
    Ry = C*Rx_t1*C' + Rv;
    K = Rx_t1 * C' / Ry;
    yhat(t) = C*x_t1;
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
%%

n = 3; % Number of state variables
T = length(rain_t); % Total number of time steps

% Initialize storage for state estimates
x_est_all = zeros(n, T);
% Define the state transition matrix A and measurement matrix H
A = [model.A(2), 0, 0; 1, 1, 0; 0, 1, 1]; % Example structure, adjust based on your AR(1) model
H = [1, 1, 1]; % Based on y_t = x_t + x_{t-1} + x_{t-2}
large_value = 10
% Initialize state estimate and error covariance
x_est = [1; 1; 1];
P_est = eye(3) * large_value; % Adjust as necessary

% Define process noise covariance Q and measurement noise covariance R
Q = eye(3);
R = 1;

% Assuming you have a series of y measurements
Y = rain_t;

% Implement the Kalman Filter
for i = 1:length(Y)
    % Prediction
    x_pred = A * x_est;
    P_pred = A * P_est * A' + Q;

    % Update
    K = P_pred * H' / (H * P_pred * H' + R);
    x_est = x_pred + K * (Y(i) - H * x_pred);
    P_est = (eye(3) - K * H) * P_pred;

    % x_est contains the updated estimates of x1, x2, and x3
    x_est_all(:, i) = x_est;
end
plot(x_est_all(:))
%%
list1 = x_est_all(:)
summed_list = [];
for i = 1:3:length(list1)
    summed_list(end + 1) = sum(list1(i:min(i+2, end)));
end