%%%% Project %%%%
addpath('C:\Users\anton\Desktop\University\Time series analysis\Computer Exercises\Exercise 3')
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


%% Run this if you wanna simulate AR(1) data etc and try Kalman on simulated data
plott = true


rng(0)                                          % Set the seed (just done for the lecture!)
extraN = 100;
N  = 999;
A0 = [1 -0.95];
C0 = [1];                            
e  = randn( N+extraN, 1 );
x_real  = filter( C0, A0, e );   x_real = x_real(extraN+1:end);
tt = 1:N;
y = zeros(N/3,1)
%Create measurements. 
for i = 3:3:N
    y(i/3)=x_real(i) + x_real(i-1) + x_real(i-2);
end



% Initialize a new array with space for additional points
datat = y
y = nan(1, 3 * length(y));


% Copy existing data and insert empty points
y(3:3:end) = datat;
y(2:3:end) = nan;
y(1:3:end) = nan
% figure;
% set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
if plott == true
    figure(1)
    subplot(211)
    plot(datat)
    title("Summed test")

    subplot(212)
    plot(x_real)
    title("AR(1) i.e. not summed - process (test)")

end
%% Let us try on real data:
load proj23.mat
% We need the rain (summed over 3 data points, i.e. without nan:
rain_3 = ElGeneina.rain_org
rain_int = ElGeneina.rain
%include the nan
datat = rain_3 %The data
y = nan(1, 3 * length(rain_3));


%Copy existing data and insert empty points
y(3:3:end) = datat;
y(2:3:end) = nan;
y(1:3:end) = nan;

% y(1:3:end) = datat;
% y(2:3:end) = nan;
% y(3:3:end) = nan
%% kolla om transform behövs
% figure; 
% lambda_max = bcNormPlot(rain_org,1);        % log-transform är lämplig
% rain_org_trans = log(rain_org+1);
% figure
% plot(rain_org_trans)







%% sum last three data points to check if rain_3_est is the same as rain_org
close all
%define stuff:
y = y; %The one with nans
data = datat; %y without nans
% Known parameters

Q = 1;  % Process noise variance
R = 1;   % Measurement noise variance
% State Transition Matrix
a1s = linspace(-0.999, 0.999,1000);
resids=zeros(1,length(a1s));
y_sums = zeros(length(datat),1) ;
for i = 1:length(a1s)
    a1  = a1s(i);
    A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
    A = [-a1, 0, 0; 0, -a1, 0; 0, 0, -a1];
    
    % Measurement Matrix
    H = [1, 1,1];
    
    % Initial State Estimate and Covariance
    x_est = zeros(3, length(y)); % State estimates for each time step
    P = eye(3);                 % Initial state covariance, 3x3 identity matrix
    
    % Process Noise Covariance Matrix
    Q_matrix = Q * eye(3);
    
    % Iterate over y
    for k = 1:length(y)
        % Prediction Step
        if k == 1
            x_pred = zeros(3,1); % Initial prediction
        else
            x_pred = A * x_est(:, k-1);
        end
        P_pred = A * P * A' + Q_matrix;
    
        % Update Step
        if ~isnan(y(k))
            % Compute Kalman Gain
            K = P_pred * H' / (H * P_pred * H' + R);
    
            % Update state estimate and covariance
            x_est(:, k) = x_pred + K * (y(k) - H * x_pred);
            P = (eye(3) - K * H) * P_pred;
        else
            % If y(k) is NaN, use the prediction as the estimate
            x_est(:, k) = x_pred;
            P = P_pred;
        end
    end
    x_st = x_est(:,1:3:end);
    x_st_vec = x_st(:); % This looks about right. 
    for j = 1:length(data)
        y_sums(j) = x_st_vec(3*j) + x_st_vec(3*j-1)+ x_st_vec(3*j-2);
    end
    
    %create total mean squared errors for sums:
    resids(i) = mean((y_sums-data).^2);
    %resids_x(i) = mean((x_st_vec-x_real).^2);
end
% Output
[minValue, index] = min(resids);
title(fprintf('The minimum value is for a1 %d and its index is %d.\n', a1s(index), index));



a1 = a1s(index);

%Check errors for Re also, is this good idea or should this be based on
%intuition?


scales=linspace(0.01,100, 1000);
resids_s=zeros(1,length(scales));
for i = 1:length(scales)
    % a1  = a1s(index);
    Q = scales(i);
    A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
    A = [-a1, 0, 0; 0, -a1, 0; 0, 0, -a1];
    
    % Measurement Matrix
    H = [1, 1,1];
    
    % Initial State Estimate and Covariance
    x_est = zeros(3, length(y)); % State estimates for each time step
    P = eye(3);                 % Initial state covariance, 3x3 identity matrix
    
    % Process Noise Covariance Matrix
    Q_matrix = Q * eye(3);
    
    % Iterate over y
    for k = 1:length(y)
        % Prediction Step
        if k == 1
            x_pred = zeros(3,1); % Initial prediction
        else
            x_pred = A * x_est(:, k-1);
        end
        P_pred = A * P * A' + Q_matrix;
    
        % Update Step
        if ~isnan(y(k))
            % Compute Kalman Gain
            K = P_pred * H' / (H * P_pred * H' + R);
    
            % Update state estimate and covariance
            x_est(:, k) = x_pred + K * (y(k) - H * x_pred);
            P = (eye(3) - K * H) * P_pred;
        else
            % If y(k) is NaN, use the prediction as the estimate
            x_est(:, k) = x_pred;
            P = P_pred;
        end
    end
    x_st = x_est(:,1:3:end);
    x_st_vec = x_st(:); % This looks about right. 
    for j = 1:length(data)
        y_sums(j) = x_st_vec(3*j) + x_st_vec(3*j-1)+ x_st_vec(3*j-2);
    end
    
    %create total mean squared errors for sums:
    resids_s(i) = mean((y_sums-data).^2);
    %resids_x(i) = mean((x_st_vec-x_real).^2);
end




[minValue, index] = min(resids_s);
title(fprintf('The minimum value is for scale %d and its index is %d.\n', scales(index), index));
% CONTINUE WITH THE "GOOD" A1 and scale.
scale=scales(index);

%Try R
Rs = linspace(0.01, 50, 1000);
resids_r = zeros(1,length(Rs));
for i = 1:length(Rs)
    a1  = a1s(index);
    Q = scale;
    R = Rs(i);
    %test a thing
    % R=0.1 Using this i get negative rain, but the sum does look a little
    % bit nicer, maybe that's ok?
    A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
    A = [-a1, 0, 0; 0, -a1, 0; 0, 0, -a1];
    % Measurement Matrix
    H = [1, 1,1];
    
    % Initial State Estimate and Covariance
    x_est = zeros(3, length(y)); % State estimates for each time step
    P = eye(3);                 % Initial state covariance, 3x3 identity matrix
    
    % Process Noise Covariance Matrix
    Q_matrix = Q * eye(3);
    
    % Iterate over y
    for k = 1:length(y)
        % Prediction Step
        if k == 1
            x_pred = zeros(3,1); % Initial prediction
        else
            x_pred = A * x_est(:, k-1);
        end
        P_pred = A * P * A' + Q_matrix;
    
        % Update Step
        if ~isnan(y(k))
            % Compute Kalman Gain
            K = P_pred * H' / (H * P_pred * H' + R);
    
            % Update state estimate and covariance
            x_est(:, k) = x_pred + K * (y(k) - H * x_pred);
            P = (eye(3) - K * H) * P_pred;
        else
            % If y(k) is NaN, use the prediction as the estimate
            x_est(:, k) = x_pred;
            P = P_pred;
        end
    end
    x_st = x_est(:,1:3:end);
    x_st_vec = x_st(:); % This looks about right. 
    for j = 1:length(data)
        y_sums(j) = x_st_vec(3*j) + x_st_vec(3*j-1)+ x_st_vec(3*j-2);
    end
    
    x_st = x_est(:,1:3:end);
    x_st_vec = x_st(:); % This looks about right. 
    for j = 1:length(data)
        y_sums(j) = x_st_vec(3*j) + x_st_vec(3*j-1)+ x_st_vec(3*j-2);
    end
    %create total mean squared errors for sums:
    resids_r(i) = mean((y_sums-data).^2);
    %resids_x(i) = mean((x_st_vec-x_real).^2);
end





[minValue, index] = min(resids_r);
title(fprintf('The minimum value is for R %d and its index is %d.\n', Rs(index), index));
% CONTINUE WITH THE "GOOD" A1.
R=Rs(index);











% Output

A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
A = [-a1, 0, 0; 0, -a1, 0; 0, 0, -a1];
Q=scale;
% Measurement Matrix
H = [1, 1,1];

% Initial State Estimate and Covariance
x_est = zeros(3, length(y)); % State estimates for each time step
P = eye(3);                 % Initial state covariance, 3x3 identity matrix

% Process Noise Covariance Matrix
Q_matrix = Q * eye(3);

% Iterate over y
for k = 1:length(y)
  
    % Prediction Step
    if k == 1
        x_pred = zeros(3,1); % Initial prediction
    else
        x_pred = A * x_est(:, k-1);
    end
    P_pred = A * P * A' + Q_matrix;

    % Update Step
    if ~isnan(y(k))
        % Compute Kalman Gain
        K = P_pred * H' / (H * P_pred * H' + R);

        % Update state estimate and covariance
        x_est(:, k) = x_pred + K * (y(k) - H * x_pred);
        P = (eye(3) - K * H) * P_pred;
    else
        % If y(k) is NaN, use the prediction as the estimate
        x_est(:, k) = x_pred;
        P = P_pred;
    end
end

%% This will plot the rains etc. 
x_st = x_est(:,3:3:end);
x_st_vec = x_st(:); % This looks about right. 
rain_int_sum = zeros(1, length(rain_int)/3)
% create sums
y_sums=zeros(length(datat),1)
for i = 1:length(datat)
    y_sums(i) = x_st_vec(3*i) + x_st_vec(3*i-1)+ x_st_vec(3*i-2);
    rain_int_sum(i) = rain_int(3*i) + rain_int(3*i-1)+ rain_int(3*i-2);
end
%Plot the stuff:

%Create plots:
subplot(221) 
plot(x_st_vec)
title("Estimated daily rain")
subplot(222)
%plot(x_real)
%title("Actual daily rain")
subplot(223)
plot(y_sums)
title("Estimated summed rain")
subplot(224)
plot(datat)
title("Actual summed rain")

%% In same figure

figure(1)
plot( y_sums )
hold on
plot(rain_int_sum/3)
title("estimated and real summed rain")
hold off

rain=x_st_vec

%% Look only at rain, does it need transform etc?
% plot(x_st_vec)
figure; 
lambda_max = bcNormPlot(rain,1);        % log-transform är lämplig
rain_log= log(rain+1);
figure
plot(rain_log)
%% Let us start to try to build model
% Do we want to add the periodicity?
S = 27
AS = [1 zeros(1,S-1) -1]
%%
% Try 1
%acf_pacf_norm(rain_log) % Not super Normal dist, but close enough. Let us try with an AR(4) first.
A = [1 0 ]
A = conv(A, AS)
modelinit = idpoly(A, [], [1]) 
modelinit.Structure.a.Free=[1 1 zeros(1, 25) 1 1 ]
model =pem( rain_log, modelinit);

error = resid(rain_log, model ); error = error(length(A):end);
acf_pacf_norm(error)
present(model)
%% Try 2 We now try to add an MA part (3)
C = [1 0 0 1]
modelinit = idpoly(A,[],C)
modelinit.Structure.c.Free=[zeros(1,3) 1]
modelinit.Structure.a.Free=[1 1 zeros(1, 26) 1 ]
model = pem(rain_log, modelinit)
error = resid(rain_log, model ); error = error(length(A):end);
acf_pacf_norm(error) % seems to be 28 trend in data, I'll add a large A-term
present(model)
%% Try 3 We now try to add an additional AR part, maybe 6?
%A = [1 0 0 0 0 0 0]
C = [1 0 0 1 zeros(1, 23) 1]
modelinit = idpoly(A,[],C)
modelinit.Structure.c.Free=[zeros(1,3) 1 zeros(1, 23) 1]
modelinit.Structure.a.Free=[1 1 zeros(1, 26) 1 ]
model = pem(rain_log, modelinit);
present(model)
error = resid(rain_log, model ); error = error(length(A):end);
acf_pacf_norm(error);
% Dont know if adding 6 helped, it is significant but there is still a very
% clear periodicity

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