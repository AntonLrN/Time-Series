% For plotting on 2 monitors:
monitorPositions = get(0, 'MonitorPositions');

% Assuming your second monitor is the second row
secondMonitor = monitorPositions(2, :);

% Define the size of your figure [width, height]
figSize = [1500, 850];

% Calculate position for the new figure
% The figure is placed at the center of the second monitor
figX = secondMonitor(1) + (secondMonitor(3) - figSize(1)) / 2;
figY = secondMonitor(2) + (secondMonitor(4) - figSize(2)) / 2;

% Create a new figure and set its position
figure;
set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
%% Run this if you wanna simulate AR(1) data etc and try Kalman on simulated data
plott = true %Do i want plots in this cell?
realdat = 0 % this is done so that later in plots we know if we are using proj. data or sim. data
rng(0)                                         
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
realdat = 1
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


%%
close all
%define stuff:
y = y; %The data with nans
data = datat; %y without nans

Q = 2;  % process noise variance                        %%% (Q,R) = (2,1), (10, 10) looks pretty good. 
R = 20;   % measurement noise variance
% state Transition Matrix
a1s = linspace(-0.999, 0.999,1000);
resids=zeros(1,length(a1s));
y_sums = zeros(length(datat),1);

% Basically doing a loop for a lot of a1-coeff, to see which one gives 
%best result when comparing sums
for i = 1:length(a1s)
    a1  = a1s(i);
    A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
    
    % measurement Matrix
    H = [1, 1,1]; %"this is C
    
    % initial State Estimate and Covariance
    x_est = zeros(3, length(y)); % State estimates for each time step
    P = eye(3);                 % Initial state covariance, 3x3 identity matrix
    
    % process Noise Covariance Matrix
    Q_matrix = Q * eye(3);
    
    % Iterate over y
    %KALMAN FILTER - starts here
    for k = 1:length(y)
        % Prediction Step
        if k == 1 %I.e. if we 
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
    %Kalman filter ends here, rest is residual calc. 
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
fprintf('The minimum value is for a1 %d and its index is %d.\n', a1s(index), index);

%Set the a1 coeff to the "correct" a1 and rerun the kalman
a1 = a1s(index);
A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];


% Iterate over y
for k = 1:length(y)
  
    % prediction Step
    if k == 1
        x_pred = zeros(3,1); % Initial prediction
        %x_pred =[; % Initial prediction
    else
        x_pred = A * x_est(:, k-1);
    end
    P_pred = A * P * A' + Q_matrix;

    % update Step
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


%% plot results
x_st = x_est(:,3:3:end); % Take every third state (becuase we only have measurements on every third iteration
x_st_vec = x_st(:); % This looks about right. 
if realdat==1 %this creates sum for the interpolated rain, to compare with if we are using the project data.
    rain_int_sum = zeros(1, length(rain_int)/3)
end
% create sums
y_sums=zeros(length(datat),1)
for i = 1:length(datat)
    y_sums(i) = x_st_vec(3*i) + x_st_vec(3*i-1)+ x_st_vec(3*i-2);
    if realdat==1 %This only runs if we are using project data
         rain_int_sum(i) = rain_int(3*i) + rain_int(3*i-1)+ rain_int(3*i-2);
    end
end


%Create plots:
subplot(221) 
plot(x_st_vec)
title("Estimated daily rain")
hold on
%subplot(222)
if realdat == 0
    plot(x_real)
    title("Actual daily rain")
end
hold off
subplot(223)

plot(y_sums)
title("Estimated summed rain")
subplot(224)
plot(datat)
title("Actual summed rain")
set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
%% Compare the sums in same graph
if realdat == 1
    y_sum_compare = rain_int_sum/3
elseif realdat == 0
    y_sum_compare = datat
end
figure(1)
plot( y_sums )
hold on
plot(y_sum_compare)
title("estimated and real summed rain")
hold off

rain=x_st_vec
