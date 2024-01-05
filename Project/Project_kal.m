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
alpha =-0.6
A0 = [1 alpha];
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
    h = figure(1)
    subplot(211)
    plot(x_real)
    title("x_t (Simulation)")
 
    subplot(212)
    plot(datat)
    title("Sum of x_t, x_{t-1} and x_{t-2} (Simulation)")

   

end
%set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', [pos(3), pos(4)])
print(h,"AR_sim" + string(abs(alpha*100)),'-dpng','-r900')
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
R = 1;   % measurement noise variance
% state Transition Matrix
a1s = linspace(-1,1 ,1000);
resids=zeros(1,length(a1s));
y_sums = zeros(length(datat),1);
x_sums = zeros(length(datat),1); %remove if real data
residsx=zeros(1,length(a1s));%remove if real data
% Basically doing a loop for a lot of a1-coeff, to see which one gives 
%best result when comparing sums
for i = 1:length(a1s)
    a1  = a1s(i);
    A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];
    % A = [-a1, 0, 0; 0, -a1, 0; 0, 0, -a1];
    
    % measurement Matrix
    H = [1, 1,1]; %"this is C
    
    % initial State Estimate and Covariance
    x_est = zeros(3, length(y)); % State estimates for each time step
    P = eye(3);                 % Initial state covariance, 3x3 identity matrix
    
    % process Noise Covariance Matrix
    % Q_matrix = Q * eye(3);
    Q_matrix = Q * [1 0 0;
                    0 1 0;
                    0 0 1];
    
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
    x_st = x_est(:,2:3:end);
    x_st_vec = x_st(:); % This looks about right. 
    for j = 1:length(data)
        y_sums(j) = x_st_vec(3*j) + x_st_vec(3*j-1)+ x_st_vec(3*j-2);
    end
    
    %create total mean squared errors for sums:
    resids(i) = mean((y_sums-data).^2);

end



% Output
[minValue, index] = min(resids);

fprintf('The minimum value is for a1 %d and its index is %d.\n', a1s(index), index);
fprintf('The min MSE was %d', min(resids))
%Set the a1 coeff to the "correct" a1 and rerun the kalman
a1 = a1s(index);
A = [-a1, 0, 0; 1, 0, 0; 0, 1, 0];

P = eye(3);  
% Iterate over y
for k = 1:length(y)
  
    % prediction Step
    if k == 1
        x_pred = zeros(3,1); % Initial prediction
        
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
%% Plot the residuals

h=figure
plot(a1s, resids,"Displayname", "Residuals")
title("Residual for different a_1 - coefficients",'FontSize',14)
% xline(alpha, 'Color', 'red', "Displayname", "Real a_1","LineStyle","--", "LineWidth",1.2)
xline(a1, 'Color', 'blue', "Displayname", "Estimated a_1","LineStyle","--", "LineWidth",1.2)
xlabel("a_1 coefficient",'FontSize',14)
ylabel("Residual",'FontSize',14)
legend("show", 'FontSize',14,'Location','North')
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', [pos(3), pos(4)])
print(h,"a1_resid_rain",'-dpng','-r900')

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

MSE = mean((y_sums-data).^2)
%Create plots:
h = figure
subplot(211) 
plot(x_st_vec, "DisplayName", "Estimated x_t")
title("Estimated x_t")
xlim([1 length(x_st_vec)])
xlabel("Time")
legend("show","Fontsize", 10)
hold on
if realdat == 0
    % plot(x_real,"DisplayName", "Real x_t")
% legend("show","Fontsize", 10)
end
hold off
subplot(212)

plot(y_sums,"DisplayName", "Estimated sum")
xlim([1 length(y_sums)])
hold on
title("Estimated vs real sum of x_t, x_{t-1} and x_{t-2}")
% subplot(224)
plot(datat,"DisplayName", "Real sum")
legend("show", "Fontsize", 10)
hold off
xlabel("Time")
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', [pos(3), pos(4)])
print(h,"xt_xtsum_rain",'-dpng','-r900')
% set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
%% plot results (rain)
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

MSE = mean((y_sums-data).^2)
%Create plots:
h = figure
subplot(211) 
plot(x_st_vec, "DisplayName", "Estimated x_t")
xlabel("Time")
title("Estimated x_t")
subplot(212)
plot(y_sums,"DisplayName", "Estimated sum")
xlim([1 length(y_sums)])
hold on
title("Estimated vs real sum of x_t, x_{t-1} and x_{t-2}")
% subplot(224)
plot(datat,"DisplayName", "Real sum")
legend("show", "Fontsize", 14)
hold off
xlabel("Time")
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize', [pos(3), pos(4)])
print(h,"xt_xtsum_rain",'-dpng','-r900')
%% Compare the sums in same graph
if realdat == 1
    y_sum_compare = rain_int_sum/3
elseif realdat == 0
    y_sum_compare = datat
end
figure(1)

plot( y_sums )
hold on
plot(datat)
title("estimated and real summed rain")
plot(y_sum_compare)
hold off

legend("Estimated", "Real", "Interpolated")
rain=x_st_vec
%% MSE:
MSE = mean((y_sums-data).^2)
%% Compare interpolated and actual
figure(1)
%real
plot(datat)
hold on
plot(y_sum_compare)
%inter
%% Remove negative values in x_st
for i = 1:length(x_st_vec)
    if x_st_vec(i) < 0
        x_st_vec(i) = 0
    end
end