%% Use reconstructed rain as input model to predict NVDI
% clear; 
%run Project;
% close all; clc;
load proj23.mat
% load rain_3_est.mat
noLags = 333;
%% Transform to -1 to 1
nvdi = ElGeneina.nvdi;

min_nvdi = min(nvdi);
max_nvdi = max(nvdi);

normalized_nvdi = (nvdi - min_nvdi) / (max_nvdi - min_nvdi);


nvdi_s= 2 * normalized_nvdi - 1;



%%

%% Should we transform the data?
x = x_st_vec(22*36+1:end)';
%This to train on 
cutoff= 300 %this I just did because the first few hundred values in the rain seemed to not help.
x_input = x_st_vec(cutoff:end)';
% x and y need to be in phase
y = nvdi_s;           

figure
lambda_max = bcNormPlot( x_input );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
figure
lambda_max = bcNormPlot( y );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
xlabel('Lambda')
ylabel('Log-likelihood')
%% Transform and make stationary. They should now be stationary (maybe y has a trend)

x_offset = abs(min(x_input))+1
% % Both seems to indicate that a log-transform could be helpful.
x_input = log(x_input+x_offset);
% % x=log(x+abs(min(x))+1);
y=log(y+2);
% take away the means, which i suppose we have to do?
x_input = x_input - mean(x_input)
meany = mean(y)
y = y - meany
% y might have slight trend. 
%% extract modeling and validation sets and test
modelLim = 36*11%36*13;        % We select 13 years. Does the data seem stable enough?
% sets for y
modely = y(1:modelLim);  %modely = modely-mean(modely); 
valy = y(1:36*14); 
testy = y(1:36*16)
% If yd, change from 1 to 2
modelx = x(1:modelLim);  modelx = modelx-mean(modelx);
valx = x(1:36*14);  valx = valx-mean(valx);
testx = x(1:36*16); testx = testx-mean(testx);
% %% if do diff
% modelLim = 36*11%36*13;        % We select 13 years. Does the data seem stable enough?
% x=x(2:end)
% x = x_input(22*36+1 - (cutoff-1):end)';
% % sets for y
% modely = yd(1:modelLim);  %modely = modely-mean(modely); 
% valy = yd(1:36*14); 
% testy = yd(1:36*16)
% % If yd, change from 1 to 2
% modelx = x(1:modelLim);  %modelx = modelx-mean(modelx);
% valx = x(1:36*14);
% testx = x(1:36*16);
%% Create a model for the input.
plotACFnPACF(x_input, noLags, 'Input, x_t' ); %From this we can see a very clear season of 36
set(gcf, 'Position', [figX, figY, figSize(1), figSize(2)]);
%% There seems to be a strong periodicity at 36, suggesting that a differentiation might help.
%estimateARMA( modelx, [ 1 zeros(1,35) 1], [ 1 ], 'Differentiated input, version 2', noLags );



%% Lets try add a36 first this time
As = [1 zeros(1,35) -1]
A = As
C = [1]
inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets try add c1
C = [1 1]
inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets to add a1 a2
A = conv(As, [1 1 1])

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);
% The residual is white when doing whiteness test! Both Monti and Spectrum.
%% Lets remove C1, not sig
C = [1 ]

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);
% The residual is white when doing whiteness test! Both Monti and Spectrum.
%% add c12
C = [1 zeros(1,11) 1]

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% add c36
C = [1 zeros(1,11) 1 zeros(1,23) 1]

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% add a3
A = conv(As, [1 1 1 zeros(1, 30) 1])

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);
%% LET US CONTINUE TO GET K = 7 GOOD. 
% remove C3, it was not sig. 
C = [1 0 1]

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%%
% add a2
A(13) = 1

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);
% How we are super white, and all coeffs are significant. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the BJ model.
%inputModel.A = conv([1 zeros(1, 36-1) -1], inputModel.A);   % Add the differentiation to the model.
%ex = filter( inputModel.A, inputModel.C, modelx );              % Lets look at the residaual

figure
plot( ex )                                                  % As there are lots of roots on the unit circle, the ringing lasts longer than expected!
ylabel('Pre-filtered input signal')
xlabel('Time')
title('Prolonged ringing due to root on the unit circle')


%% Lets remove some additional samples just to be safe...
cutoff = 10
ey = filter( inputModel.A, inputModel.C, modely );   
ey = ey(length(inputModel.A)+cutoff:end );
% ex = ex(length(inputModel.A)+cutoff:end );                      % Remove some more samples given the ringing (this is much more than needed).
ex =   filter(inputModel.A, inputModel.C, modelx); 
ex = ex(length(inputModel.A)+cutoff:end );
%%
var_ex = var(ex);

figure;
[Cxy,lags] = xcorr( ey, ex, noLags, 'coeff' );
stem( lags, Cxy )
hold on
condInt = 2*ones(1,length(lags))./sqrt( length(ey) );
plot( lags, condInt,'r--' )
plot( lags, -condInt,'r--' )
hold off
xlabel('Lag')
ylabel('Amplitude')
title('Crosscorrelation between filtered in- and output')
% d = 2 or potentially 0?, maybe r = 0 and s = 2

%% Lets form an initial model.
% The function call is estimateBJ( y, x, C1, A1, B(d+s), A2(r), titleStr, noLags )
[ foundModel, ey, ~, pacfEst ] = estimateBJ( modely, modelx, [1], [1 1], [0 0 0 1 1], [1  1 1], 'BJ model 1', noLags );
var_ey = var(ey);

%% Lets form an initial model if we are using yd
% The function call is estimateBJ( y, x, C1, A1, B(d+s), A2(r), titleStr, noLags )
[ foundModel, ey, ~, pacfEst ] = estimateBJ( modely, modelx, [1 1], [1 1], [0 1 1], [1], 'BJ model 1', noLags );
var_ey = var(ey);
%% Better... Maybe add a c2 term?
% [ foundModel, ey, ~, pacfEst ] = estimateBJ( modely, modelx, [1 0 1], [1 1], [0 0 0 0 1], [1], 'BJ model 3', noLags );
% var_ey = var(ey);


%% We now have a white residual; can we trust the Monti test?
checkIfNormal( pacfEst(2:end), 'PACF' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lets predict the input first.
close all
k = 7;
[Fx, Gx] = polydiv( inputModel.C, inputModel.A, k );
xhatk = filter(Gx, inputModel.C, valx);

% -----------inverse transform before plotting anything-------------- %
% xhat = exp(xhatk)-1;   x = exp(valx)-1;   
xhat = xhatk;   x = valx  

ehat = x - xhat;
ehat = ehat(modelLim+1:end);

% plot prediction
plot([x xhat])
xline(modelLim,'--');
title('Prediction, k = 1, model 2')
legend('valx','xhatk','validation start')
figure
plot(ehat)
legend('ehat')
MSE = mean(ehat.^2)
%% Naive predictor k =1: 
k=1 % month i rains as much as in month i-1
% k=7 %, month i rains as much as in month i-36, i.e. last year
naive_x = zeros(length(x),1);
for i=2:length(x)
    naive_x(i) = x(i-1);
end
e_naive = x - naive_x;
e_naive = e_naive(length(modely)+1:end);
figure
plot([x naive_x])
legend('valx','naive x')
MSEN = mean(e_naive.^2)
%% Naive predictor k = 7: 
k=1 % month i rains as much as in month i-1
k=7 %, month i rains as much as in month i-36, i.e. last year
naive_x = zeros(length(x),1);
for i=37:length(x)
    naive_x(i) = x(i-36);
end
e_naive = x - naive_x;
e_naive = e_naive(length(modely)+1:end);
figure
plot([x naive_x])
legend('valx','naive x')
MSEN = mean(e_naive.^2)
%% analyse input prediction
clc; close all
% calculate variances for prediction residual
var_ehat = var(ehat)
% Normalized prediction variance, should be < 1
norm_var = var(ehat)/var(x(modelLim+1:end))

% calculate variances for prediction
var_naive = var(e_naive)
norm_var_naive = var(e_naive)/var(x(modelLim+1:end))

% plot the acf of prediction and naive prediction
figure
acf( ehat, 100, 0.05,1,0,false); title('ACF of prediction residual, model 2, k=7')
figure
acf( e_naive, 100, 0.05,1,0,false); title('ACF of naive prediction residual, model 2, k=7')

% Check if normal with D'Agostino-Pearson's K2 test
acfEst = acf( ehat, 100, 0.05 );
checkIfNormal( acfEst(k+1:end), 'ACF' );

%% Predict NVDI, yhatk, using both x, xhatk, and y
close all
k = 7
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, valx)+ filter(Gy, KC, valy);

% -----------inverse transform before plotting anything-------------- %
% yhat = exp(yhatk)-1;   y = exp(valy)-1; % Dont really need to transform
% back yet (and when we do, we need to transform back the naive one also i
% think, to make comparisons. 
yhat = yhatk;   y = valy
% check prediction error
ehat_y = y - yhat;
ehat_y = ehat_y(modelLim+1:end);

figure
plot([y yhat])
xline(length(modelx),'--');
legend('y', 'yhat', 'validation start')

figure
plot(ehat_y)
legend('ehat_y')
MSE = mean(ehat_y.^2)
%% Naive y predictor: 
k=1%, month i rains as much as in month i-1
% k=7, month i rains as much as in month i-36, i.e. last year
naive_y = zeros(length(y),1);
for i=2:length(y)
    naive_y(i) = y(i-1);
end
e_naive_y = y - naive_y;
e_naive_y = e_naive_y(modelLim+1:end);
plot([y naive_y])
legend('valy','naive y')
MSEN = mean(e_naive_y.^2)
%% Naive y pred, k = 7
% k=1, month i rains as much as in month i-1
k=7 %, month i rains as much as in month i-36, i.e. last year
naive_y = zeros(length(y),1);
for i=37:length(y)
    naive_y(i) = y(i-36);
end
e_naive_y = y - naive_y;
e_naive_y = e_naive_y(modelLim+1:end);
plot([y naive_y])
legend('valy','naive y')
MSEN = mean(e_naive_y.^2)
%% analyse prediction
clc; close all
% calculate variances for prediction residual
var_ehat_y = var(ehat_y)
% Normalized prediction variance, should be < 1
norm_var_y = var(ehat_y)/var(y(length(modely)+1:end))

% calculate variances for 
var_naive_y = var(e_naive_y)
norm_var_naive_y = var(e_naive_y)/var(y(length(modely)+1:end))

% plot the acf of prediction and naive prediction
figure
acf( ehat_y, 100, 0.05,1,1,false); title('ACF of prediction residual, model 2, k=7')
figure
acf( e_naive_y, 100, 0.05,1,0,false); title('ACF of naive prediction residual, model 2, k=7')

% Check if normal with D'Agostino-Pearson's K2 test
acfEst = acf( ehat_y, 100, 0.05 );
checkIfNormal( acfEst(k+1:end), 'ACF' );


%%
%Restore
yres = filter([1], [1 1], yhat)
figure(1)
%%
plot([1:length(y)], [1:length(y)]*a.p1 + a.p2)
hold on
plot(y)
hold off