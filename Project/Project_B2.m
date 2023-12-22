%% Use reconstructed rain as input model to predict NVDI
% clear; 
%run Project;
% close all; clc;
load proj23.mat
% load rain_3_est.mat
noLags = 140;
%%
nvdi = ElGeneina.nvdi;
nvdi = nvdi/255*2 - 1;      % rescale to [-1,1]
%% Should we transform the data?
x = x_st_vec(22*36+1:end)';
%This to train on 
cutoff= 300
x_input = x_st_vec(cutoff:end)';
% x and y need to be in phase
y = nvdi;           

figure
lambda_max = bcNormPlot( x_input );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
figure
lambda_max = bcNormPlot( y );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
xlabel('Lambda')
ylabel('Log-likelihood')

% Both seems to indicate that a log-transform could be helpful.
x_input = log(x_input+abs(min(x_input))+1);
% x=log(x+abs(min(x))+1);
y=log(y+1);

%% extract modeling and validation sets and test
modelLim = 36*11%36*13;        % We select 13 years. Does the data seem stable enough?
x = x_input(22*36+1 - (cutoff-1):end)';
% sets for y
modely = y(1:modelLim);  %modely = modely-mean(modely); 
valy = y(1:36*14); 
testy = y(1:36*16)

modelx = x(1:modelLim);  %modelx = modelx-mean(modelx);
valx = x(1:36*14);
testx = x(1:36*16);
%% Create a model for the input.
plotACFnPACF(x_input, noLags, 'Input, x_t' ); %From this we can see a very clear season of 36

%% There seems to be a strong periodicity at 36, suggesting that a differentiation might help.
%estimateARMA( modelx, [ 1 zeros(1,35) 1], [ 1 ], 'Differentiated input, version 2', noLags );



%% Lets try add a1 since we know that should be there
A = [1 1]
C = [1]
inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets try add c 2 c3
C = [1 0 1 1]
inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets try add to conv with a36
A = conv(A, [1 zeros(1,35) -1])

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);
% The residual is white when doing whiteness test! Both Monti and Spectrum.
% 
%% Also add C36?
C = [1 0 0 1 zeros(1,32) 1]

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets try add a3?
A(4) = 1 

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%% Lets remove C3 since it is almost not sig. 
C(4) = 0

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% Add C1

C(2) = 1 

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% Add A4

A(5) = 1 

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% Add C12

C(13) = 1

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);

%% Add A5

A(6) = 1

inputModel = estimateARMA( x_input', A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, x_input' ); ex(length(A):end);
acf_pacf_norm(ex);


%White accoring to spectrum test
%White accoring to spectrum test
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

%% Lets form an initial model.
% The function call is estimateBJ( y, x, C1, A1, B(d+s), A2(r), titleStr, noLags )
estimateBJ( modely, modelx, [1], [1 1], [0 0 1 1], [1 1], 'BJ model 1', noLags );

%% Better... Maybe add a c2 term?
[ foundModel, ey, ~, pacfEst ] = estimateBJ( modely, modelx, [1 0 1], [1 1], [0 0 0 0 1], [1], 'BJ model 3', noLags );
var_ey = var(ey);


%% We now have a white residual; can we trust the Monti test?
checkIfNormal( pacfEst(2:end), 'PACF' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lets predict the input first.
close all
k = 7;
[Fx, Gx] = polydiv( inputModel.C, inputModel.A, k );
xhatk = filter(Gx, inputModel.C, valx);

% -----------inverse transform before plotting anything-------------- %
xhat = exp(xhatk)-1;   x = exp(valx)-1;    

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

%% Naive predictor: 
% k=1 % month i rains as much as in month i-1
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
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, valx) + filter(Gy, KC, valy);

% -----------inverse transform before plotting anything-------------- %
yhat = exp(yhatk)-1;   y = exp(valy)-1;

% check prediction error
ehat_y = y - yhat;
ehat_y = ehat_y(k+20:end);

figure
plot([y yhat])
xline(length(modelx),'--');
legend('y', 'yhat', 'validation start')

figure
plot(ehat_y)
legend('ehat_y')
%% Naive y predictor: 
% k=1, month i rains as much as in month i-1
% k=7, month i rains as much as in month i-36, i.e. last year
naive_y = zeros(length(y),1);
for i=2:length(y)
    naive_y(i) = y(i-1);
end
e_naive_y = y - naive_y;
e_naive_y = e_naive_y(modelLim+1:end);
plot([y naive_y])
legend('valy','naive y')
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