 %% Use reconstructed rain as input model to predict NVDI
clear; 
%run Project;
close all; clc;
load proj23.mat
load rain_3_est.mat
noLags = 140;
%%
nvdi = ElGeneina.nvdi;
nvdi = nvdi/255*2 - 1;      % rescale to [-1,1]
x = rain_3_est(22*36:end)';   % x and y need to be in phase
y = nvdi;           
%% Should we transform the data?
figure
lambda_max = bcNormPlot( x );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
figure
lambda_max = bcNormPlot( y );
fprintf('The Box-Cox curve is maximized at %4.2f.\n', lambda_max)
xlabel('Lambda')
ylabel('Log-likelihood')

% Both seems to indicate that a log-transform could be helpful.
x=log(x+1);
y=log(y+1);

%% extract modeling and validation sets
modelLim = 32*13;        % We select 13 years. Does the data seem stable enough?
%y = y-mean(y);
modely = y(1:modelLim);  modely = modely-mean(modely); 
valy = y(1:32*16); 

%x = x-mean(x);
modelx = x(1:modelLim); 
valx = x(1:32*16);
%% Create a model for the input.
plotACFnPACF( modelx, noLags, 'Input, x_t' );


%% There seems to be a strong periodicity at 36, suggesting that a differentiation might help.
AS = [1 zeros(1,35) -1]
A = [1 1]
A = conv(A, AS)
A(1) = 1
inputModel = estimateARMA( modelx, A, [ 1 ], 'Differentiated input, version 2', noLags );


%% Try dumb model
scaleA = 0.1
scaleB=0.2
inputModel = estimateARMA( modelx, [ 1/scaleA 1  ]*scaleA, [ 1/scaleB 0 0  1]* scaleB, 'Differentiated input, version 2', noLags );
present(inputModel)




%% Lets try add a1.
scaleA = 0.1
scaleB=0.2
inputModel = estimateARMA( modelx, [ 1/scaleA 1  1 1 zeros(1,31) 1 1 1 ]*scaleA, [ 1/scaleB 1 1 1 zeros(1,32) 1]*scaleB, 'Differentiated input, version 2', noLags );
present(inputModel)



%%  Try pem instead:
% scaleA = 0.1
% scaleB=0.2
% A = [ 1/scaleA 1  1 1 zeros(1,31) 1 1 1 ]*scaleA
% C = [ 1/scaleB 1 1 1 zeros(1,32) 1]*scaleB
% modelinit = idpoly(A, [],C)
% modelinit.Structure.A.Free = [ 1 1  1 1 zeros(1,31) 1 1 1 ]
% modelinit.Structure.C.Free = [ 1 1  zeros(1,34) 1]
% model = pem(modelx, modelinit)
% present(model)
%% We need c3, and perhaps c2?
%estimateARMA( diff_xM, [ 1 1 ], [1 0 0 1], 'Differentiated input, version 3', noLags );


%% Ok, add c2 too... Maybe an a2 term as well...
% Good, now it is white and all coefficients are significant.
%inputModel = estimateARMA( diff_xM, [ 1 1 1 ], [1 0 1 1], 'Differentiated input, version 4', noLags );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create the BJ model.
%inputModel.A = conv([1 zeros(1, 36-1) -1], inputModel.A);   % Add the differentiation to the model.
ex = filter( inputModel.A, inputModel.C, modelx );              % Lets look at the residaual

figure
plot( ex )                                                  % As there are lots of roots on the unit circle, the ringing lasts longer than expected!
ylabel('Pre-filtered input signal')
xlabel('Time')
title('Prolonged ringing due to root on the unit circle')


%% Lets remove some additional samples just to be safe...
ey = filter( inputModel.A, inputModel.C, modely );   
ex = ex(length(inputModel.A)+30:end );                      % Remove some more samples given the ringing (this is much more than needed).
ey = ey(length(inputModel.A)+30:end );
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
% Seems like we only need b0.


%% Lets form an initial model.
% The function call is estimateBJ( y, x, C1, A1, B, A2, titleStr, noLags )
estimateBJ( modely, modelx, [1], [1], [1], [1], 'BJ model 1', noLags );


%% Lets try to add a1 and a2
estimateBJ( modely, modelx, [1], [1 1], [1], [1], 'BJ model 2', noLags );


%% Better... Maybe add a c2 term?
% Yes, now it is white.
[ foundModel, ey, ~, pacfEst ] = estimateBJ( modely, modelx, [1], [1 1], [1], [1], 'BJ model 3', noLags );
var_ey = var(ey);


%% We now have a white residual; can we trust the Monti test?
checkIfNormal( pacfEst(2:end), 'PACF' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lets predict the input first.
k = 1;
N = length(valx);
[Fx, Gx] = polydiv( inputModel.C, inputModel.A, k );
xhatk = filter(Gx, inputModel.C, valx);

%% Naive predictor: 
% k=1, month i rains as much as in month i-1
naive_x = zeros(length(valx),1);
for i=2:length(valx)
    naive_x(i) = valx(i-1);
end

%% plot input prediction and naive prediction
figure
plot([valx xhatk naive_x])
line( [modelLim modelLim], [-1e6 1e6 ], 'Color','red','LineStyle',':' )
legend('Input signal', 'Predicted input', 'naive predictor','Prediction starts')
title( sprintf('Predicted input signal, x_{t+%i|t}', k) )
axis([1 N min(valx)*1.5 max(valx)*1.5])

std_xk = sqrt( sum( Fx.^2 )*var_ex );
fprintf( 'The theoretical std of the %i-step prediction error is %4.2f.\n', k, std_xk)


%% Form the residual. Is it behaving as expected? Recall, no shift here!
ehat = valx - xhatk;
ehat = ehat(modelLim:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step input prediction residual', k) )
fprintf('This is a %i-step prediction. Ideally, the residual should be an MA(%i) process.\n', k, k-1)
checkIfWhite( ehat );
pacfEst = pacf( ehat, 70, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );

close all
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

yhatk_naive  = filter(Fhh, 1, naive_x) + filter(Ghh, KC, valx) + filter(Gy, KC, valy);

figure
plot(exp([valy yhatk])-1)
xline(length(modelx),'--');
legend('Output signal', 'Predicted output', 'validation start')
%% check prediction error
ehat = valy - yhatk;
ehat = ehat(k+20:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );

% Examine the variance of the prediction residual in comparison to the
% variance of the data.
fprintf('Prediction the signal %i-steps ahead.\n', k)
fprintf('  The variance of original signal is         %5.2f.\n', var(valy)')
fprintf('  The variance of the prediction residual is %5.2f.\n', var(ehat)')
if var(ehat)<var(valy)
    fprintf('  Amount of signal that was predicted is     %5.2f %%.\n', (1-var(ehat)/var(valy))*100)
else
    fprintf('  **** BEWARE: the prediction is not accurate!!! ****\n')
end

%% Predict NVDI using both input and past data
close all
KA = conv( foundModel.D, foundModel.F );
KB = conv( foundModel.D, foundModel.B );
KC = conv( foundModel.F, foundModel.C );

[Fy, Gy] = polydiv( foundModel.C, foundModel.D, k );
[Fhh, Ghh] = polydiv( conv(Fy, KB), KC, k );

yhatk  = filter(Fhh, 1, xhatk) + filter(Ghh, KC, valx) + filter(Gy, KC, valy);

figure
plot(exp([valy yhatk yhatk_naive])-1)
xline(length(modelx),'--');
legend('Output signal', 'Predicted output', 'Naive', 'validation start')
%% check prediction error
ehat = valy - yhatk;
ehat = ehat(k+20:end);

figure
acf( ehat, noLags, 0.05, 1 );
title( sprintf('ACF of the %i-step output prediction residual', k) )
checkIfWhite( ehat );
pacfEst = pacf( ehat, noLags, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );

% Examine the variance of the prediction residual in comparison to the
% variance of the data.
fprintf('Prediction the signal %i-steps ahead.\n', k)
fprintf('  The variance of original signal is         %5.2f.\n', var(valy)')
fprintf('  The variance of the prediction residual is %5.2f.\n', var(ehat)')
if var(ehat)<var(valy)
    fprintf('  Amount of signal that was predicted is     %5.2f %%.\n', (1-var(ehat)/var(valy))*100)
else
    fprintf('  **** BEWARE: the prediction is not accurate!!! ****\n')
end



