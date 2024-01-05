clear; clc; close all
load proj23.mat
nvdi = ElGeneina.nvdi;

min_nvdi = min(nvdi);
max_nvdi = max(nvdi);

normalized_nvdi = (nvdi - min_nvdi) / (max_nvdi - min_nvdi);


nvdi_s= 2 * normalized_nvdi - 1;

nvdi=nvdi_s
plot(nvdi,'g')
title('rescaled NDVI in El Geneina')
%% Check if transform is needed
figure; 
lambda_max = bcNormPlot(nvdi,1);    % log-transform kan vara lämplig
original = nvdi
%%
nvdi_log = log(nvdi + 1 - min(nvdi))
nvdi_log = nvdi_log - mean(nvdi_log)

% figure
% normplot( original )                % verkligen inte gausian
figure
normplot(nvdi_log)                      % mycket mer gausian!
figure
plot(nvdi_log)

%% Extract a modelling set.
              
modely = nvdi(1:36*12);   % We select 13 years. Does the data seem stable enough?
modely = modely - mean(modely);
figure
plot(modely)		
title('log-transformed nvdi');

valy = nvdi(1:36*14);

%% analyse spectra: season of 36, 18, 12,
close all
Padd = 1024;                                     % Try increasing the zero-padding to 1024.
N = length(modely);
X = fftshift( abs( fft(modely, Padd) ).^2 / N );

ff = (0:Padd-1)'/Padd-0.5;

% Notice the loss of power, the increasing with the mainlobe, as well as
% the weaker sidelobes when using the windowing. 
Y = X;
Y2 = fftshift( abs( fft(modely.*hamming(N), Padd) ).^2 / N );
Y3 = fftshift( abs( fft(modely.*blackman(N), Padd) ).^2 / N );

figure
subplot(211)
plot( ff, [Y Y2 Y3])                
title('Frequency-domain')
legend('Rectangular','Hamming','Blackman')
subplot(212)
semilogy( ff, [Y Y2 Y3])
ylabel('Log scale')
xlabel('Frequency')

%% Examine ACF and PACF
% yearly season, i.e. 36 lags
close all
noLags = length(modely)/4;
plotACFnPACF( modely, noLags, 'Data' );
%%
A = [1 1]
C = [1]
model = estimateARMA( modely, A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, modely ); ex(length(A):end);
acf_pacf_norm(ex);

%% add c36
A = [1 1]
C = [1 zeros(1,35) 1]
model = estimateARMA( modely, A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, modely ); ex(length(A):end);
acf_pacf_norm(ex);

%% a36
A = [1 1 zeros(1, 34) 1]
C = [1 zeros(1,35) 1]
model = estimateARMA( modely, A, C, 'Differentiated input, version 2', noLags );
ex = filter( inputModel.A, inputModel.C, modely ); ex(length(A):end);
acf_pacf_norm(ex);



%% Differentiate to remove the periodicities.
% close all
% s36 = 36;
% yearPoly  = [1 zeros(1,s36-1) -1];
% modelData = filter(yearPoly,1,modelData);      modelData = modelData(length(s36):end);
% 
% plotACFnPACF( modelData, noLags, 'Differentiated data' );
% %% Make a first model
% close all
% A = [1 zeros(1,37)];
% C = [1 zeros(1,37)];
% y = iddata(modely);
% model_init = idpoly(A,[],C);
% model_init.Structure.A.Free = [1 1 zeros(1,34) 1 1];
% model_init.Structure.C.Free = [1 0 0 0 zeros(1,31) 0 0 0];
% model = pem(y,model_init);
% e = resid(model,y);      e = e(length(model.A):end );
% present(model)
% 
% % plotter(e.y,noLags)
% acf_pacf_norm(e.y)
% %figure
% %plot(e.y)
%% predict validation
close all;
y = valy;
k = 7;
A = model.A;
C = model.C;
[F, G] = polydiv( C, A, k );

yhatk = filter(G,C,y);
                                
ehat = y-yhatk;     % prediction residual
ehat = ehat(length(modely)+1:end);

% plot prediction
plot([y yhatk] )
xline(length(modely),'--');
title('Prediction, k = 1, model 2')
legend('y','yhat','validation start')
figure
plot(ehat)
legend('ehat')

figure
acf( ehat, noLags, 0.05,1);
title('ACF')
MSE = mean(ehat.^2)
%% Naive predictor: 
% k=1 %, month i rains as much as in month i-1
k=7%, month i rains as much as in month i-36, i.e. last year
naive_y = zeros(length(valy),1);
for i=37:length(valy)
    naive_y(i) = valy(i-36);
end
e_naive = y - naive_y;
e_naive = e_naive(length(modely)+1:end);
fprintf('  The variance of original signal is         %5.2f.\n', var(e_naive)')
MSE = mean(e_naive.^2)
%% analyse prediction
fprintf('Prediction the signal %i-steps ahead.\n', k)
fprintf('  The variance of original signal is         %5.2f.\n', var(y)')
fprintf('  The variance of the prediction residual is %5.2f.\n', var(ehat)')
if var(ehat)<var(y)
    fprintf('  Amount of signal that was predicted is     %5.2f %%.\n', (1-var(ehat)/var(y))*100)
else
    fprintf('  **** BEWARE: the prediction is not accurate!!! ****\n')
end

fprintf('Theoretical prediction error variance %5.2f.\n', sum(F.^2)*var(e.y))

pacfEst = pacf( ehat, 100, 0.05 );
checkIfNormal( pacfEst(k+1:end), 'PACF' );
