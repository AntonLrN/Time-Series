rng(0)
n = 500; % Number of samples
A3 = [1 .5];
C3 = [1, -0.3, 0.2];
w = sqrt(2) * randn(n + 100, 1);
x = filter(C3, A3, w); % Create the input
A1 = [1, -0.65];
A2 = [1, 0.90, 0.78];
C = 1;
B = [0, 0, 0, 0, 0.4];
e = sqrt(1.5) * randn(n + 100, 1);
y = filter(C, A1, e) + filter(B, A2, x); % Create the output
x = x(101:end); 
y = y(101:end); % Omit initial samples
clear A1 A2 C B e w A3 C3
%% 2.1
% We start with creating a model for the input. I.e finding A3 and C3. 
arma_input = armax(x, [1 2]) %but it kinda works with an AR(1).. ? 
w_t = filter(arma_input.a,arma_input.c,x) ; w_t = w_t(length(arma_input.a):end) %
acf_pacf_norm(w_t) 
whitenessTest(w_t)
%% Now we look at the filtered y_t to find B and A2 (ie. (d,r,s) to find the order, then optimize. 
eps_t =  filter(arma_input.a,arma_input.c,y); eps_t = eps_t(length(arma_input.a):end)
M = 40; 
stem(-M:M, crosscorr(w_t,eps_t, M)); 
title('Cross_correlation_function')
xlabel('lag')
hold on
plot(-M:M, 2/sqrt(n)*ones(1,2*M+1), '--')
plot(-M:M, -2/sqrt(n)*ones(1,2*M+1), '--')
%From this we can see that d = 4, and r = 2. But why is s = 0?
%% Given A2 and B, we can create the model for the output and input, which should leave the C1/A1 * et

A2 = [1, 0, 0];
B = [0 0 0 0 0.4];
Mi = idpoly([1], B, [], [], A2); %Create the polynomial using B/A2. 
z = iddata(y, x); % Creates a time series thing with both input and output
Mba2 = pem(z, Mi); %Takes the actual poly with defined input and output relation and optimizes the params. 
present(Mba2);
etilde = resid(Mba2, z); etilde = etilde(length(A2):end)%Calculates the resid (including input and everything)
crosscorr(etilde.y,x,40); %Insanely white. So the output y, and the input x are modelled. 
%% Now we have etilde, which should be equal to C1/A1 * e, and should therefore be able to get fileterd out.
arma_etilde = armax(etilde.y, [1 0])
error_etilde = resid(arma_etilde, etilde.y)
acf_pacf_norm(error_etilde(length(A2):end))
whitenessTest(error_etilde(length(A2):end))

%%
% 
A1 = arma_etilde.a;
A2 = Mba2.f
B = Mba2.b
C = [];
Mi = idpoly(1, B, C, A1, A2);
z = iddata(y, x);
MboxJ = pem(z, Mi);
present(MboxJ);
ehat = resid(MboxJ, z);

%% Test the ehat, is it white?
figure(1)
acf_pacf_norm(ehat.y)
hold on 
figure(2)
whitenessTest(ehat.y)
hold off

%Very white!
%% Cross correlation.
M = 100; 
stem(-M:M, crosscorr(x,ehat.y, M)); 
title('Cross_correlation_function')
xlabel('lag')
hold on
plot(-M:M, 2/sqrt(n)*ones(1,2*M+1), '--')
plot(-M:M, -2/sqrt(n)*ones(1,2*M+1), '--')
%Somewhat white.

%% 2.2 Hairdryer data
%%%%%%%%%%%%%%%%%%%%%%

load('tork.dat')
tork = tork - repmat(mean(tork), length(tork), 1);
y = tork(:, 1); x = tork(:, 2);
z = iddata(y, x);
plot(z(1:300))
n=length(y)

%%
% We start with creating a model for the input. I.e finding A3 and C3. 
arma_input = armax(x, [1 0]) %but it kinda works with an AR(1).. ? 
w_t = resid(arma_input,x) ; w_t = w_t(length(arma_input.a):end) %
figure(1)
acf_pacf_norm(w_t) 
hold on
figure(2)
whitenessTest(w_t)
hold off
A3=arma_input.a
C3=arma_input.c
%% eps_t = A3/C3 *y = B/A * z-d * w_t
eps_t = resid(arma_input, y); eps_t = eps_t(length(A3):end)
M = 100; 
stem(-M:M, crosscorr(w_t,eps_t,M)); 
title('Cross correlation function')
xlabel('lag')
hold on
plot(-M:M, 2/sqrt(n)*ones(1,2*M+1), '--')
plot(-M:M, -2/sqrt(n)*ones(1,2*M+1), '--')
%From this, d=2, r=1, s=2. 
%% define B and A2
B = [0 0 0 1  1] % d = 2, s = 2
A2 = [1 0 0]
Mi = idpoly([1], B, [], [], A2); %Create the polynomial using B/A2. 
z = iddata(y, x); % Creates a time series thing with both input and output
Mba2 = pem(z, Mi); %Takes the actual poly with defined input and output relation and optimizes the params. 
present(Mba2);
etilde = resid(Mba2, z); etilde = etilde(length(A2):end) %Calculates the resid (including input and everything)
crosscorr(etilde.y, x,40); %Least white shit ive seen in my life
%% Get C1 and A1 (based on etilde, which is based on (d,r,s)
arma_etilde = armax(etilde.y, [1 0]);
error_etilde = resid(arma_etilde, etilde.y);
figure(1)
acf_pacf_norm(error_etilde(length(A2):end));
hold on
figure(2)
whitenessTest(error_etilde(length(A2):end));
hold off
%% 
A1 = arma_etilde.a;
A2 = Mba2.f
B = Mba2.b
C = arma_etilde.c;
Mi = idpoly(1, B, C, A1, A2);
z = iddata(y, x);
MboxJ = pem(z, Mi);
present(MboxJ);
ehat = resid(MboxJ, z);
%% Is ehat white? seems so. And all coeffs highly significant. 
figure(1)
whitenessTest(ehat.y);
hold on
figure(2)
acf_pacf_norm(ehat.y);
%% 2.3 Prediction of ARMA-proc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load svedala
y = svedala;
%% prediction for k = 3
A = [ 1 -1.79 0.84 ] ;
C = [ 1 -0.18 -0.11 ] ;
% First we need to find the G_k through the Diophantine equation:
k=3 %put k = 1 here to get the variance of the noise. 
cutoff = 4 %
[ Fk , Gk ] = polydiv ( C, A, k) ; %Divide C/A, to get an F part and a G part (times z^-k). Cant use resid here?
yhat_k = filter(Gk ,C, y );% Get the k step prediction
pred_error_k = y - yhat_k; pred_error_k = pred_error_k(cutoff:end)
acfEst = acf( pred_error_k, 20, 0.05, 1 );
checkIfNormal( acfEst(k+1:end), 'ACF' ); %aprntly normal.
sigma2 = var(pred_error_k) % for k = 1, 0.3754
sigmae_2 = 0.3754

%1. mean of pred error:
PE_mean = mean(pred_error_k)
%2. Theoretical prediction error variance: sigma*e * sqrt(Fk^2
PE_var_est = sigma2
PE_var_theo = sigmae_2*(sum(Fk.^2)) % For k = 1 this is simply sigma2
% not the same, model not perfect so noise is not the true white noise. 
%3. Confidence intervals



%% for k = 26
A = [ 1 -1.79 0.84 ] ;
C = [ 1 -0.18 -0.11 ] ;
% First we need to find the G_k through the Diophantine equation:
k=26
[ Fk , Gk ] = polydiv ( C, A, k) ; %Divide C/A, to get an F part and a G part (times z^-k)
yhat_k = filter(Gk ,C, y ); %Get the k step prediction
pred_error_k = y - yhat_k; pred_error_k = pred_error_k(k+1:end) 
acfEst = acf( pred_error_k, 40, 0.05, 1 );
checkIfNormal( acfEst(k+1:end), 'ACF' ); %aprntly normal ?!



