% Computer exercise 1


load("svedala.mat")
load("data.dat")
load("noise.dat")
%% Polynomials 
A1 = [ 1 -1.79 0.84 ];
C1 = [ 1 -0.18 -0.11 ];

A2 = [ 1 -1.79 ];
C2 = [ 1 -0.18 -0.11 ];

%% #Create the ARMA using A1 and C1, (and A2 and C2 respectively). 
ARMA1 = idpoly( A1, [ ] , C1 );
ARMA2 = idpoly( A2, [ ] , C2 );

%% plot
pzmap(ARMA2)
%% set seed 

rng(0)

%% this was used to create arma_signal
sigma2=1.5
N=1000 
e = sqrt ( sigma2 ) * randn( N, 1 ) ;
%Find the process as filtered noise
y = filter(ARMA1.c, ARMA1.a, e ) ;
%remove the first 100
y = y(101:end)

plot(y)
%% simulate signal 1 
sigma2=1.5
N=300
y1 = ARMA_signal(N,sigma2,ARMA1,1);
%% Simulate signal 2
y2 = ARMA_signal(N,sigma2,ARMA2,1)

%% Plot: Question 1. We can see that the second process has root outside unit circle.
subplot (211)
plot ( y1 )
subplot (212)
plot ( y2 )
%% Real and estimated covariance arma1 /
m =20
[real_cov1, tau1] = kovarians(ARMA1.c,ARMA1.a, m)
%[real_cov2, tau2] = kovarians(ARMA2.c,ARMA2.a, m)

r_est1 = covf(y1,m+1)

%% Plot kov and est kov. Becuase we assume unit variance of the driving noise, combined with too little data. 
stem( 0:m, real_cov1*sigma2 )
hold on
stem( 0:m, r_est1, 'r' )

%% 
[acf1,pafc1]= acf_pacf_norm(y1)
%%
%normplot(acf1)

%% Try Different order ARs:
data_y1= iddata(y1);
ar2 = arx(data_y1, 2);
error_ar1 =MyFilter(ar2.A,ar2.C,y1);
[autocorrAR,pautocorrAR]=acf_pacf_norm(error_ar1);
%%
present(ar2)
%FPE 1.587 for AR3
%% Whiteness test
whitenessTest(error_ar1)
% It seems white enough
%% Arma model
arma_model_1 = armax ( y1 , [ 2, 1] ) ;
error_arma1 = MyFilter(arma_model_1.A, arma_model_1.C, y1)
[autoARMA,pautoARMA]=acf_pacf_norm(error_arma1)
present(arma_model_1)

%% Whiteness test
whitenessTest(error_arma1)
%% 2.2

data = iddata(data)
noise = iddata(noise)

acf_pacf_norm(data.y)

%% Create ar(1) to ar(5) models
armod1 = arx(data,1)
armod2 = arx(data,2)
armod3 = arx(data,3)
armod4 = arx(data,4)
armod5 = arx(data,5)

%% Create residual
rar1=resid(armod1, data);
rar2=resid(armod2,data); 
rar3=resid(armod3,data);
rar4=resid(armod4,data); 
rar5=resid(armod5,data);
%% (we can see that AR2 is not quite white, but AR3 is). 
acf_pacf_norm(rar3.y)
whitenessTest(rar3.y)
%% Plot noise vs residuals:
subplot(611)
plot(noise.y(cutoff:end))
subplot(612)
plot(rar1.y(cutoff:end))
subplot(613)
plot(rar2.y(cutoff:end))    
subplot(614)
plot(rar3.y(cutoff:end))
subplot(615)
plot(rar4.y(cutoff:end))
subplot(616)
plot(rar5.y(cutoff:end))

%% Now try with arma models:
arma11 = armax(data, [ 1 1 ] )
arma12 = armax(data, [ 1 2 ] )
arma21 = armax(data, [ 2 1 ] )
arma22 = armax(data, [ 2 2 ] )

%% resid of arma
rarm1=resid(arma11, data);
rarm2=resid(arma12,data);
rarm3=resid(arma21, data);
rarm4=resid(arma22,data); 
%QUESTION: Should i have a cutoff here?
acf_pacf_norm(rarm1.y(2:end));
whitenessTest(rarm1.y(2:end));

% Conclusion: Since both arma11 and armod 3 are valid, id choose the one
% with the lowest FPE (arma11)

%% 2.3 Simulate arma data and try to deal with seasonality
rng( 0 )
N=600
A = [ 1 -1.5 0.7 ] ;
C = [ 1 zeros(1 , 11) -0.5] ;
A12 = [ 1 zeros( 1 , 11 ) -1] ;
A_star = conv(A,A12) ;
e = randn( N , 1 ) ;
y = filter(C, A_star , e ) ;
y = y(101:end) ;
plot( y )
acf_pacf_norm(y)

%% Remove the 12the seasonal dependance
y_s = MyFilter (A12 , 1 , y ); 
datays = iddata(y_s );
acf_pacf_norm(datays.y)
%Can see 1st and 2nd coeff in the PACF, add these to A-poly
%% N 
modelinit=idpoly([1 0 0],[],[]);
model_armax = pem(datays ,modelinit)
resid1 = resid(model_armax,datays)
acf_pacf_norm(resid1.y(length(modelinit.A):end))
% at lag 12 in both ACF and PACF
%% Try adding the 12th MA coeff:
modelinit=idpoly([ 1 0 0 ],[],[1 zeros( 1 , 12 ) ] ) ;
modelinit.Structure.c.Free = [ zeros( 1 , 12 ) 1 ] ;
model_armax = pem(datays , modelinit )
resid2 = resid(model_armax,datays)
acf_pacf_norm(resid2.y(length(modelinit.A):end))
present(model_armax)
whitenessTest(resid2.y(length(modelinit.A):end))
%It did remove the season and the reisduals are white. The parameters are significant.   

%NN

%% Create an init model (AR) - not remoed seasonal dependence (rel)
datay=iddata(y)
A = [1,  zeros(1 , 14)] % Had to include alot of coeffecicients to get a good FPE, but all coeffs signif.
B=[]
C = [ 1 zeros(1 , 12) ]
%%
model_init = idpoly(A,B,C) ;
model_init.structure.a.Free = [ 1  1 1 zeros(1 , 9) 1 1 1 ] ;
model_init.structure.c.Free = [ zeros(1 , 12) 1] ;
model_armax_nrs = pem( datay , model_init )


%% Residual of armax (above) (rel
error_armax_nrs = resid(model_armax_nrs,datay)
acf_pacf_norm(error_armax_nrs.y(length(A):end)) % can see the 12th seasonal factor.
present(model_armax_nrs)

resid3 = resid(model_armax_nrs, datay)
resid3_test = filter(model_armax_nrs.a,model_armax_nrs.c, datay.y) %This gives diff from resid..?
whitenessTest(resid3.y(length(A):end)) %white and all coefficients are significant!


%% 2.4

datas= iddata(svedala)
[autos,pautos] = acf_pacf_norm(datas.y)
A24 = [ 1 zeros( 1 , 23 ) -1 ] ;
%% Remove season DIDNT SEEM TO DO THAT MUCH?
datas_s = filter(A24 , 1 , datas.y ) ;
datas_s=datas_s( length(A24 ) : end );

[autos2,pautos2] = acf_pacf_norm(datas_s)


%% First try with ar up to 2
% modelinit = idpoly ( [1 0 0] , [ ] , [] ) ;
% modelfin = pem( datas_s , modelinit )
% error_modelfin_season = resid(modelfin, datas_s)
% [autoseas,pautoseas]=acf_pacf_norm(error_modelfin_season(length(A):end))
% %Didnt work completely, add a factor for the seasonality


%% try a few things 
modelinit=idpoly([1 zeros(1,7)], [], [1 zeros(1, 24)])
modelinit.Structure.a.Free=[1 1 1 0 0 1 0  1]
modelinit.Structure.c.Free=[1 zeros(1,23) 1]

modelfin = pem(datas_s, modelinit)
error_modelfin_season = resid(modelfin, datas_s)
error_modelfin_season = filter( modelfin.a, modelfin.c, datas_s ); error_modelfin_season = error_modelfin_season(length(modelfin.a):end)
figure(1)
[autoseas,pautoseas]=acf_pacf_norm(error_modelfin_season);
present(modelfin)
figure(2)
whitenessTest(error_modelfin_season)


%% Next try svedala (non diff data)
A = [1 zeros(1, 17)]
C= [1 zeros(1,24)]
modelinit=idpoly(A, [], C)
modelinit.Structure.a.Free=[1 1 1 1 1 zeros(1,12) 1 ]
modelinit.Structure.c.Free=[1 0 0 1 zeros(1, 6) 1 zeros(1, 12) 1 1]
modelfin = pem(datas, modelinit)
error_modelfin_season = resid(modelfin, datas); error_modelfin_season(15:end);
%error_modelfin_season = filter( modelfin.a, modelfin.c, datas.y); error_modelfin_season = error_modelfin_season(length(modelfin.a):end)
[autoseas,pautoseas]=acf_pacf_norm(error_modelfin_season.y)
present(modelfin)
figure(4)
whitenessTest(error_modelfin_season.y)

