
addpath("C:\Users\anton\Desktop\University\Time series analysis\Computer Exercises\Exercise 1\CourseMaterial\CourseMaterial\Code\data")
addpath("C:\Users\anton\Desktop\University\Time series analysis\Computer Exercises\Exercise 1\CourseMaterial\CourseMaterial\Code\functions")
addpath("C:\Users\anton\Desktop\University\Time series analysis\Computer Exercises\Exercise 1")
load('tork.dat')
%%
rng ( 0 )

n = 500; % Number of samples
A3 = [ 1 .5 ] ;
C3 = [ 1 -.3 .2 ] ;
w = sqrt( 2 ) * randn(n + 100 , 1 ) ;
x = filter(C3 ,A3 ,w) ; % Create the input
A1 = [ 1 -.65] ;
A2 = [ 1 .90 .78 ] ;
C = 1 ;
B = [ 0 0 0 0 .4 ] ;
e = sqrt ( 1.5 ) * randn(n + 100 , 1 ) ;
y = filter (C,A1 , e ) + filter (B,A2 , x ) ; % Create the output
x = x (101 : end) , y = y ( 101 : end) % Omit i n i t i a l samples
clear A1 , A2 , C, B, e , w, A3 , C3

%% Creating ARMA model for input: (Question 1)
arma_input= armax ( x , [ 1 0 ] ) ;
w_t = filter(arma_input.a,arma_input.c,x); w_t = w_t(length(arma_input.a):end) %x = -1 might be a thing, but the coeff right after is not significant. Then we have 2 significant then 2 not significant then ringing again. Prolly 4,2,0
acf_pacf_norm(w_t)


eps_t =  filter(arma_input.a,arma_input.c,y); eps_t = eps_t(length(arma_input.a):end)
whitenessTest(w_t)
%%
M = 40; 
stem(-M:M, crosscorr(w_t,eps_t, M)); 
 title('Cross_correlation_function')

 xlabel('lag')
 hold on
 plot(-M:M, 2/sqrt(n)*ones(1,2*M+1), '--')
 plot(-M:M, -2/sqrt(n)*ones(1,2*M+1), '--')

%% Use d r s coeffs in A2 and B:
A2=[1 1]
B = [0 0 0 0 ]
Mi = idpoly( [ 1 ] , [B] , [ ] , [ ] , [ A2 ] ) ;
z = iddata (y , x ) ;
Mba2 = pem( z ,Mi ) ; present(Mba2)
etilde= resid(Mba2 , z ); etilde = etilde(length(A2):end)
crosscorr(etilde.y,x,40) % not quite uncorrelated but quite cloe
%% 

arma_etilde= armax(etilde.y, [3 2])
error_etilde = filter(arma_etilde.a, arma_etilde.c, etilde.y)
acf_pacf_norm(error_etilde(length(A2):end))
whitenessTest(error_etilde(length(A2):end))


%% asd asd asd asd

tork = tork - repmat(mean( tork ) , length ( tork ) , 1 ) ;
y = tork( : , 1 ) ; x = tork ( : , 2 ) ;
z = iddata(y , x ) ;
plot( z ( 1 : 300 ) )

%%
arma_input2= armax ( x , [ 1 0 ] ) ;
w_t = filter(arma_input2.a,arma_input2.c,x); w_t = w_t(length(arma_input2.a):end)%x = -1 might be a thing, but the coeff right after is not significant. Then we have 2 significant then 2 not significant then ringing again. Prolly 4,2,0
acf_pacf_norm(w_t)


eps_t =  filter(arma_input2.a,arma_input2.c,y); eps_t = eps_t(length(arma_input2.a):end)
whitenessTest(w_t)

%%
M = 40; 
stem(-M:M, crosscorr(w_t,eps_t, M)); 
 title('Cross_correlation_function')

 xlabel('lag')
 hold on
 plot(-M:M, 2/sqrt(n)*ones(1,2*M+1), '--')
 plot(-M:M, -2/sqrt(n)*ones(1,2*M+1), '--') % d r s = 2 1 2 it seems from the book. 

 %%
A2=[1 ]
B = [0 0 1 1 ]
Mi = idpoly( [ 1 ] , [B] , [ ] , [ ] , [ A2 ] ) ;
z = iddata (y , x ) ;
Mba2 = pem( z ,Mi ) ; present(Mba2)
etilde= resid(Mba2 , z ); etilde = etilde(length(A2):end)
%% 
arma_etilde= armax(etilde.y, [1 1])
error_etilde = filter(arma_etilde.a, arma_etilde.c, etilde.y)
acf_pacf_norm(error_etilde(length(A2):end))
whitenessTest(error_etilde(length(A2):end))

%%  seems like d=4 not d=3 works better. And increasing the A polys with 1 each. Resid is not uncorr from input data, but resid is kidna white. 
A1 = [1 1]
A2=[1 1]
C = [1]
B = [0 0 0 0 1 1 ]
Mi = idpoly ( 1 ,B,C,A1 ,A2 ) ;
z = iddata(y , x ) ;
MboxJ = pem( z ,Mi ) ;

ehat = resid(MboxJ , z ) ;


acf_pacf_norm(ehat.y(length(A2):end))
present(MboxJ)





