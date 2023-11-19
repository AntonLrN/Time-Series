%% Functions
function f = ARMA_signal(N,sigma2, ARMA_model, plotting)
    rng ( 0 )
    e = sqrt ( sigma2 ) * randn( N, 1 ) ;
    y = filter(ARMA_model.c, ARMA_model.a, e );
    f=y(101:end)
    if plotting ==1
        plot(f)
    end
end
