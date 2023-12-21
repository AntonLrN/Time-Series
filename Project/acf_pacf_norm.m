%% ACF function and plot
function [autocorr,pautocorr] = acf_pacf_norm(y1)
    lag=55
    autocorr = acf(y1,lag);
    
    pautocorr = pacf(y1,lag);
    subplot(311)
    normplot(y1)
    subplot(312)
    stem([0:lag],autocorr)
    title("ACF")
    ylim([-1 1])
    xlabel('Lag')
    ylabel('Amplitude')
    yline(0)
    yline(2/sqrt(length(y1)))
    yline(-2/sqrt(length(y1)))
    subplot(313)
    stem([0:lag],pautocorr)
    title("PACF")
    ylim([-1 1])
    xlabel('Lag')
    ylabel('Amplitude')
    yline(0)
    yline(2/sqrt(length(y1)))
    yline(-2/sqrt(length(y1)))
end