%% Filter process that omits number of samples (equal to order of ar-part of model). Runs the inverse filter of my model. 
function Filtered = MyFilter(poly_a,poly_c, y)
    e_hat=filter(poly_a,poly_c,y)
    Filtered=e_hat(length(poly_a):end)
end