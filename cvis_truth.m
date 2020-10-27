function cvis_truth()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    w0 = 3;
    l0 = mu + w0*std;
    
    Q0 = @(x) l0-x;
    
    ns = 10000000;
    samples = randn(ns,1);
    model = Q0(samples(:))<0;
    m = mean(model);
    v = var(model);
    disp(1-normcdf(l0));
    disp(m);
    disp(v);
        
end