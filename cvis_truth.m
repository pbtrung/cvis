function cvis_truth()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2;
    l0 = mu+w0*std;
    
    Q0 = @(x) l0-x;
    
    vs = 100000;
    mc(1:vs) = 0;
    parfor i = 1:vs
        s = mvnrnd(mu,std,1000000);
        mc(i) = mean(Q0(s(:))<0);
    end
    
    var(mc)
    mean(mc)
    prob0 = 1-normcdf(w0)
    prob1 = 1-normcdf(w1)
    
%     v = 1.352806663625048e-09;
%     prob0 = 0.001349898031630;
%     prob1 = 0.022750131948179;
    
end