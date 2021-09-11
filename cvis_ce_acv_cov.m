function [cov_Q0Q1, rho] = cvis_ce_acv_cov()

    format long;
    rng('default');

    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2.8;
    l0 = mu + w0*std;
    l1 = mu + w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;
    
%     prob0 = 1-normcdf(l0);
%     prob1 = 1-normcdf(l1);
        
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 3000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
        
    %%
    n0_CV = 10^6;
    n1_CV = n0_CV;
        
    s_MC = random(gm, n0_CV);
    Q0s_CV = Q0(s_MC(1:n0_CV))<0;
    w0_CV = mvnpdf(s_MC(1:n0_CV),mu,std)./qce(s_MC(1:n0_CV));
    wQ0s_CV = w0_CV.*Q0s_CV;
    Q1s_CV = Q1(s_MC(1:n1_CV))<0;
    w1_CV = mvnpdf(s_MC(1:n1_CV),mu,std)./qce(s_MC(1:n1_CV));
    wQ1s_CV = w1_CV.*Q1s_CV;
    cov_Q0Q1_MC = cov(wQ0s_CV,wQ1s_CV);
    cov_Q0Q1 = cov_Q0Q1_MC(1,2)
    
    rho = cov_Q0Q1_MC(1,2)/sqrt(cov_Q0Q1_MC(1,1)*cov_Q0Q1_MC(2,2))
    
end