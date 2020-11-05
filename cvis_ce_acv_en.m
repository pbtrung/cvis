function cvis_ce_acv_en()

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
    
    prob0 = 1-normcdf(l0);
    prob1 = 1-normcdf(l1);
    
    a = linspace(-1.5,0.5,33);
    
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
    
    n_MC = 5*10^5;
    s_MC = random(gm, n_MC);
    Q0s_MC = Q0(s_MC(:))<0;
    w_MC = mvnpdf(s_MC,mu,std)./qce(s_MC);
    wQ0s_MC = w_MC.*Q0s_MC;
    
    nbatches = 1000;
    
    % c0 = c*c1
    c = 10;
    
    %%
    n0_CV = fix(n_MC*c/(c+1));
    n1_CV = n0_CV;
    
    n0_CV_EN = fix(n0_CV/nbatches);
    n1_CV_EN = fix(n1_CV/nbatches);
    
    CV_EN(1:nbatches) = 0;
    for i = 1:nbatches
        s_CV_EN = random(gm, n0_CV_EN);
        Q0s_CV_EN = Q0(s_CV_EN(:))<0;
        w0_CV_EN = mvnpdf(s_CV_EN,mu,std)./qce(s_CV_EN);
        wQ0s_CV_EN = w0_CV_EN.*Q0s_CV_EN;
        Q1s_CV_EN = Q1(s_CV_EN(:))<0;
        w1_CV_EN = mvnpdf(s_CV_EN,mu,std)./qce(s_CV_EN);
        wQ1s_CV_EN = w1_CV_EN.*Q1s_CV_EN;
        
        cov_Q0Q1 = cov(wQ0s_CV_EN,wQ1s_CV_EN);
        cov_MC_Q0Q1 = cov_Q0Q1/n0_CV_EN;
        astar = -cov_MC_Q0Q1(1,2)/cov_MC_Q0Q1(2,2);
        CV_EN(i) = mean(wQ0s_CV_EN) + astar*(mean(wQ1s_CV_EN)-prob1);
    end
    
    v0_MC = var(wQ0s_MC)/n_MC;
    v_CV_EN = var(CV_EN)/nbatches;
    disp('v_CV_EN/v0_MC');
    disp(v_CV_EN/v0_MC);
    
    %%
    r = 2.5;
    n0_IS = fix(n_MC*c/(c+r+1));
    n1_IS = n0_IS;
    m1_IS = fix(r*n0_IS);
    r1 = (n1_IS+m1_IS)/n0_IS;
    
    n0_IS_EN = fix(n0_IS/nbatches);
    n1_IS_EN = fix(n1_IS/nbatches);
    m1_IS_EN = fix(m1_IS/nbatches);
    
    IS_EN(1:nbatches) = 0;
    for i = 1:nbatches
        s_IS_EN = random(gm, n0_IS_EN);
        Q0s_IS_EN = Q0(s_IS_EN(:))<0;
        w0_IS_EN = mvnpdf(s_IS_EN,mu,std)./qce(s_IS_EN);
        wQ0s_IS_EN = w0_IS_EN.*Q0s_IS_EN;
        Q1s_IS_EN = Q1(s_IS_EN)<0;
        w1_IS_EN = mvnpdf(s_IS_EN,mu,std)./qce(s_IS_EN);
        wQ1s_IS_EN = w1_IS_EN.*Q1s_IS_EN;
        
        s_mu_IS_EN = random(gm, m1_IS_EN);
        mu1_IS_EN = Q1(s_mu_IS_EN)<0;
        w1_mu_IS_EN = mvnpdf(s_mu_IS_EN,mu,std)./qce(s_mu_IS_EN);
        wmu_IS_EN = w1_mu_IS_EN.*mu1_IS_EN;
        mu1_IS = (sum(wQ1s_IS_EN) + sum(wmu_IS_EN))/(n1_IS_EN + m1_IS_EN);
        
        cov_Q0Q1 = cov(wQ0s_IS_EN,wQ1s_IS_EN);
        cov_MC_Q0Q1 = cov_Q0Q1/n0_IS_EN;
        astar = -cov_MC_Q0Q1(1,2)/cov_MC_Q0Q1(2,2);
        IS_EN(i) = mean(wQ0s_IS_EN) + astar*(mean(wQ1s_IS_EN)-mu1_IS);
    end
    
    v_IS_EN = var(IS_EN)/nbatches;
    disp('v_IS_EN/v0_MC');
    disp(v_IS_EN/v0_MC);
    
end