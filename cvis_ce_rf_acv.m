function cvis_ce_rf_acv()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx0 = 60;
    nely0 = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 60;
    Ly = 20;
    lower = 1;
    upper = 2;
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx0,nely0,2);

    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx0*nely0,neig);
    
    l0 = 118.923;
    l1 = 108.510;
    
    dof0 = 2*(nely0+1)*nelx0+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    force = 1;
    
    Q0 = @(x) l0-FE0(lx,ly,nelx0,nely0,dof0,force,x,P,L,lower,upper);
    Q1 = @(x) l1-FE1(lx,ly,nelx1,nely1,dof1,force,x,PsiE1,L,lower,upper);
        
    a = linspace(-1.5,0.5,33);
    
    % definition of the random variables
    d      = neig;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 5;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    mu = zeros(1,neig);
    stdd = eye(neig);
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    n_MC = 10^5;
    s_Q0_MC = readmatrix(fullfile(repath,'s_MC.txt'));
    Q0s_MC = readmatrix(fullfile(repath,'Q0s_MC.txt'));
    wQ0s_MC = readmatrix(fullfile(repath,'wQ0s_MC.txt'));
    
    s_Q1_MC = readmatrix(fullfile(repath,'s_mu_IS.txt'));
    Q1s_MC = readmatrix(fullfile(repath,'mu1_IS.txt'));
    wQ1s_MC = readmatrix(fullfile(repath,'wmu_IS.txt'));
    
    % c0 = c*c1
    c = 11;
    
    %%
    n0_CV = fix(n_MC*c/(c+1));
    n1_CV = n0_CV;
    
    wQ0s_CV = wQ0s_MC(1:n0_CV);
    Q1s_CV = Q1(s_MC(1:n1_CV,:)')<0;
    w1_CV = w_MC(1:n1_CV);
    wQ1s_CV = w1_CV.*Q1s_CV;
    
    v0_MC = var(wQ0s_MC)/n_MC;
    cov_Q0Q1_CV = cov(wQ0s_CV,wQ1s_CV);
    cov_MC_Q0Q1_CV = cov_Q0Q1_CV/n0_CV;
    v_MC_Q0_CV = cov_MC_Q0Q1_CV(1,1);
    v_MC_Q1_CV = cov_MC_Q0Q1_CV(2,2);
    v_CV = v_MC_Q0_CV + a.^2*v_MC_Q1_CV + 2*a*cov_MC_Q0Q1_CV(1,2);
    
    %%
    r = 2.5;
    n0_IS = fix(n_MC*c/(c+r+1));
    n1_IS = n0_IS;
    m1_IS = fix(r*n0_IS);
    r1 = (n1_IS+m1_IS)/n0_IS;
    
    Q0s_IS = Q0s_MC(1:n0_IS);
    w0_IS = w_MC(1:n0_IS);
    wQ0s_IS = w0_IS.*Q0s_IS;
    Q1s_IS = Q1(s_MC(1:n1_IS))<0;
    w1_IS = w_MC(1:n1_IS);
    wQ1s_IS = w1_IS.*Q1s_IS;
    
    s_mu_IS = random(gm, m1_IS);
    mu1_IS = Q1(s_mu_IS')<0;
    w1_mu_IS = mvnpdf(s_mu_IS,mu,stdd)./qce(s_mu_IS);
    wmu_IS = w1_mu_IS.*mu1_IS;
    
    %%
    cov_Q0Q1_IS = cov(wQ0s_IS,wQ1s_IS);
    cov_MC_Q0Q1_IS = cov_Q0Q1_IS/n0_IS;
    v_MC_Q0_IS = cov_MC_Q0Q1_IS(1,1);
    v_MC_Q1_IS = cov_MC_Q0Q1_IS(2,2);
    v_IS = v_MC_Q0_IS + a.^2*(v_MC_Q1_IS*((r1-1)/r1)) + 2*a*(cov_MC_Q0Q1_IS(1,2)*((r1-1)/r1));
    
    m0s_IS = mean(wQ0s_IS);
    m1s_IS = mean(wQ1s_IS);
    mu1s_IS = (sum(wQ1s_IS) + sum(wmu_IS))/(n1_IS + m1_IS);
    m_IS = m0s_IS + a*(m1s_IS-mu1s_IS);
    
end