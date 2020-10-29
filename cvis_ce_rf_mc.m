function cvis_ce_rf_mc()
    
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
    s_MC = random(gm, n_MC);
    writematrix(s_MC, fullfile(repath,'s_MC.txt'));
    Q0s_MC = Q0(s_MC')<0;
    writematrix(Q0s_MC, fullfile(repath,'Q0s_MC.txt'));
    w0_MC = mvnpdf(s_MC,mu,stdd)./qce(s_MC);
    writematrix(w0_MC, fullfile(repath,'w0_MC.txt'));
    wQ0s_MC = w0_MC.*Q0s_MC;
    writematrix(wQ0s_MC, fullfile(repath,'wQ0s_MC.txt'));
    Q1s_MC = Q1(s_MC')<0;
    writematrix(Q1s_MC, fullfile(repath,'Q1s_MC.txt'));
    w1_MC = mvnpdf(s_MC,mu,stdd)./qce(s_MC);
    writematrix(w1_MC, fullfile(repath,'w1_MC.txt'));
    wQ1s_MC = w1_MC.*Q1s_MC;
    writematrix(wQ1s_MC, fullfile(repath,'wQ1s_MC.txt'));
    
    disp(mean(wQ0s_MC));
    disp(mean(wQ1s_MC));
    
    m1_IS = 10^6;
    s_IS = random(gm, m1_IS);
    writematrix(s_IS, fullfile(repath,'s_IS.txt'));
    Q1s_IS = Q1(s_IS')<0;
    writematrix(Q1s_IS, fullfile(repath,'Q1s_IS.txt'));
    w1_IS = mvnpdf(s_IS,mu,stdd)./qce(s_IS);
    writematrix(w1_IS, fullfile(repath,'w1_IS.txt'));
    wQ1s_IS = w1_IS.*Q1s_IS;
    writematrix(wQ1s_IS, fullfile(repath,'wQ1s_IS.txt'));
    
end