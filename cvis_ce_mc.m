function cvis_ce_mc()
    
    format long;
    rng('default');
    
    %% materials
    E = 10000;
    poisson = 0.30;
    kapa = 5/6;
    L = 1;
    
    %%
    nelx0 = 30;
    nely0 = 30;
    nelx1 = 10;
    nely1 = 10;
    if nelx0 ~= nely0 || mod(nelx0,2) ~= 0 || nelx1 ~= nely1 || mod(nelx1,2) ~= 0
        error("Check inputs")
    end
    
    l0 = 0.013165;
    l1 = 0.011284;
    Q0 = @(x) l0-FE_plate(nelx0,nely0,E,poisson,kapa,L,x);
    Q1 = @(x) l1-FE_plate(nelx1,nely1,E,poisson,kapa,L,x);
    
    %% definition of the random variables
    d = 8;
    lower = 1;
    upper = 2;
    pi_pdf = repmat(ERADist('uniform','PAR',[lower upper]),d,1);
    if any(strcmp('Marginals',fieldnames(pi_pdf))) == 1
        % use Nataf transform (dependence)
        u2x = @(u) pi_pdf.U2X(u);
    else
        % use distribution information for the transformation (independence)
        u2x = @(u) pi_pdf(1).icdf(normcdf(u));
    end
    
    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 5;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    %%
    n_MC = 5*10^5;
    
    s_MC = random(gm, n_MC);
    writematrix(s_MC, fullfile(repath,'s_MC.txt'));
    
    mu = zeros(1,d);
    stdd = eye(d);
    
    w0_MC = mvnpdf(s_MC,mu,stdd)./qce(s_MC);
    writematrix(w0_MC, fullfile(repath,'w0_MC.txt'));
    w1_MC = mvnpdf(s_MC,mu,stdd)./qce(s_MC);
    writematrix(w1_MC, fullfile(repath,'w1_MC.txt'));
    
    s_MC = u2x(s_MC);
    
    Q0s_MC = Q0(s_MC)<0;
    writematrix(Q0s_MC, fullfile(repath,'Q0s_MC.txt'));
    wQ0s_MC = w0_MC.*Q0s_MC;
    writematrix(wQ0s_MC, fullfile(repath,'wQ0s_MC.txt'));
    
    Q1s_MC = Q1(s_MC)<0;
    writematrix(Q1s_MC, fullfile(repath,'Q1s_MC.txt'));
    wQ1s_MC = w1_MC.*Q1s_MC;
    writematrix(wQ1s_MC, fullfile(repath,'wQ1s_MC.txt'));
    
    %%
    m1_IS = 2*10^6;
    
    s_IS = random(gm, m1_IS);
    writematrix(s_IS, fullfile(repath,'s_IS.txt'));
    w1_IS = mvnpdf(s_IS,mu,stdd)./qce(s_IS);
    writematrix(w1_IS, fullfile(repath,'w1_IS.txt'));
    
    s_IS = u2x(s_IS);
    
    Q1s_IS = Q1(s_IS)<0;
    writematrix(Q1s_IS, fullfile(repath,'Q1s_IS.txt'));
    wQ1s_IS = w1_IS.*Q1s_IS;
    writematrix(wQ1s_IS, fullfile(repath,'wQ1s_IS.txt'));
    
    disp(mean(wQ0s_MC));
    disp(mean(wQ1s_MC));
    disp(mean(wQ1s_IS));
    
end