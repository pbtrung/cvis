function cvis_ce_acv()    

    format long;
    rng('default');

    nelx = 20;
    nely = 20;
    nelx1 = 10;
    nely1 = 10;
    if nelx ~= nely || mod(nelx,2) ~= 0 || nelx1 ~= nely1 || mod(nelx1,2) ~= 0
        error("Check inputs")
    end
        
    l0 = 0.013165;
    l1 = 0.011284;
%     EQ0 = 
%     VQ0 = 
%     EQ1 = 0.011165;
    
    % x(nsamples,8)
    Q0 = @(x) l0-FE_plate(nelx,nely,x);
    Q1 = @(x) l1-FE_plate(nelx1,nely1,x);
    
    a = linspace(-1.5,0.5,33);
    K = 1000;
    
    wQ0s(1:K) = 0;
    wQ0s_cv(1:K) = 0;
    wQ1s_cv(1:K) = 0;
    wQ0s_acv1(1:K) = 0;
    wQ1s_acv1(1:K) = 0;
    mu_acv1(1:K) = 0;
    wQ0s_acv2(1:K) = 0;
    wQ1s_acv2(1:K) = 0;
    mu_acv2(1:K) = 0;
    
    % definition of the random variables
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
    
    ns_mfis = 200;
    ns_Q0_cv = 165;
    ns_Q1_cv = 165;
    
    ns_Q0_acv1 = 112;
    ns_Q1_acv1 = 112;
    ns_mu_acv1 = 300;
    
    ns_Q0_acv2 = 147;
    ns_Q1_acv2 = 147;
    ns_mu_acv2 = 100;
    
    % limit state function
    g = @(x) Q1(x');
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    for j = 1:K
        fprintf('K: %d\n',j);
        tic
        ss = random(gm,ns_mfis);
        w = mvnpdf(ss,zeros(1,d),eye(d))./qce(ss);
        ss = u2x(ss);
        Q0s = Q0(ss)<0;
        wQ0s(j) = mean(w.*Q0s);
        
        wQ0s_cv(j) = mean(w(1:ns_Q0_cv).*Q0s(1:ns_Q0_cv));
        Q1s = Q1(ss(1:ns_Q1_cv,:))<0;
        wQ1s_cv(j) = mean(w(1:ns_Q1_cv).*Q1s(1:ns_Q1_cv));
        
        wQ0s_acv1(j) = mean(w(1:ns_Q0_acv1).*Q0s(1:ns_Q0_acv1));
        Q1s_acv1 = Q1(ss(1:ns_Q1_acv1,:))<0;
        wQ1s_acv1(j) = mean(w(1:ns_Q1_acv1).*Q1s_acv1(1:ns_Q1_acv1));
        ss_mu_acv1 = random(gm,ns_mu_acv1);
        w_mu_acv1 = mvnpdf(ss_mu_acv1,zeros(1,d),eye(d))./qce(ss_mu_acv1);
        ss_mu_acv1 = u2x(ss_mu_acv1);
        Q1s_mu_acv1 = Q1(ss_mu_acv1)<0;
        mu_acv1(j) = mean(w_mu_acv1(1:ns_mu_acv1).*Q1s_mu_acv1(1:ns_mu_acv1));
        
        wQ0s_acv2(j) = mean(w(1:ns_Q0_acv2).*Q0s(1:ns_Q0_acv2));
        Q1s_acv2 = Q1(ss(1:ns_Q1_acv2,:))<0;
        wQ1s_acv2(j) = mean(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2));
        ss_mu_acv2 = random(gm,ns_mu_acv2);
        w_mu_acv2 = mvnpdf(ss_mu_acv2,zeros(1,d),eye(d))./qce(ss_mu_acv2);
        ss_mu_acv2 = u2x(ss_mu_acv2);
        Q1s_mu_acv2 = Q1(ss_mu_acv2)<0;
        mu_acv2(j) = (sum(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2))+sum(w_mu_acv2(1:ns_mu_acv2).*Q1s_mu_acv2(1:ns_mu_acv2)))/(ns_Q1_acv2+ns_mu_acv2);
        toc
    end
    writematrix(wQ0s,'wQ0s.txt');
    
    writematrix(wQ0s_cv,'wQ0s_cv.txt');
    writematrix(wQ1s_cv,'wQ1s_cv.txt');
    
    writematrix(wQ0s_acv1,'wQ0s_acv1.txt');
    writematrix(wQ1s_acv1,'wQ1s_acv1.txt');
    writematrix(mu_acv1,'mu_acv1.txt');
    
    writematrix(wQ0s_acv2,'wQ0s_acv2.txt');
    writematrix(wQ1s_acv2,'wQ1s_acv2.txt');
    writematrix(mu_acv2,'mu_acv2.txt');
    
    v0 = var(wQ0s);
    cov_cv = cov(wQ0s_cv,wQ1s_cv);
    v0_cv = cov_cv(1,1);
    v1_cv = cov_cv(2,2);
    ve = v0_cv + a.^2*v1_cv + 2*a*cov_cv(1,2);
    
    cov_acv1 = cov(wQ0s_acv1,wQ1s_acv1);
    v0_acv1 = cov_acv1(1,1);
    v1_acv1 = cov_acv1(2,2);
    cov_Q1_mu_acv1 = cov(wQ1s_acv1,mu_acv1);
    v_mu_acv1 = cov_Q1_mu_acv1(2,2);
    ve_acv1 = v0_acv1 + a.^2*(v1_acv1+v_mu_acv1) + 2*a*cov_acv1(1,2);
    
    cov_acv2 = cov(wQ0s_acv2,wQ1s_acv2);
    v0_acv2 = cov_acv2(1,1);
    v1_acv2 = cov_acv2(2,2);
    cov_Q1_mu_acv2 = cov(wQ1s_acv2,mu_acv2);
    v_mu_acv2 = cov_Q1_mu_acv2(2,2);
    cov_0acv2 = cov(wQ0s_acv2,mu_acv2);
    ve_acv2 = v0_acv2 + a.^2*(v1_acv2+v_mu_acv2-2*cov_Q1_mu_acv2(1,2)) + 2*a*(cov_acv2(1,2)-cov_0acv2(1,2));
      
%     display('VQ0/v0')
%     VQ0/v0
%     display('VQ0./ve')
%     VQ0./ve'
    display('ve./v0')
    ve'./v0
    display('ve_acv1./v0')
    ve_acv1'./v0
    display('ve_acv2./v0')
    ve_acv2'./v0
    
    min(ve'./v0)
    min(ve_acv1'./v0)
    min(ve_acv2'./v0)
    
%     figure(1)
%     hold on
%     plot(a,ve,'-o',a,v0*ones(1,length(a)),'--*')
%     legend('v_e','v_0')
%     xlabel('alpha')
%     hold off
    
    figure(2)
    hold on
    plot(a,ve./v0,'-s',a,ve_acv1./v0,'-v',a,ve_acv2./v0,'-d',a,ones(1,length(a)),'-o')
    legend('v_e/v_0','v_1/v_0','v_2/v_0')
    xlabel('alpha')
    hold off
    
end