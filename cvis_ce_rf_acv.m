function cvis_ce_rf_acv()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 60;
    Ly = 20;
    lower = 1;
    upper = 2;
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx,nely,2);

    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);

    l0 = 118.923;
    l1 = 108.510;
%     EQ0 = 0.001173;
%     VQ0 = 1.186257257257260e-06;
%     EQ1 = 0.022782;
    
    dof0 = 2*(nely+1)*nelx+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    force = 1;
    
    Q0 = @(x) l0-FEQ0(lx,ly,nelx,nely,dof0,force,x,P,L,lower,upper)';
    Q1 = @(x) l1-FEQ1(lx,ly,nelx1,nely1,dof1,force,x,PsiE1,L,lower,upper)';
        
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
    d      = neig;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    ns_mfis = 200;
    ns_Q0_cv = 190;
    ns_Q1_cv = 100;
    
    ns_Q0_acv1 = 180;
    ns_Q1_acv1 = 100;
    ns_mu_acv1 = 100;
    
    ns_Q0_acv2 = 180;
    ns_Q1_acv2 = 100;
    ns_mu_acv2 = 100;
     
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    mu = zeros(1,neig);
    stdd = eye(neig);
    for j = 1:K
        fprintf('K: %d\n',j);
        tic
        ss = random(gm,ns_mfis);
        Q0s = Q0(ss')<0;
        w = mvnpdf(ss,mu,stdd)./qce(ss);
        wQ0s(j) = mean(w.*Q0s);
        toc
        
        wQ0s_cv(j) = mean(w(1:ns_Q0_cv).*Q0s(1:ns_Q0_cv));
        Q1s = Q1(ss(1:ns_Q1_cv,:)')<0;
        wQ1s_cv(j) = mean(w(1:ns_Q1_cv).*Q1s(1:ns_Q1_cv));
        
        wQ0s_acv1(j) = mean(w(1:ns_Q0_acv1).*Q0s(1:ns_Q0_acv1));
        Q1s_acv1 = Q1(ss(1:ns_Q1_acv1,:)')<0;
        wQ1s_acv1(j) = mean(w(1:ns_Q1_acv1).*Q1s_acv1(1:ns_Q1_acv1));
        ss_mu_acv1 = random(gm,ns_mu_acv1);
        Q1s_mu_acv1 = Q1(ss_mu_acv1')<0;
        w_mu_acv1 = mvnpdf(ss_mu_acv1,mu,stdd)./qce(ss_mu_acv1);
        mu_acv1(j) = mean(w_mu_acv1(1:ns_mu_acv1).*Q1s_mu_acv1(1:ns_mu_acv1));
        
        wQ0s_acv2(j) = mean(w(1:ns_Q0_acv2).*Q0s(1:ns_Q0_acv2));
        if ns_Q1_acv1 == ns_Q1_acv2
            Q1s_acv2 = Q1s_acv1;
        else
            error('ns_Q1_acv1 != ns_Q1_acv2')
        end
        wQ1s_acv2(j) = mean(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2));
        ss_mu_acv2 = random(gm,ns_mu_acv2);
        Q1s_mu_acv2 = Q1(ss_mu_acv2')<0;
        w_mu_acv2 = mvnpdf(ss_mu_acv2,mu,stdd)./qce(ss_mu_acv2);
        mu_acv2(j) = (sum(w(1:ns_Q1_acv2).*Q1s_acv2(1:ns_Q1_acv2))+sum(w_mu_acv2(1:ns_mu_acv2).*Q1s_mu_acv2(1:ns_mu_acv2)))/(ns_Q1_acv2+ns_mu_acv2);
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
    
%     m0 = mean(wQ0s);
%     m1 = mean(wQ1s);
%     m1_acv1 = mean(wQ1s_acv1);
%     m1_acv2 = mean(wQ1s_acv2);
%     me = m0+a*(m1-EQ1);
%     me_acv1 = m0+a*(m1-m1_acv1);
%     me_acv2 = m0+a*(m1-m1_acv2);
%       
%     EQ0
%     EQ1
%     
%     display('EQ0./m0')
%     EQ0./m0
%     display('EQ0./me')
%     EQ0./me'
%     display('EQ0./me_acv1')
%     EQ0./me_acv1'
%     display('EQ0./me_acv2')
%     EQ0./me_acv2'
    
%     covar = cov(wQ0s,wQ1s);
%     covar_acv1 = cov(wQ1s,wQ1s_acv1);
%     covar_acv2 = cov(wQ1s,wQ1s_acv2);
%     covar_0acv2 = cov(wQ0s,wQ1s_acv2);
%     v0 = covar(1,1);
%     v1 = covar(2,2);
%     v_acv1 = covar_acv1(2,2);
%     v_acv2 = covar_acv2(2,2);
%     ve = v0 + a.^2*v1 + 2*a*covar(1,2);
%     ve_acv1 = v0 + a.^2*(v1+v_acv1) + 2*a*covar(1,2);
%     ve_acv2 = v0 + a.^2*(v1+v_acv2-2*covar_acv2(1,2)) + 2*a*(covar(1,2)-covar_0acv2(1,2));
      
    display('ve./v0')
    ve'./v0
    display('ve_acv1./v0')
    ve_acv1'./v0
    display('ve_acv2./v0')
    ve_acv2'./v0
    
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