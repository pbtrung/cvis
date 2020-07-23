function cvis_ce_acv()    

    format long;
    rng('default');

    nelx = 30;
    nely = 30;
    nelx1 = 20;
    nely1 = 20;
    if nelx ~= nely || mod(nelx,2) ~= 0 || nelx1 ~= nely1 || mod(nelx1,2) ~= 0
        error("Check inputs")
    end
        
    l0 = 0.013165;
    l1 = 0.011284;
%     EQ0 = 
%     VQ0 = 
    EQ1 = 0.011165;
    
    % x(nsamples,8)
    Q0 = @(x) l0-FE_plate(nelx,nely,x);
    Q1 = @(x) l1-FE_plate(nelx1,nely1,x);
    
    a = linspace(-1.5,0.5,33);
    ansamples = 1000;
    wQ0s(1:ansamples) = 0;
    wQ1s(1:ansamples) = 0;
    wQ1s_acv1(1:ansamples) = 0;
    wQ1s_acv2(1:ansamples) = 0;
    
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
    nsamples = 100;
    acv1_nsamples = 100;
    acv2_nsamples = 100;
    
    % limit state function
    g = @(x) Q1(x');
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);
    
    for j = 1:ansamples
        fprintf('ansamples: %d\n',j);
        tic
        samples = random(gm,nsamples);
        
        w = mvnpdf(samples,zeros(1,d),eye(d))./qce(samples);
        samples = u2x(samples);
        Q0s = Q0(samples)<0;
        Q1s = Q1(samples)<0;

        wQ0s(j) = mean(w.*Q0s);
        wQ1s(j) = mean(w.*Q1s);
        toc
        
        samples_acv1 = random(gm,acv1_nsamples);
        w_acv1 = mvnpdf(samples_acv1,mu,stdd)./qce(samples_acv1);
        samples_acv1 = u2x(samples_acv1);
        Q1s_acv1 = Q1(samples_acv1')<0;
        wQ1s_acv1(j) = mean(w_acv1.*Q1s_acv1);
        
        samples_acv2 = random(gm,acv2_nsamples);
        w_acv2 = mvnpdf(samples_acv2,mu,stdd)./qce(samples_acv2);
        samples_acv2 = u2x(samples_acv2);
        Q1s_acv2 = Q1(samples_acv2')<0;
        wQ1s_acv2(j) = (sum(w.*Q1s)+sum(w_acv2.*Q1s_acv2))/(nsamples+acv2_nsamples);
    end
    writematrix(wQ0s,'wQ0s.txt');
    writematrix(wQ1s,'wQ1s.txt');
    writematrix(wQ1s_acv1,'wQ1s_acv1.txt');
    writematrix(wQ1s_acv2,'wQ1s_acv2.txt');
    
    m0 = mean(wQ0s);
    m1 = mean(wQ1s);
    m1_acv1 = mean(wQ1s_acv1);
    m1_acv2 = mean(wQ1s_acv2);
    me = m0+a*(m1-EQ1);
    me_acv1 = m0+a*(m1-m1_acv1);
    me_acv2 = m0+a*(m1-m1_acv2);
      
    EQ0
    EQ1
    
    display('EQ0./m0')
    EQ0./m0
    display('EQ0./me')
    EQ0./me'
    display('EQ0./me_acv1')
    EQ0./me_acv1'
    display('EQ0./me_acv2')
    EQ0./me_acv2'
    
    covar = cov(wQ0s,wQ1s);
    covar_acv1 = cov(wQ1s,wQ1s_acv1);
    covar_acv2 = cov(wQ1s,wQ1s_acv2);
    covar_0acv2 = cov(wQ0s,wQ1s_acv2);
    v0 = covar(1,1);
    v1 = covar(2,2);
    v_acv1 = covar_acv1(2,2);
    v_acv2 = covar_acv2(2,2);
    ve = v0 + a.^2*v1 + 2*a*covar(1,2);
    ve_acv1 = v0 + a.^2*(v1+v_acv1) + 2*a*covar(1,2);
    ve_acv2 = v0 + a.^2*(v1+v_acv2-2*covar_acv2(1,2)) + 2*a*(covar(1,2)-covar_0acv2(1,2));
      
    display('VQ0/v0')
    VQ0/v0
    display('VQ0./ve')
    VQ0./ve'
    display('ve./v0')
    ve./v0
    display('ve_acv1./v0')
    ve_acv1./v0
    display('ve_acv2./v0')
    ve_acv2./v0
    
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