function cvis_ce()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2.8;
    l0 = mu+w0*std;
    l1 = mu+w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;

    v = 1.352806663625048e-09;
    prob0 = 0.001349898031630;
    prob1 = 1-normcdf(l1);
    
    a = linspace(-1.5,0.5,33);
    ansamples = 10000;
    wQ0s(1:ansamples) = 0;
    wQ1s(1:ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    % CE method
    N      = 3000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 10000;
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    qce = @(x) pdf(gm,x);

    for j = 1:ansamples
        samples = random(gm,nsamples);

        Q0s = Q0(samples(:))<0;
        Q1s = Q1(samples(:))<0;
        w = mvnpdf(samples,mu,std)./qce(samples);

        wQ0s(j) = mean(w.*Q0s);
        wQ1s(j) = mean(w.*Q1s);
    end
    
    m0 = mean(wQ0s);
    m1 = mean(wQ1s);
    me = m0+a*(m1-prob1);
      
    prob0
    prob1
    
    display('prob0./m0')
    prob0./m0
    display('prob0./me')
    prob0./me'
    
    covar = cov(wQ0s,wQ1s);
    v0 = covar(1,1);
    v1 = covar(2,2);
    ve = v0+a.^2*v1+2*a*covar(1,2);
      
    display('v./v0')
    v/v0
    display('v/ve')
    v./ve'
    display('ve./v0')
    ve'./v0
    
    figure(1)
    hold on
    plot(a,ve,'-o',a,v0*ones(1,length(a)),'--*')
    legend('v_e','v_0')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,ve./v0,'-o',a,ones(1,length(a)),'--*')
    legend('v_e/v_0')
    xlabel('alpha')
    hold off
    
    q0 = @(x) ((Q0(x)<0).*mvnpdf(x,mu,std))/prob0;
    q1 = @(x) ((Q1(x)<0).*mvnpdf(x,mu,std))/prob1;
    x = linspace(2,5.5,1000)';
    figure(3)
    hold on
    plot(x,q0(x),'-',x,q1(x),'-',x,qce(x),'--')
    l = legend('$q_0$','$q_1$','$\hat{q}$');
    set(l,'interpreter','latex')
    xlabel('z')
    hold off
    
end