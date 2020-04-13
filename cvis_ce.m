function cvis_ce()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2;
    l0 = mu+w0*std;
    l1 = mu+w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;

    v = 1.352806663625048e-09;
    prob0 = 0.001349898031630;
    prob1 = 1-normcdf(w1);
    
    a = linspace(-2,2,33);
    ansamples = 100000;
    wQ0s(ansamples) = 0;
    wQ1s(ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    n = [100 300 500 1000 2000 4000];
    kldiv(1:length(n)) = 0;
    minvev0(1:length(n)) = 0;
    
    for r = 1:length(n)
    % CE method
    N      = n(r);     % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 1000;
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init,mu,std);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);

    for j = 1:ansamples
        samples = random(gm,nsamples);
        q = pdf(gm,samples);

        Q0s = Q0(samples(:))<0;
        Q1s = Q1(samples(:))<0;
        w = mvnpdf(samples,mu,std)./q;

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
     
%     figure(1)
%     hold on
%     plot(a,v0*ones(1,length(a)),'-o',a,ve,'--*')
%     legend('v_0','v_e')
%     xlabel('alpha')
%     hold off
%     
%     figure(2)
%     hold on
%     plot(a,ve./v0,'-o',a,ones(1,length(a)),'--')
%     legend('v_e/v_0')
%     xlabel('alpha')
%     hold off
    
    astar = -covar(1,2)/covar(2,2);
    min_ve = v0+astar^2*v1+2*astar*covar(1,2);
    display('astar')
    astar
    display('covar')
    covar
    display('min(ve)/v0')
    minvev0(r) = min(ve)/v0
    display('min_ve/v0')
    min_ve/v0
    
    qce = @(x) pdf(gm,x);
    n0 = @(x) normpdf(x,mu,std);
    n1 = @(x) normpdf(x,mu,std);
    q0 = @(x) ((Q0(x)<0).*n0(x))./prob0;
    q1 = @(x) ((Q1(x)<0).*n1(x))./prob1;
    
    nrs = 10000000;
    umin = 0;
    umax = 7;
    u = umin+(umax-umin)*rand(nrs,1);
    sample_value = q1(u);
    max_value = max(sample_value);
    accepted = rand(nrs,1)<(sample_value/max_value);
    s = u(accepted,:);
    length(s)
    kldiv(r) = mean(log(q1(s(:)))-log(qce(s(:))))
    end
    plot(n,minvev0,'-o',n,kldiv,'--*')
    
%     qce = @(x) pdf(gm,x);
%     n0 = @(x) normpdf(x,mu,std);
%     n1 = @(x) normpdf(x,mu,std);
%     q0 = @(x) ((Q0(x)<0).*n0(x))./prob0;
%     q1 = @(x) ((Q1(x)<0).*n1(x))./prob1;
%     figure(3)
%     hold on
%     X = linspace(0,5,1000)';
%     plot(X,q0(X),'-',X,q1(X),'--',X,qce(X))
%     legend('q_0','q_1','q_{CE}')
%     hold off
    
end