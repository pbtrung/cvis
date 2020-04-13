function test()
    
    format long;
    rng('default');
    
    mu = 0;
    std = 1;
    b = 3;
    b1 = 2;
    l = mu+b*std;
    l1 = mu+b1*std;
    
    Q = @(x) l-x;
    Q1 = @(y) l1-y;
    
%     vs = 100000;
%     mc(1:vs) = 0;
%     parfor i = 1:vs
%         s = mvnrnd(mu,std,1000000);
%         mc(i) = mean(Q(s(:))<0);
%     end
%     v = var(mc)
%     prob = 1-normcdf(b)
%     prob1 = 1-normcdf(b1)

    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);

    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init,mu,std);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    
    q = @(y) pdf(gm,y);
    nx = @(x) normpdf(x,mu,std);
    ny = @(y) normpdf(y,mu,std);
    probx = integral(nx,l,Inf);
    proby = integral(ny,l1,Inf);
    qx = @(x) ((Q(x)<0).*nx(x))./probx;
    qy = @(y) ((Q1(y)<0).*ny(y))./proby;
    
    figure(1)
    hold on
    X = linspace(0,5,10000)';
    ss = random(gm,5000);
    plot(X,qx(X),'-',X,qy(X),'--',X,q(X))
%     histogram(ss+1,100,'Normalization','pdf');
    legend('q-high','q-low','q-CE')
    hold off
    
end