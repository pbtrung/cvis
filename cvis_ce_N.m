function cvis_ce_N()
    
    dfile = '~/results.txt';
    if exist(dfile, 'file') 
        delete(dfile); 
    end
    diary(dfile);
    
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
    q1 = @(x) ((Q1(x)<0).*mvnpdf(x,mu,std))/prob1;
    
    a = linspace(-1.5,0.5,33);
    ansamples = 10000;
    wQ0s(1:ansamples) = 0;
    wQ1s(1:ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);
    
    n = 500:500:3000;
    kldiv(1:length(n)) = 0;
    min_ve_v0(1:length(n)) = 0;
    prob1_m1(1:length(n)) = 0;
    astar(1:length(n)) = 0;
    umin = 0;
    umax = 7;
    nkl = 10000000;
    nsamples = 100000;
    
    for r = 1:length(n)
        % CE method
        N      = n(r);    % total number of samples for each level
        p      = 0.1;     % quantile value to select samples for parameter update
        k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)

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
        
        astar(r) = -covar(1,2)/v1
        min_ve_v0(r) = (v0+astar(r)^2*v1+2*astar(r)*covar(1,2))/v0
        
        prob1_m1(r) = prob1./m1
        kldiv(r) = kl(nkl,umin,umax,q1,qce)
                    
    end
    
    figure(1)
    hold on
    plot(n,min_ve_v0,'-o')
    legend('min(v_e/v_0)')
    xlabel('N')
    hold off
    
    figure(2)
    hold on
    plot(n,kldiv,'--*')
    legend('KL Divergence')
    xlabel('N')
    hold off

end

function kldiv = kl(nkl,umin,umax,q1,qce)
    u = umin+(umax-umin)*rand(nkl,1);
    sample_value = q1(u);
    max_value = max(sample_value);
    accepted = rand(nkl,1)<(sample_value/max_value) & qce(u) ~= 0;
    s = u(accepted,:);
    if length(s) < 10^5
        error('kldiv: length(s) < 10^5');
    end
    kldiv = mean(log(q1(s(:)))-log(qce(s(:))));
end