function cvis_rs_tols()
    
    format long;
    diary '~/results.txt';
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
    
    nsamples = 10000;
    umin = 0;
    umax = 7;
    
    tol = [0 0.1 0.2 0.4 0.6 0.8 1 1.2];
    min_ve_v0(1:length(tol)) = 0;
    
    for r = 1:length(tol)
        
        fprintf('tol: %f\n',tol(r));
        prob2 = 1-normcdf(l1-tol(r));
        Q2 = @(x) l1-tol(r)-x;
        q = @(x) ((Q2(x)<0).*normpdf(x))/prob2;

        for j = 1:ansamples
            % rejection sampling
            u = umin+(umax-umin)*rand(nsamples,1);
            sample_value = q(u);
            max_value = max(sample_value);
            accepted = rand(nsamples,1)<(sample_value/max_value);
            samples = u(accepted,:);
            if length(samples) < 200
                error('length(samples) < 200');
            end

            Q0s = Q0(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = normpdf(samples)./q(samples);

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
        
        min_ve_v0(r) = min(ve./v0)
    
    end
    
    figure(1)
    hold on
    plot(tol,min_ve_v0,'-o')
    legend('min(v_e/v_0)')
    xlabel('tol')
    hold off
    
end