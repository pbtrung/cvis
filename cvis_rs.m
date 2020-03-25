function cvis_rs()
    
    format long;
    % rng('default');
    
    mu = 0;
    std = 1;
    l = mu+3*std;
    l1 = mu+2*std;
    
    Q = @(x) l-x;
    Q1 = @(y) l1-y;
    
    vs = 1000;
    mc(1:vs) = 0;
    mc1(1:vs) = 0;
    parfor i = 1:vs
        s = mvnrnd(mu,std,1000000);
        mc(i) = mean(Q(s(:))<0);
        mc1(i) = mean(Q1(s(:))<0);
    end
    prob = mean(mc);
    v = var(mc);
    prob1 = mean(mc1);
    
    a = linspace(-1,0,3);
    ansamples = 10000;
    e1(1:length(a),ansamples) = 0;
    e2(1:length(a),ansamples) = 0;
    wQs(1:length(a),ansamples) = 0;
    wQ1s(1:length(a),ansamples) = 0;
    
    nsamples = 100000;
    umin = 0;
    umax = 6;
    q = @(x) (Q1(x)<0).*mvnpdf(x,mu,std)./prob1;
    
    for i = 1:length(a)
        parfor j = 1:ansamples
            % rejection sampling
            u = umin+(umax-umin)*rand(nsamples,1);
            sample_value = q(u);
            max_value = max(sample_value);
            accepted = rand(nsamples,1)<(sample_value/max_value);
            samples = u(accepted,:);
            
            Qs = Q(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = mvnpdf(samples,mu,std)./q(samples);

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
        end
    end
    
    e2 = wQs;
    for i = 1:length(a)
        e1(i,:) = wQs(i,:)+a(i)*(wQ1s(i,:)-prob1);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
    
    prob
    prob1
    
    m1 = mean(e1,2);
    m2 = mean(e2,2);
    prob./m1
    prob./m2
    
    v1 = var(e1,0,2)
    v2 = var(e2,0,2)
    v./v1
    v./v2
    v1./v2
    
    figure(1)
    hold on
    plot(a,v1,'-o',a,v2,'--*')
    legend('v1','v2')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,log(v1),'-o',a,log(v2),'--*')
    legend('log(v1)','log(v2)')
    xlabel('alpha')
    hold off
    
end