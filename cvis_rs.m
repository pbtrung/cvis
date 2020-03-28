function cvis_rs()
    
    format long;
    % rng('default');
    
    mu = 0;
    std = 1;
    l = mu+3*std;
    l1 = mu+2*std;
    
    Q = @(x) l-x;
    Q1 = @(y) l1-y;

    v = 1.352806663625048e-09;
    prob = 0.001349898031630;
    prob1 = 0.022750131948179;
    
    a = linspace(-4,4,10);
    ansamples = 1000;
    e1(1:length(a),ansamples) = 0;
    e2(1:length(a),ansamples) = 0;
    wQs(1:length(a),ansamples) = 0;
    wQ1s(1:length(a),ansamples) = 0;
    vals(1:length(a),ansamples) = 0;
    
    nsamples = 100000;
    umin = 0;
    umax = 7;
    q = @(x) ((Q1(x)<0).*normpdf(x))/prob1;
    
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
            w = normpdf(samples)./q(samples);

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
            
            vals(i,j) = wQs(i,j)+a(i)*(wQ1s(i,j)-prob1);
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
    legend('v_0','v_1')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(a,log(v1),'-o',a,log(v2),'--*')
    legend('log(v_0)','log(v_1)')
    xlabel('alpha')
    hold off
    
    figure(3)
    hold on
    plot(a,var(vals,0,2))
    xlabel('alpha')
    hold off
    
end