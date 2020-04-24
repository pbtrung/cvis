function fig(p0,p1)
    
    mu = 0;
    std = 1;
    w1 = 2.8;
    l1 = mu+w1*std;

    prob0 = 0.001349898031630;
    prob1 = 1-normcdf(l1);
    
    a = linspace(-1.5,0.5,33);
    
    wQ0s = writematrix(p0);
    wQ1s = writematrix(p1);
      
    covar = cov(wQ0s,wQ1s);
    v0 = covar(1,1);
    v1 = covar(2,2);
    ve = v0+a.^2*v1+2*a*covar(1,2);
    
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