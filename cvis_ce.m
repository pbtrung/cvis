function cvis_ce()
    
    format long;
    rng('default');
    
    lx = 60;
    ly = 20;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    mu = 0;
    sigma = 1;
    
%     outer = 1000;
%     inner = 1000;
%     [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,mu,sigma)

    EQ = 0.001251;
    EQ1 = 0.03141;
    VQ = 1.219218218218235e-06;
    l = 3.705226467613969e+02;
    l1 = 2.246457497067920e+02;
    
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    KE1 = ElmStiffnessMatrix(lx,ly,nelx1,nely1);
    dof = 2*(nely+1)*nelx+2;
    dof1 = 2*(nely1+1)*nelx1+2;
    
    Qx = @(x) l-FE3(nelx,nely,dof,KE,x);
    Q1y = @(y) l1-FE3(nelx1,nely1,dof1,KE1,y);
    select = @(v,idx) v(idx,:);
    Q = @(x) select(Qx(x),dof)';
    Q1 = @(y) select(Q1y(y),dof1)';
    
    alpha = linspace(-4,4,33);
    ansamples = 1000;
    e1(1:length(alpha),ansamples) = 0;
    e2(1:length(alpha),ansamples) = 0;
    wQs(1:length(alpha),ansamples) = 0;
    wQ1s(1:length(alpha),ansamples) = 0;
    
    % definition of the random variables
    d      = 1;
    pi_pdf = repmat(ERADist('standardnormal','PAR'),d,1);

    % CE method
    N      = 5000;    % total number of samples for each level
    p      = 0.1;     % quantile value to select samples for parameter update
    k_init = 3;       % initial number of distributions in the Mixture models (GM/vMFNM)
    nsamples = 1000; 
    
    % limit state function
    g = @(x) Q1(x);
    [~,~,~,~,~,~,~,mu_hat,Si_hat,Pi_hat] = CEIS_GM(N,p,g,pi_pdf,k_init);
    gm = gmdistribution(mu_hat,Si_hat,Pi_hat);
    
    tic
    for i = 1:length(alpha)
        parfor j = 1:ansamples
            fprintf('alpha: %f, iter: %d\n',alpha(i),j);
            samples = random(gm,nsamples);
            q = pdf(gm,samples);

            Qs = Q(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = mvnpdf(samples,mu,sigma)./q;

            wQs(i,j) = mean(w.*Qs);
            wQ1s(i,j) = mean(w.*Q1s);
        end
        fprintf('alpha: %f, time: %f\n',alpha(i),toc);
    end
    fprintf('finish: %f\n',toc);
    
    e2 = wQs;
    for i = 1:length(alpha)
        e1(i,:) = wQs(i,:)+alpha(i)*(wQ1s(i,:)-EQ1);
    end
    
    EQ
    EQ1
    
    m1 = mean(e1,2);
    m2 = mean(e2,2);
    display('EQ./m1')
    EQ./m1
    display('EQ./m2')
    EQ./m2

    v2 = var(e2,0,2);
    v1(1:length(alpha)) = 0;
    for i = 1:length(alpha)
        covar = cov(wQs(i,:),wQ1s(i,:));
        v1(i) = v2(i)+alpha(i)^2*var(wQ1s(i,:))+2*alpha(i)*covar(1,2);
        corrcoef(wQs(i,:),wQ1s(i,:))
    end
    
    v1 = v1'
    v2
    display('VQ./v1')
    VQ./v1
    display('VQ./v2')
    VQ./v2
    display('v1./v2')
    v1./v2
    
    figure(1)
    hold on
    plot(alpha,v1,'-o',alpha,v2,'--*')
    legend('v_0','v_1')
    xlabel('alpha')
    hold off
    
    figure(2)
    hold on
    plot(alpha,log(v1),'-o',alpha,log(v2),'--*')
    legend('ln(v_0)','ln(v_1)')
    xlabel('alpha')
    hold off
    
    figure(3)
    hold on
    plot(alpha,v1./v2,'-o',alpha,ones(1,length(alpha)),'--')
    legend('v_0/v_1')
    xlabel('alpha')
    hold off
    
    astar(1:length(alpha)) = 0;
    parfor i = 1:length(alpha)
        cov01 = cov(wQs(i,:),wQ1s(i,:));
        astar(i) = -cov01(1,2)/v2(i);
    end
    astar
    
end

function [EQ,EQ1,VQ] = EQ_EQ1(lx,ly,nelx,nely,nelx1,nely1,outer,inner,mu,sigma)

    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely+1)*nelx+2;
    KE = ElmStiffnessMatrix(lx,ly,nelx,nely);
    tic
    for i = 1:outer
        forces = mvnrnd(mu,sigma,inner);
        parfor j = 1:inner
            U = FE2(nelx,nely,dof,KE,forces(j));
            Uy(i,j) = U(dof);
        end
    end
    toc
    m = max(Uy(:));
    my = mean(mean(Uy,2));
    l = 0.65*m
    EQ = mean(mean(l-Uy<0,2));
    VQ = var(mean(l-Uy<0,2));
    
    Uy(1:outer,1:inner) = 0;
    dof = 2*(nely1+1)*nelx1+2;
    KE = ElmStiffnessMatrix(lx,ly,nelx1,nely1);
    tic
    for i = 1:outer
        forces = mvnrnd(mu,sigma,inner);
        parfor j = 1:inner
            U = FE2(nelx1,nely1,dof,KE,forces(j));
            Uy(i,j) = U(dof);
        end
    end
    toc
    m = max(Uy(:));
    my = mean(mean(Uy,2));
    l = 0.4*m
    EQ1 = mean(mean(l-Uy<0,2));
end

function KE = ElmStiffnessMatrix(lx,ly,nelx,nely)
    E = 1.;
    nu = 0.3;
    h = 1;
    a = lx/nelx;
    b = ly/nely;
    c = b/a;
    
    ka = (2/c)*(1-nu);
    kb = (3/2)*(1+nu);
    kc = (3/2)*(1-3*nu);
    kd = (2*c)*(1-nu);
    KE = (E*h)/(12*(1-nu^2))*[ ...
        4*c+ka kb -4*c+ka/2 -kc -2*c-ka/2 -kb 2*c-ka kc
        kb 4/c+kd kc 2/c-kd -kb -2/c-kd/2 -kc -4/c+kd/2
        -4*c+ka/2 kc 4*c+ka -kb 2*c-ka -kc -2*c-ka/2 kb
        -kc 2/c-kd -kb 4/c+kd kc -4/c+kd/2 kb -2/c-kd/2
        -2*c-ka/2 -kb 2*c-ka kc 4*c+ka kb -4*c+ka/2 -kc
        -kb -2/c-kd/2 -kc -4/c+kd/2 kb 4/c+kd kc 2/c-kd
        2*c-ka -kc -2*c-ka/2 kb -4*c+ka/2 kc 4*c+ka -kb
        kc -4/c+kd/2 kb -2/c-kd/2 -kc 2/c-kd -kb 4/c+kd
        ];
end