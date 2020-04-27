function cvis_rs_tols()
    
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
    q1 = @(x) ((Q1(x)<0).*normpdf(x))/prob1;
    
    a = linspace(-1.5,0.5,33);
    ansamples = 10000;
    wQ0s(1:ansamples) = 0;
    wQ1s(1:ansamples) = 0;
    
    nsamples = 10000;
    umin = 0;
    umax = 7;
    
    tol = [0 0.1 0.2 0.4 0.6 0.8 1 1.2];
    min_ve_v0(1:length(tol)) = 0;
    astar(1:length(tol)) = 0;
    kldiv(1:length(tol)) = 0;
    covar(1:length(tol)) = 0;
    nkl = 10000000;
    
    for r = 1:length(tol)
        
        fprintf('tol: %f\n',tol(r));
        prob2 = 1-normcdf(l1-tol(r));
        Q2 = @(x) l1-tol(r)-x;
        q2 = @(x) ((Q2(x)<0).*normpdf(x))/prob2;

        for j = 1:ansamples
            % rejection sampling
            u = umin+(umax-umin)*rand(nsamples,1);
            sample_value = q2(u);
            max_value = max(sample_value);
            accepted = rand(nsamples,1)<(sample_value/max_value);
            samples = u(accepted,:);
            if length(samples) < 200
                error('length(samples) < 200');
            end

            Q0s = Q0(samples(:))<0;
            Q1s = Q1(samples(:))<0;
            w = normpdf(samples)./q2(samples);

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

        co = cov(wQ0s,wQ1s);
        v0 = co(1,1);
        v1 = co(2,2);
        covar(r) = co(1,2)
        ve = v0+a.^2*v1+2*a*covar(r);

        display('v./v0')
        v/v0
        display('v/ve')
        v./ve'
        display('ve./v0')
        ve'./v0
        
        astar(r) = -covar(r)/v1;
        astar(2:end)
        min_ve_v0(r) = (v0+astar(r)^2*v1+2*astar(r)*covar(r))/v0
        
        kldiv(r) = kl(nkl,umin,umax,q1,q2)
    
    end
    
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 3;
    height = 3;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(l1-tol,min_ve_v0,'-o',l1-tol,kldiv,'-o',l1-tol,covar*10^7,'-o')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    legend({'$\min(v_e/v_0)$','KL Divergence','Covariance$\times 10^7$'},...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'Location','NorthEast')
    ylabel('',...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontSize',fontsize,...
        'FontName','Times')
    xlabel({'Threshold'},...
        'FontUnits','points',...
        'interpreter','latex',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times')
    fn(1) = "tol_min.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', sprintf('%s', fn(k)))
    end
    
end

function kldiv = kl(nkl,umin,umax,q1,q2)
    u = umin+(umax-umin)*rand(nkl,1);
    sample_value = q1(u);
    max_value = max(sample_value);
    accepted = rand(nkl,1)<(sample_value/max_value);
    s = u(accepted,:);
    if length(s) < 10^5
        error('kldiv: length(s) < 10^5');
    end
    kldiv = mean(log(q1(s(:)))-log(q2(s(:))));
end