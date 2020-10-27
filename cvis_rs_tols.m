function cvis_rs_tols(path)

    mu = 0;
    std = 1;
    w0 = 3;
    w1 = 2.8;
    l0 = mu + w0*std;
    l1 = mu + w1*std;
    
    Q0 = @(x) l0-x;
    Q1 = @(x) l1-x;
    
    % prob0 = 1-normcdf(l0);
    prob1 = 1-normcdf(l1);
    q1 = @(x) ((Q1(x)<0).*normpdf(x))/prob1;
    
    % a = linspace(-1.5,0.5,33);
    
    tols = [0 0.1 0.2 0.4 0.6 0.8 1 1.2];
    min_ratio(1:length(tols)) = 0;
    astar(1:length(tols)) = 0;
    kldiv(1:length(tols)) = 0;
    nkl = 10000000;
    umin = 0;
    umax = 7;
    
    n_MC = 100000;
    % c1 = c*c2
    c = 10;
    n0 = fix(n_MC*c/(c+1));
    n1 = n0;
    
    for t = 1:length(tols)
        fprintf('tol: %f\n', tols(t));
        prob2 = 1-normcdf(l1-tols(t));
        Q2 = @(x) l1-tols(t)-x;
        q2 = @(x) ((Q2(x)<0).*normpdf(x))/prob2;
        
        s_MC = rs(umin,umax,q2,n_MC);
        Q0s_MC = Q0(s_MC(:))<0;
        w_MC = normpdf(s_MC)./q2(s_MC);
        wQ0s_MC = w_MC.*Q0s_MC;
        
        Q0s = Q0s_MC(1:n0);
        w0 = normpdf(s_MC(1:n0))./q2(s_MC(1:n0));
        wQ0s = w0.*Q0s;
        Q1s = Q1(s_MC(1:n1))<0;
        w1 = normpdf(s_MC(1:n1))./q2(s_MC(1:n1));
        wQ1s = w1.*Q1s;

        v0_MC = var(wQ0s_MC)/n_MC;
        cov_Q0Q1 = cov(wQ0s,wQ1s);
        cov_MC_Q0Q1 = cov_Q0Q1/n0;
        v_MC_Q0 = cov_MC_Q0Q1(1,1);
        v_MC_Q1 = cov_MC_Q0Q1(2,2);
        % v_CV = v_MC_Q0 + a.^2*v_MC_Q1 + 2*a*cov_MC_Q0Q1(1,2);
        
        astar(t) = -cov_MC_Q0Q1(1,2)/v_MC_Q1;
        min_v_CV = v_MC_Q0 + astar(t)^2*v_MC_Q1 + 2*astar(t)*cov_MC_Q0Q1(1,2);
        min_ratio(t) = min_v_CV/v0_MC;
        
        kldiv(t) = kl(nkl,umin,umax,q1,q2);
    end
    
    set(0,'defaultLineLineWidth',0.7);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','on',...
        'PaperPositionMode','auto');
    hold on
    plot(l1-tols,min_ratio,'-s',l1-tols,kldiv,'-d')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times','Box','on')
    legend({'$\min(v_{CV}/v_0)$','KL Divergence'},...
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
    fn(1) = "tol_min_rs.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file
        print(figs(k), '-depsc2', sprintf('%s', append(path,fn(k))))
    end
    
end

function kldiv = kl(nkl,umin,umax,q1,q2)
    u = umin + (umax-umin)*rand(nkl,1);
    sample_value = q1(u);
    max_value = max(sample_value);
    accepted = rand(nkl,1)<(sample_value/max_value);
    s = u(accepted,:);
    if length(s) < 10^5
        error('kldiv: length(s) < 10^5');
    end
    kldiv = mean(log(q1(s(:)))-log(q2(s(:))));
end

function s = rs(umin,umax,q,n)
    rnsamples = n*50;
    u = umin + (umax-umin)*rand(rnsamples,1);
    sample_value = q(u);
    max_value = max(sample_value);
    accepted = rand(rnsamples,1) < (sample_value/max_value);
    s = u(accepted,:);
    if length(s) < n
        error('length(s) < n');
    else
        s = s(1:n);
    end
end