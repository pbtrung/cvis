function cvis_rs_2()
    
    format long;
    % rng('default');
    
    mux = 0;
    muy = 0;
    stdx = 1;
    stdy = 1;
    lx = mux+3*stdx;
    ly = muy+2*stdy;
    
    nx = @(x) normpdf(x,mux,stdx);
    ny = @(y) normpdf(y,muy,stdy);
    probx = integral(nx,lx,Inf)
    proby = integral(ny,ly,Inf)
    
    a = linspace(-4,4,10);
    ansamples = 1000;
    wQx(1:length(a),1:ansamples) = 0;
    wQy(1:length(a),1:ansamples) = 0;
    Q(1:length(a),1:ansamples) = 0;
    nsamples = 100;
    ss(1:length(a),1:ansamples,1:nsamples) = 0;
    umin = 0;
    umax = 7;
    m = 10;
    
    tol = -0.2;
    proby2 = integral(ny,ly-tol,Inf)
    q = @(y) ((ly-tol-y<0).*ny(y))./proby2;
    
    for i = 1:length(a)
        parfor j = 1:ansamples
            % rejection sampling
            samples = rs(nsamples,umin,umax,m,q);
            ss(i,j,:) = samples;
            
            Qx = (lx-samples)<0;
            w = nx(samples)./q(samples);
            wQx(i,j) = mean(w.*Qx);
            Qy = (ly-samples)<0;
            wQy(i,j) = mean(w.*Qy);
            Q(i,j) = wQx(i,j)+a(i)*(wQy(i,j)-proby);
        end
    end
    
    probx./mean(Q,2)
    t1 = var(wQx,0,2)
    t2 = var(wQy,0,2)
    t3 = cov(wQx(2,:),wQy(2,:))
        
    figure(1)
    hold on
    plot(a,var(Q,0,2),'-o',a,t1(2)+a.^2*t2(2)+2*a*t3(1,2),'--*')
    xlabel('alpha')
    hold off

%     qx = @(x) ((lx-x<0).*nx(x))./probx;
%     qy = @(y) ((ly-y<0).*ny(y))./proby;
%     figure(1)
%     hold on
%     X = linspace(0,5,1000);
%     plot(X,qx(X),'-',X,qy(X),'--',X,q(X))
%     histogram(ss(1,1,:),100,'Normalization','pdf');
%     legend('qx','qy','q')
%     hold off
    
end

function s = rs(nsamples_in,a,b,m,q)
    nsamples = 500*nsamples_in;
    y = a+(b-a)*rand(nsamples,1);
    fy = q(y);
    u = m*rand(nsamples,1);
    accepted = u <= fy;
    ss = y(accepted);
    s = ss(1:nsamples_in);
end