function cvis_ce_rf_sample()
    
    format long;
    rng('default');
    
    lx = 0.6;
    ly = 0.2;
    nelx = 60;
    nely = 20;
    nelx1 = 30;
    nely1 = 10;
    Lx = 60;
    Ly = 20;
    a = 1;
    b = 2;
    
    [Psi,lambda,PsiE1] = KL(lx,ly,Lx,Ly,nelx,nely,2);
    
    force = 1;
    dof = 2*(nely+1)*nelx+2;
    
    neig = length(lambda);
    L = diag(sqrt(lambda(1:neig)));
    P = reshape(permute(Psi,[2 1 3]),nelx*nely,neig);
        
    Z = ones(neig,1)/2;
    LZ = L*Z;
    PZ = dot(P',repmat(LZ,1,nelx*nely));
    SS = reshape(PZ,nelx,nely)';
    E = a+(b-a)*normcdf(SS);
    tic
    U = FErf(lx,ly,nelx,nely,dof,force,E);
    t = toc;
    
    dof = 2*(nely1+1)*nelx1+2;
    SS1 = PsiE1*LZ;
    E1 = a+(b-a)*normcdf(SS1);
    tic
    U1 = FErf1(lx,ly,nelx1,nely1,dof,force,E1);
    t1 = toc;
    
    t
    t1
    t/t1
     
    lelx = lx/nelx;
    lely = ly/nely;
    x = 0:lelx:lx;
    y = 0:lely:ly;
    [X,Y] = meshgrid(x,y);
    Ux = reshape(U(1:2:end),nely+1,[]);
    Uy = reshape(U(2:2:end),nely+1,[]);
%     max(Ux,[],'all')
    max(Uy,[],'all')
    
    lelx1 = lx/nelx1;
    lely1 = ly/nely1;
    x1 = 0:lelx1:lx;
    y1 = 0:lely1:ly;
    [X1,Y1] = meshgrid(x1,y1);
    Ux1 = reshape(U1(1:2:end),nely1+1,[]);
    Uy1 = reshape(U1(2:2:end),nely1+1,[]);
%     max(Ux1,[],'all')
    max(Uy1,[],'all')
    
    set(0,'defaultLineLineWidth',0.4);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 3;
    height = 1;
    
    figs(1) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(X,Y+Uy/1000,'-b',X',(Y+Uy/1000)','-b')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'XAxisLocation','top','YAxisLocation','left',...
        'ydir','reverse','YLimSpec','Tight');
    fn(1) = "ss_high.eps";
    hold off
    
    figs(2) = figure('Units','inches',...
        'Position',[0 0 width height],...
        'visible','off',...
        'PaperPositionMode','auto');
    hold on
    plot(X1,Y1+Uy1/1000,'-b',X1',(Y1+Uy1/1000)','-b')
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'XAxisLocation','top','YAxisLocation','left',...
        'ydir','reverse','YLimSpec','Tight');
    fn(2) = "ss_low.eps";
    hold off
    
    for k = 1:length(figs)
        % print each figure in figs to a separate .eps file  
        print(figs(k), '-depsc2', sprintf('%s', fn(k)))
    end
    
end