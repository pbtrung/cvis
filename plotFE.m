function plotFE(lx,ly,nelx,nely,nelx1,nely1,U,U1)
    
    lelx = lx/nelx;
    lely = ly/nely;
    x = 0:lelx:lx;
    y = 0:lely:ly;
    [X,Y] = meshgrid(x,y);
    
    figure
    plot(X,Y,'*')
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
    
end