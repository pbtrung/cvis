function plot_plate()

    format long;

    % materials
    E = 10000;
    poisson = 0.30;
    kapa = 5/6;

    % mesh generation
    L = 1;
    numberElementsX = 30;
    numberElementsY = 30;
    numberElements = numberElementsX*numberElementsY;
    %
    [nodeCoordinates, elementNodes] = ...
        rectangularMesh(L,L,numberElementsX,numberElementsY,'Q4');
    xx = nodeCoordinates(:,1);
    yy = nodeCoordinates(:,2);

    numberNodes = size(xx,1);
    nn = flipud(reshape(1:numberElements,numberElementsY,numberElementsX)');
    nnn = mat2cell(nn,[numberElementsY/2 numberElementsX/2],[numberElementsY/2 numberElementsX/2]);
    elemNum(1:4,(numberElementsX/2)^2) = 0;
    elemNum(1,:) = sort(nnn{2,1}(:));
    elemNum(2,:) = sort(nnn{2,2}(:));
    elemNum(3,:) = sort(nnn{1,1}(:));
    elemNum(4,:) = sort(nnn{1,2}(:));
    % GDof: global number of degrees of freedom
    GDof = 3*numberNodes;

    P = ones(1,4);
    h = 0.05*ones(1,4);
    % boundary conditions
    [prescribedDof,~] = ...
            EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);
    % computation of the system stiffness matrix and force vector
    stiffness = ...
        formStiffnessMatrixMindlin_R(GDof,...
        elementNodes,numberNodes,nodeCoordinates,...
        'Q4','complete','reduced',E,poisson,kapa,h,elemNum);
    force = ...
        formForceVectorMindlin_R(GDof,...
        elementNodes,nodeCoordinates,P,'Q4','reduced',elemNum);

    % solution
    displacements = solution(GDof,prescribedDof,stiffness,force);

    set(0,'defaultLineLineWidth',0.4);
    set(0,'defaultLineMarkerSize',2);
    fontsize = 8;
    width = 2.5;
    height = 2.5;

    fig = figure('Units','inches',...
        'Position',[0 0 width height],...
        'PaperPositionMode','auto');
    % plot3(xx,yy,displacements(1:numberNodes),'-.b',yy',xx',displacements(1:numberNodes),'-.b')
    xx = reshape(xx,numberElementsY+1,numberElementsX+1);
    yy = reshape(yy,numberElementsY+1,numberElementsX+1);
    d = reshape(displacements(1:numberNodes),numberElementsY+1,numberElementsX+1);
    surf(xx,yy,d,'FaceColor',[0.3010 0.7450 0.9330])
    set(gca,...
        'Units','normalized',...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',fontsize,...
        'FontName','Times',...
        'ZDir','reverse','ZLimSpec','Tight');
    print(fig, '-depsc2', "C:\Users\Nathan\Data\ss_high.eps")

    max(displacements(1:numberNodes))

end