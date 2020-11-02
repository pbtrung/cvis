function wmax = FE_plate(nelx,nely,E,poisson,kapa,L,x)
    
    % mesh generation
    numberElements = nelx*nely;
    %
    [nodeCoordinates, elementNodes] = ...
        rectangularMesh(L,L,nelx,nely,'Q4');
    xx = nodeCoordinates(:,1);
    yy = nodeCoordinates(:,2);

    numberNodes = size(xx,1);
    nn = flipud(reshape(1:numberElements,nely,nelx)');
    nnn = mat2cell(nn,[nely/2 nelx/2],[nely/2 nelx/2]);
    elemNum(1:4,(nelx/2)^2) = 0;
    elemNum(1,:) = sort(nnn{2,1}(:));
    elemNum(2,:) = sort(nnn{2,2}(:));
    elemNum(3,:) = sort(nnn{1,1}(:));
    elemNum(4,:) = sort(nnn{1,2}(:));
    % GDof: global number of degrees of freedom
    GDof = 3*numberNodes;
    
    nsamples = size(x,1);
    midNode = (numberNodes+1)/2;
    P = x(:,1:4);
    h = x(:,5:8)/20;
    
    % boundary conditions
    [prescribedDof,~] = ...
        EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);
    
    wmax(1:nsamples,1) = 0;
    parfor i = 1:nsamples
        tic
        % computation of the system stiffness matrix and force vector
        stiffness = ...
            formStiffnessMatrixMindlin_R(GDof,...
            elementNodes,numberNodes,nodeCoordinates,...
            'Q4','complete','reduced',E,poisson,kapa,h(i,:),elemNum);

        force = ...
            formForceVectorMindlin_R(GDof,...
            elementNodes,nodeCoordinates,P(i,:),'Q4','reduced',elemNum);

        % solution
        displacements = solution(GDof,prescribedDof,stiffness,force);
        wmax(i) = displacements(midNode);
        t = toc;
        fprintf('iter: %d, time: %f\n',i,t);
    end
    
end
