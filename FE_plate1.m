function wmax = FE_plate1(nelx,nely,x)
    E = 10000;
    poisson = 0.30;
    kapa = 5 / 6;
    L = 1;
    
    numberElements = nelx * nely;
    [nodeCoordinates, elementNodes] = ...
        rectangularMesh(L, L, nelx, nely);
    xx = nodeCoordinates(:, 1);
    yy = nodeCoordinates(:, 2);
    numberNodes = size(xx, 1);
    nn = flipud(reshape(1:numberElements,nely,nelx)');
    nnn = mat2cell(nn,[nely/2 nelx/2],[nely/2 nelx/2]);
    elemNum(1:4,(nelx/2)^2) = 0;
    elemNum(1,:) = sort(nnn{2,1}(:));
    elemNum(2,:) = sort(nnn{2,2}(:));
    elemNum(3,:) = sort(nnn{1,1}(:));
    elemNum(4,:) = sort(nnn{1,2}(:));
    % GDof: global number of degrees of freedom
    GDof = 3 * numberNodes;
    
    nsamples = size(x,1);
    h = mean(x(:,1:4)/100,2)*ones(1,4);
    P = mean(x(:,5:8),2)*ones(1,4);
    
    % boundary conditions
    [prescribedDof, ~] = ...
        EssentialBC('cccc', GDof, xx, yy, nodeCoordinates, numberNodes);
    
    wmax(1:nsamples,1) = 0;
    parfor i = 1:nsamples
        [stiffness] = ...
            formStiffnessMatrixMindlinQ4(GDof, numberElements, ...
            elementNodes, numberNodes, nodeCoordinates,E,poisson,kapa,h(i,:),elemNum);
        [force] = ...
            formForceVectorMindlinQ4(GDof, numberElements, ...
            elementNodes, nodeCoordinates,P(i,:),elemNum);
        % solution
        displacements = solution(GDof, prescribedDof, stiffness, force);
        wmax(i) = max(displacements(1:numberNodes));
    end
    
end

