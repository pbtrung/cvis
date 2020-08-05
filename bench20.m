function bench20()

    format long;
    rng('default');

    % materials
    E = 10000;
    poisson = 0.30;
    kapa = 5 / 6;

    % Mesh generation
    L = 1;
    numberElementsX = 20;
    numberElementsY = 20;
    if numberElementsX ~= numberElementsY || mod(numberElementsX,2) ~= 0
        error("Check numberElementsX and numberElementsY")
    end
    numberElements = numberElementsX * numberElementsY;

    [nodeCoordinates, elementNodes] = ...
        rectangularMesh(L, L, numberElementsX, numberElementsY);
    xx = nodeCoordinates(:, 1);
    yy = nodeCoordinates(:, 2);
    numberNodes = size(xx, 1);
    nn = flipud(reshape(1:numberElements,numberElementsY,numberElementsX)');
    nnn = mat2cell(nn,[numberElementsY/2 numberElementsX/2],[numberElementsY/2 numberElementsX/2]);
    elemNum(1:4,(numberElementsX/2)^2) = 0;
    elemNum(1,:) = sort(nnn{2,1}(:));
    elemNum(2,:) = sort(nnn{2,2}(:));
    elemNum(3,:) = sort(nnn{1,1}(:));
    elemNum(4,:) = sort(nnn{1,2}(:));
    % GDof: global number of degrees of freedom
    GDof = 3 * numberNodes;

    Pmin = 1;
    Pmax = 2;
    hmin = 1;
    hmax = 2;
    % % boundary conditions
    [prescribedDof, ~] = ...
            EssentialBC('cccc', GDof, xx, yy, nodeCoordinates, numberNodes);

    outer = 10;
    inner = 10;
    t(1:outer) = 0;

    for i = 1:outer
        randP = rand(inner,4);
        randh = rand(inner,4);
        
        tic
        for j = 1:inner
            P = Pmin+(Pmax-Pmin)*randP(j,:);
            h = (hmin+(hmax-hmin)*randh(j,:))/20;

            % computation of the system stiffness matrix and force vector
            [stiffness] = ...
                formStiffnessMatrixMindlinQ4(GDof, numberElements, ...
                elementNodes, numberNodes, nodeCoordinates,E,poisson,kapa,h,elemNum);

            [force] = ...
                formForceVectorMindlinQ4(GDof, numberElements, ...
                elementNodes, nodeCoordinates,P,elemNum);

            % solution
            displacements = solution(GDof, prescribedDof, stiffness, force);
        end
        t(i) = toc
    end
    
    mean(t)

end