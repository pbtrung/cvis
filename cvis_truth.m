function cvis_truth(model)

    format long;
    rng('default');

    % materials
    E = 10920; 
    poisson = 0.30; 
    kapa = 5/6; 
    thickness = 10^-4; 
    I = thickness^3/12;

    % constitutive matrix
    % bending part
    C_bending = I*E/(1-poisson^2)* ...
        [1 poisson 0;poisson 1 0;0 0 (1-poisson)/2];
    % shear part
    C_shear = kapa*thickness*E/2/(1+poisson)*eye(2);

    % load
    P = -1;

    % mesh generation
    L = 1;
    if model == 0
        numberElementsX = 30;
        numberElementsY = 30;
    elseif model == 1
        numberElementsX = 10;
        numberElementsY = 10;
    end
    if numberElementsX ~= numberElementsY || mod(numberElementsX,2) ~= 0
        error("Check numberElementsX and numberElementsY")
    end
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

    Pmin = 1;
    Pmax = 2;
    hmin = 1;
    hmax = 2;
    % boundary conditions
    [prescribedDof,~] = ...
        EssentialBC('cccc',GDof,xx,yy,nodeCoordinates,numberNodes);

    nsamples = 10^6;
    randP = rand(nsamples,4);
    randh = rand(nsamples,4);
    midNode = (numberNodes+1)/2;
    umax(1:nsamples) = 0;

    parfor i = 1:nsamples
        P = Pmin+(Pmax-Pmin)*randP(i,:);
        h = (hmin+(hmax-hmin)*randh(i,:))/20;

        % computation of the system stiffness matrix and force vector
        [stiffness] = ...
            formStiffnessMatrixMindlin(GDof,numberElements, ...
            elementNodes,numberNodes,nodeCoordinates,C_shear, ...
            C_bending,'Q4','complete','reduced');

        [force] = ...
            formForceVectorMindlin(GDof,numberElements, ...
            elementNodes,numberNodes,nodeCoordinates,P,'Q4','reduced');

        % solution
        displacements = solution(GDof,prescribedDof,stiffness,force);

        [~,idx] = max(displacements(1:numberNodes));
        umax(i) = displacements(midNode);
        if idx ~= midNode
            error('idx ~= midNode')
        end
    end
    
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    writematrix(umax,fullfile(repath,'umax0.txt'));

end