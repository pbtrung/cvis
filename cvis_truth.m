function cvis_truth(model)

    format long;
    rng('default');

    % materials
    E = 10000; 
    poisson = 0.30;
    kapa = 5/6;
    
    % mesh generation
    L = 1;
    if model == 0
        numberElementsX = 30;
        numberElementsY = 30;
    elseif model == 1
        numberElementsX = 10;
        numberElementsY = 10;
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
    t(1:nsamples) = 0;
    midNode = (numberNodes+1)/2;
    umax(1:nsamples) = 0;
    [filepath,~,~] = fileparts(matlab.desktop.editor.getActiveFilename);
    repath = fullfile(filepath,'results');
    
    parfor i = 1:nsamples
        tic
        P = Pmin+(Pmax-Pmin)*randP(i,:);
        h = (hmin+(hmax-hmin)*randh(i,:))/20;
        
        tic
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
        t(i) = toc;
        
        [~,idx] = max(displacements(1:numberNodes));
        umax(i) = displacements(midNode);
        if idx ~= midNode
            fprintf('iter: %d, time: %f, idx ~= midNode\n',i,t(i));
        else
            fprintf('iter: %d, time: %f, idx = midNode\n',i,t(i));
        end
    end
    
    if model == 0
        writematrix(umax,fullfile(repath,'umax_30x30_midNode.txt'));
        writematrix(t,fullfile(repath,'t_30x30_midNode.txt'));
    elseif model == 1
        writematrix(umax,fullfile(repath,'umax_10x10_midNode.txt'));
        writematrix(t,fullfile(repath,'t_10x10_midNode.txt'));
    end
    
    m = max(umax(:));
    l = 0.05:0.05:1;
    lm(1:length(l)) = 0;
    Umean(1:length(l)) = 0;

    for j = 1:length(l)
        lm(j) = l(j)*m;
        Umean(j) = mean(lm(j)-umax<0);
    end

    if model == 0
        writematrix(lm',fullfile(repath,'plate_lm_30x30_midNode.txt'));
        writematrix(Umean',fullfile(repath,'plate_Umean_30x30_midNode.txt'));
    elseif model == 1
        writematrix(lm',fullfile(repath,'plate_lm_10x10_midNode.txt'));
        writematrix(Umean',fullfile(repath,'plate_Umean_10x10_midNode.txt'));
    end

end