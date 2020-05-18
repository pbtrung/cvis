clear all;
format long;
rng('default');

% materials
E = 10000;
poisson = 0.30;
kapa = 5 / 6;

% Mesh generation
L = 1;
numberElementsX = 10;
numberElementsY = 10;
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
Pmax = 1.1;
hmin = 1;
hmax = 1.1;
% % boundary conditions
[prescribedDof, activeDof] = ...
        EssentialBC('cccc', GDof, xx, yy, nodeCoordinates, numberNodes);

outer = 1000;
inner = 1000;
maxNode = (numberNodes+1)/2;
umax(1:outer,1:inner) = 0;

for i = 1:outer
    randP = rand(inner,4);
    randh = rand(inner,4);
    tic
    parfor j = 1:inner
        fprintf('outer: %d, inner: %d\n',i,j);
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

        % deformed shape
        % figure
        % plot3(xx,yy,displacements(1:numberNodes),'.')
        umax(i,j) = max(displacements(1:numberNodes));
        if umax(i,j) ~= displacements(maxNode)
            error('Different node')
        end
    end
    toc
end
writematrix(umax,'Ex3_umax_cccc_10x10_10_11.txt');

m = max(umax(:));
fprintf('m: %f\n',m);
l = 0.05:0.05:1;
lm(1:length(l)) = 0;
Umean(1:length(l)) = 0;
Uvar(1:length(l)) = 0;

for j = 1:length(l)
    lm(j) = l(j)*m;
    Umean(j) = mean(mean(lm(j)-umax<0,2));
    Uvar(j) = var(mean(lm(j)-umax<0,2));
end

writematrix(lm','plate_lm.txt');
writematrix(Umean','plate_Umean.txt');
writematrix(Uvar','plate_Uvar.txt');