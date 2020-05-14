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
Pmax = 100;
hmin = 1;
hmax = 100;
% % boundary conditions
[prescribedDof, activeDof] = ...
        EssentialBC('cccc', GDof, xx, yy, nodeCoordinates, numberNodes);

outer = 1000;
inner = 1000;
umax(1:outer,1:inner) = 0;

for i = 1:outer
    randP = rand(inner,4);
    randh = rand(inner,4);
    tic
    parfor j = 1:inner
        fprintf('outer: %d, inner: %d\n',i,j);
        P = mean(Pmin+(Pmax-Pmin)*randP(j,:))*ones(1,4);
        h = mean((hmin+(hmax-hmin)*randh(j,:))/100)*ones(1,4);

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
    end
    toc
end
writematrix(umax,'Ex3_umax_cccc_10x10.txt');

m = max(umax(:));
fprintf('m: %f\n',m);
l = 0.5:0.05:1;
for j = 1:length(l)
    fprintf('iter: %d, l: %f, l*m: %f\n',j,l(j),l(j)*m);
    mean(mean(l(j)*m-umax<0,2))
    var(mean(l(j)*m-umax<0,2))
end