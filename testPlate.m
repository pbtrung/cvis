clear all;
format long;

% materials
E = 10000;
poisson = 0.30;
kapa = 5 / 6;

% matrix C
% bending part
% C_bending = I * E / (1 - poisson^2) * ...
%     [1, poisson, 0; poisson, 1, 0; 0, 0, (1 - poisson) / 2];
% shear part
% C_shear = kapa * thickness * E / 2 / (1 + poisson) * eye(2);

% Mesh generation
L = 1;
numberElementsX = 30;
numberElementsY = 30;
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

P = 2*ones(1,4);
h = 0.05*ones(1,4);
% boundary conditions
[prescribedDof, activeDof] = ...
        EssentialBC('cccc', GDof, xx, yy, nodeCoordinates, numberNodes);
[stiffness] = ...
        formStiffnessMatrixMindlinQ4(GDof, numberElements, ...
        elementNodes, numberNodes, nodeCoordinates,E,poisson,kapa,h,elemNum);
[force] = ...
    formForceVectorMindlinQ4(GDof, numberElements, ...
    elementNodes, nodeCoordinates,P,elemNum);
% solution
displacements = solution(GDof, prescribedDof, stiffness, force);
[M,I]=max(displacements(1:numberNodes))