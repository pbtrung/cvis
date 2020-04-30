function [nodeCoordinates, elementNodes] ...
    = rectangularMesh(Lx, Ly, numberElementsX, numberElementsY)
% Generate Simple Rectangular Mesh
% Input Data: dimensions & No. of elements in x & y direction
%lx = 20;%input('Dimension in x-direction: ' );
%ly = 5;%input('Dimension in y-direction: ' );
%nelex = 10;%input('Enter No. of Elements in x-direction: ' );
%neley = 5;%input('Enter No. of Elements in y-direction: ' );
%Calculations
%No. of elements and increment in x & y direction
nnodex = numberElementsX + 1;
nnodey = numberElementsY + 1;
dx = Lx / numberElementsX;
dy = Ly / numberElementsY;
% Generate Coordinates
nelexy = numberElementsX * numberElementsY; % total no. of elements
nnodexy = nnodex * nnodey; % total no. of nodes
% zero matrices and vectors
xc = zeros(nnodexy, 1);
yc = zeros(nnodexy, 1);
nodedata = zeros(nnodexy, 3);
for row = 1:nnodey
    for col = 1:nnodex
        nt = (row - 1) * nnodex + col;
        xc(nt) = (col - 1) * dx;
        yc(nt) = (row - 1) * dy;
    end
end
% combining node data into a matrix nodedata=[no., xc, yc]
% and elementdata=[no., 4 nodes]
for i = 1:nnodexy
    nodedata(i, 1) = i;
    nodedata(i, 2) = xc(i);
    nodedata(i, 3) = yc(i);
end
nodeCoordinates(:, 1) = nodedata(:, 2);
nodeCoordinates(:, 2) = nodedata(:, 3);
elementNodes(1:nelexy, 1:4) = 0;
for elerow = 1:numberElementsY
    for elecol = 1:numberElementsX
        eleno = (elerow - 1) * (nnodex - 1) + elecol;
        node1 = eleno + (elerow - 1);
        node2 = eleno + (elerow);
        node3 = node2 + nnodex;
        node4 = node1 + nnodex;
        elementNodes(eleno, :) = [node1 node2 node3 node4];
    end
end

% Plot of the structure
% xc1=xc(1:1:11);yc1=yc(1:1:11);
% xc2=xc(12:1:22);yc2=yc(12:1:22);
% xc3=xc(23:1:33);yc3=yc(23:1:33);
% xc4=xc(34:1:44);yc4=yc(34:1:44);
% xc5=xc(45:1:55);yc5=yc(45:1:55);
% xc6=xc(56:1:66);yc6=yc(56:1:66);
%
% xc7=xc(1:11:56);yc7=yc(1:11:56);
% xc8=xc(2:11:57);yc8=yc(2:11:57);
% xc9=xc(3:11:58);yc9=yc(3:11:58);
% xc10=xc(4:11:59);yc10=yc(4:11:59);
% xc11=xc(5:11:60);yc11=yc(5:11:60);
% xc12=xc(6:11:61);yc12=yc(6:11:61);
% xc13=xc(7:11:62);yc13=yc(7:11:62);
% xc14=xc(8:11:63);yc14=yc(8:11:63);
% xc15=xc(9:11:64);yc15=yc(9:11:64);
% xc16=xc(10:11:65);yc16=yc(10:11:65);
% xc17=xc(11:11:66);yc17=yc(11:11:66);

% figure (1)
% plot (xc,yc,'-o','markersize', 2, 'markerfacecolor', 'b')
% plot (xc1,yc1,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc2,yc2,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc3,yc3,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc4,yc4,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc5,yc5,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc6,yc6,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc7,yc7,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc8,yc8,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc9,yc9,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc10,yc10,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc11,yc11,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc12,yc12,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc13,yc13,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc14,yc14,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc15,yc15,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc16,yc16,'-o','markersize', 2, 'markerfacecolor', 'b')
% hold on
% plot (xc17,yc17,'-o','markersize', 2, 'markerfacecolor', 'b')
% grid on
% axis equal
% elementdata
% nodedata
% NoOfElements = nelexy
% NoOfNodes= nnodexy

end