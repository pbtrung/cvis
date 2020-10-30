function [K] = ...
    formStiffnessMatrixMindlin_R(GDof,...
    elementNodes,numberNodes,nodeCoordinates,...
    elemType,quadTypeB,quadTypeS,E,poisson,kapa,h,elemNum)
    % elemType: type of element Q4, Q8, Q9
    % quadTypeB: type of quadrature for bending
    % quadTypeS: type of quadrature for shear

    % computation of stiffness matrix for Mindlin plate element

    % K : stiffness matrix
    K = zeros(GDof);

    % Gauss quadrature for bending part
    [gaussWeights,gaussLocations] = gaussQuadrature(quadTypeB);

    % cycle for element
    for k = 1:length(h)
        I = h(k)^3 / 12;
        C_bending = I * E / (1 - poisson^2) * ...
        [1, poisson, 0; poisson, 1, 0; 0, 0, (1 - poisson) / 2];
        elem = elemNum(k,:);
        for e = 1:length(elem)
            % indice : nodal connectivities for each element
            % elementDof: element degrees of freedom
            indice = elementNodes(e,:);
            elementDof = [indice indice+numberNodes indice+2*numberNodes];
            ndof = length(indice);

            % cycle for Gauss point
            for q = 1:size(gaussWeights,1)
                GaussPoint = gaussLocations(q,:);
                xi = GaussPoint(1);
                eta = GaussPoint(2);

                % shape functions and derivatives
                [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);

                % Jacobian matrix, inverse of Jacobian,
                % derivatives w.r.t. x,y
                [Jacob,~,XYderivatives] = ...
                    Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

                % [B] matrix bending
                B_b = zeros(3,3*ndof);
                B_b(1,ndof+1:2*ndof)   = XYderivatives(:,1)';
                B_b(2,2*ndof+1:3*ndof) = XYderivatives(:,2)';
                B_b(3,ndof+1:2*ndof)   = XYderivatives(:,2)';
                B_b(3,2*ndof+1:3*ndof) = XYderivatives(:,1)';

                % stiffness matrix bending
                K(elementDof,elementDof) = K(elementDof,elementDof) + ...
                    B_b'*C_bending*B_b*gaussWeights(q)*det(Jacob);

            end  % Gauss point
        end    % element
    end

    % shear stiffness matrix

    % Gauss quadrature for shear part
    [gaussWeights,gaussLocations] = gaussQuadrature(quadTypeS);

    % cycle for element
    for k = 1:length(h)
        C_shear = kapa * h(k) * E / 2 / (1 + poisson) * eye(2);
        elem = elemNum(k,:);
        for e = 1:length(elem)
            % indice : nodal connectivities for each element
            % elementDof: element degrees of freedom
            indice = elementNodes(e,:);
            elementDof = [indice indice+numberNodes indice+2*numberNodes];
            ndof = length(indice);

            % cycle for Gauss point
            for q = 1:size(gaussWeights,1)
                GaussPoint = gaussLocations(q,:);
                xi = GaussPoint(1);
                eta = GaussPoint(2);

                % shape functions and derivatives
                [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);

                % Jacobian matrix, inverse of Jacobian,
                % derivatives w.r.t. x,y
                [Jacob,~,XYderivatives] = ...
                    Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

                % [B] matrix shear
                B_s = zeros(2,3*ndof);
                B_s(1,1:ndof)         = XYderivatives(:,1)';
                B_s(2,1:ndof)         = XYderivatives(:,2)';
                B_s(1,ndof+1:2*ndof)  = shapeFunction;
                B_s(2,2*ndof+1:3*ndof)= shapeFunction;

                % stiffness matrix shear
                K(elementDof,elementDof) = K(elementDof,elementDof) + ...
                    B_s'*C_shear*B_s*gaussWeights(q)*det(Jacob);
            end  % gauss point
        end    % element
    end

end