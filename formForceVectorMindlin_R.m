function force = ...
    formForceVectorMindlin_R(GDof,...
    elementNodes,nodeCoordinates,P,elemType,quadType,elemNum)

    % computation of force vector for Mindlin plate element

    % force: force vector
    force = sparse(GDof,1);

    % Gauss quadrature for bending part
    [gaussWeights,gaussLocations] = gaussQuadrature(quadType);

    % cycle for element
    for k = 1:length(P)
        elem = elemNum(k,:);
        for e = 1:length(elem)
            % indice : nodal connectivities for each element
            indice = elementNodes(elem(e),:);

            % cycle for Gauss point
            for q = 1:size(gaussWeights,1)
                GaussPoint = gaussLocations(q,:);
                GaussWeight = gaussWeights(q);
                xi = GaussPoint(1);
                eta = GaussPoint(2);

                % shape functions and derivatives
                [shapeFunction,naturalDerivatives] = ...
                    shapeFunctionsQ(xi,eta,elemType);

                % Jacobian matrix, inverse of Jacobian,
                % derivatives w.r.t. x,y
                [Jacob,~,~] = ...
                    Jacobian(nodeCoordinates(indice,:),naturalDerivatives);

                % force vector
                force(indice) = force(indice) + ...
                    shapeFunction*P(k)*det(Jacob)*GaussWeight;
            end  % end Gauss point loop

        end % end element loop
    end

end
