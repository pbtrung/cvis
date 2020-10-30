function [weights,locations] = gaussQuadrature(option)
% Gauss quadrature for 2D elements
% option 'third' (3x3)
% option 'complete' (2x2)
% option 'reduced'  (1x1)
% locations: Gauss point locations
% weights: Gauss point weights

switch option
    case 'third'
        locations = [-0.774596669241483 -0.774596669241483;
                      0.                -0.774596669241483;
                      0.774596669241483 -0.774596669241483;
                     -0.774596669241483  0.;
                      0.                 0.;
                      0.774596669241483  0.;
                     -0.774596669241483  0.774596669241483;
                      0.                 0.774596669241483;
                      0.774596669241483  0.774596669241483];
        weights = [0.555555555555556*0.555555555555556;
                   0.555555555555556*0.888888888888889;
                   0.555555555555556*0.555555555555556;
                   0.888888888888889*0.555555555555556; 
                   0.888888888888889*0.888888888888889;
                   0.555555555555556*0.888888888888889;
                   0.555555555555556*0.555555555555556;
                   0.555555555555556*0.888888888888889;
                   0.555555555555556*0.555555555555556];
        
    case 'complete'
        locations = [ -0.577350269189626 -0.577350269189626;
                       0.577350269189626 -0.577350269189626;
                       0.577350269189626  0.577350269189626;
                      -0.577350269189626  0.577350269189626];
        weights = [1;1;1;1];
        
    case 'reduced'
        locations = [0 0];
        weights = 4;
end

end  % end function gaussQuadrature