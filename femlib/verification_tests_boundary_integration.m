function verification_tests_boundary_integration()
% verification_tests - run verification tests for each element 
% Each element will integrated to compute its volume. 
% All quadrature rules and element types will be tested.
% 
%  Syntax:  verification_tests_boundary_integration
% 
%  Outputs:
%     print tests results
%
%  Other m-files required: FemLib.m, BndFemLib.m  
%                          QuadratureRules.m ShapeFunction.m
%                          cprintf.m
%  MAT-files required: ShapeFunctions.mat
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 05-Oct-2020; Last revision:

  tol = 1.0e-12;
  Line.e = { [1; 4]
    [1; 2.5; 4]
    [1; 2; 3; 4]};
  Line.S = 1;
  Line.iOrder = [0, 1, 2, 3];
  Line.sideId = 1;
  
  Triangle.e = {[1, 1; 3, 1; 1, 3]
                [1, 1; 3, 1; 1, 3; 2, 1; 2, 2; 1, 2]};
  Triangle.iOrder = [0, 1, 2, 3];
  Triangle.S = sqrt(8);
  Triangle.sideId = 2;
  
  Quadrilateral.e = {[0,0; 1, 0; 2, 2; 1, 2]
                     [0,0; 1, 0; 2, 2; 1, 2; 0.5, 0; 1.5, 1; 1.5, 2; 0.5, 1]};
  Quadrilateral.iOrder = [0, 1, 2];
  Quadrilateral.S = sqrt(5);
  Quadrilateral.sideId = 2;
  Tetrahedron.e = {[0,0,0; 2,0,0; 0,2,0; 0,0,2]
                   [0,0,0; 2,0,0; 0,2,0; 0,0,2;1,0,0; 1,1,0; 0,1,0; 0, 0, 1; 1,0,1; 0, 1, 1]};
  Tetrahedron.iOrder = [0, 1, 2];
  Tetrahedron.S = 4*sin(pi*2/3);
  Tetrahedron.sideId = 3;
  
  Hexahedron.e = {[0,0,0; 2,0,0; 2,2,0; 0,2,0; 0,0,2; 2,0,2; 2,2,2; 0,2,2]
                  [0,0,0; 2,0,0; 2,2,0; 0,2,0; 0,0,2; 2,0,2; 2,2,2; 0,2,2
                   1,0,0; 2,1,0; 1,2,0; 0,1,0; 1,0,2; 2,1,2; 1,2,2; 0,1,2
                   0,0,1; 2,0,1; 2,2,1; 0,2,1]};
  Hexahedron.iOrder = [0, 1, 2];
  Hexahedron.S = 4;
  Hexahedron.sideId = 1;
  eType = {'Line', Line
           'Triangle', Triangle
           'Quadrilateral', Quadrilateral
           'Tetrahedron', Tetrahedron
           'Hexahedron', Hexahedron
           };
  fe = BndFemLib;
  for iA = 1: size(eType, 1)
    for e = 1: numel(eType{iA, 2}.e)
      for itg = 1: numel(eType{iA, 2}.iOrder)
        eOrder = e - 1;
        fe.set_an_element(eType{iA, 2}.e{e}, [], EnumElementType.(eType{iA, 1}), eType{iA, 2}.sideId, eOrder, eType{iA, 2}.iOrder(itg));
        S = 0;
        for ip = 1: fe.quadRule.nint
          fe.ElemBasis(ip);
          S = S + fe.detJxW;
        end
        err = abs(S - eType{iA, 2}.S)/eType{iA, 2}.S;
        fprintf('Test for %s element with element order = %d and interation order = %d, err = %e', ...
          eType{iA, 1}, eOrder, eType{iA, 2}.iOrder(itg), err);
        if(err < tol)
          cprintf('blue', ': Pass\n');
        else
          cprintf('red', ': Failed\n');
        end
      end
    end
  end 
end