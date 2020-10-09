function verification_tests_volume_integration()
% verification_tests - run verification tests for each element 
% Each element will integrated to compute its volume. 
% All quadrature rules and element types will be tested.
% 
%  Syntax:  verification_tests_volume_integration
% 
%  Outputs:
%     print tests results
%
%  Other m-files required: femlib.m, QuadratureRules.m ShapeFunction.m
%                          cprintf.m
%  MAT-files required: ShapeFunctions.mat
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 05-Oct-2020; Last revision:

  tol = 1.0e-12;
  Line.e = { [1; 4]
    [1; 4; 2.5]
    [1; 4; 2; 3]};
  Line.V = 3;
  Line.iOrder = [0, 1, 2, 3];
  
  Triangle.e = {[1, 1; 3, 1; 1, 3]
    [1, 1; 3, 1; 1, 3; 2, 1; 2, 2; 1, 2]};
  Triangle.iOrder = [0, 1, 2, 3];
  Triangle.V = 2;
  
  Quadrilateral.e = {[0,0; 1, 0; 2, 2; 1, 2]
                     [0,0; 1, 0; 2, 2; 1, 2; 0.5, 0; 1.5, 1; 1.5, 2; 0.5, 1]};
  Quadrilateral.iOrder = [0, 1, 2];
  Quadrilateral.V = 2;

  Tetrahedron.e = {[0,0,0; 2,0,0; 0,2,0; 0,0,2]
                   [0,0,0; 2,0,0; 0,2,0; 0,0,2;1,0,0; 1,1,0; 0,1,0; 0, 0, 1; 1,0,1; 0, 1, 1]};
  Tetrahedron.iOrder = [0, 1, 2];
  Tetrahedron.V = 8/6;
  
  Hexahedron.e = {[0,0,0; 2,0,0; 2,2,0; 0,2,0; 0,0,2; 2,0,2; 2,2,2; 0,2,2]
                  [0,0,0; 2,0,0; 2,2,0; 0,2,0; 0,0,2; 2,0,2; 2,2,2; 0,2,2
                   1,0,0; 2,1,0; 1,2,0; 0,1,0; 1,0,2; 2,1,2; 1,2,2; 0,1,2
                   0,0,1; 2,0,1; 2,2,1; 0,2,1]};
  Hexahedron.iOrder = [0, 1, 2];
  Hexahedron.V = 8;
  eType = {'Line', Line
           'Triangle', Triangle
           'Quadrilateral', Quadrilateral
           'Tetrahedron', Tetrahedron
           'Hexahedron', Hexahedron};
  fe = FemLib;
  eid = 1;
  for iA = 1: size(eType, 1)
    for e = 1: numel(eType{iA, 2}.e)
      for itg = 1: numel(eType{iA, 2}.iOrder)
        eOrder = e - 1;
        fe.set_an_element(eType{iA, 2}.e{e}, eid, EnumElementType.(eType{iA, 1}), eOrder, eType{iA, 2}.iOrder(itg));
        V = 0;
        for ip = 1: fe.quadRule.nint
          fe.ElemBasis(ip);
          V = V + fe.detJxW;
        end
        err = abs(V - eType{iA, 2}.V)/eType{iA, 2}.V;
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