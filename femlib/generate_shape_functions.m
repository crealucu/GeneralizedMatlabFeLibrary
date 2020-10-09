function generate_shape_functions()
% generate_shape_functions - generate shape functions for femlib
%
%  Syntax:  generate_shape_functions
%  Outputs:
%     \param[out] ShapeFunctions.mat file will be created
%  Other m-files required: ShapeFunction.m
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 04-Oct-2020; Last revision:
  fn = 'ShapeFunctions.mat';
  elem = {%element type, element order
          EnumElementType.Point, [0]    
          EnumElementType.Line,  [0, 1, 2]
          EnumElementType.Triangle, [0, 1]
          EnumElementType.Quadrilateral, [0, 1]
          EnumElementType.Tetrahedron, [0, 1]
          EnumElementType.Hexahedron, [0, 1]
          };

  for ia = 1: size(elem, 1)
    fprintf('Generate shape function for %s element. Order:', elem{ia, 1});
    tmpN = cell(max(elem{ia, 2}) + 1, 1);
    tmpdN = cell(max(elem{ia, 2}) + 1, 1);
    for ib = 1: numel(elem{ia, 2})
      fprintf('%d.', elem{ia, 2}(ib));
      sf = ShapeFunction(elem{ia, 1}, elem{ia, 2}(ib));
      sf.compute_shape_function;
       tmpN{ elem{ia, 2}(ib)+1} = matlabFunction( sf.N, 'Vars', sf.X);
       tmpdN{elem{ia, 2}(ib)+1} = matlabFunction(sf.dN, 'Vars', sf.X);
    end
     N.(sprintf('%s', elem{ia, 1})) = tmpN;
    dN.(sprintf('%s', elem{ia, 1})) = tmpdN;
    fprintf('\n');
  end
  save(fn, 'N', 'dN');
end