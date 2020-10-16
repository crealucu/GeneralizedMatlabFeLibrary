classdef MeshGenerator < Mesh
% MeshGenerator - class to generate mesh
%
%  Syntax: mesh = MeshGenerator
% 
%  Properties:
%     all properties from FemLib
% 
%  Methodes:
%     generate_box(xmin, xmax, nx, dof, isTet, eOrder) - generate box mesh
%       - xmin  : minimum x, y, z coordinates
%       - xmax  : maximum x, y, z coordinates
%       - nx    : number of element in each direction(x, y, and z)
%       - dof   : degree of freedom
%       - isTet : if true, triangle or tetrahedron mesh otherwize line, quadrilateral, or hexahedron
%       - eOrder: 0: linear, 1: quadratic element, currently support upto quadratic elements
% 
%  Other m-files required: Mesh.m, BoundaryElement.m, femlib
% 
%  See also: Mesh, BoundaryElement
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 16-Oct-2020; Last revision:
%
  properties
  end
  methods
    function generate_box(obj, xmin, xmax, nx, dof, isTet, eOrder)

      if(nargin < 6); isTet = true; end
      if(nargin < 7); eOrder = 0; end
      
      obj.nsd = dof;

      % generate node list
      switch(dof)
        case 1
          obj.node = generate_node_1D(xmax, xmin, nx);
        case 2
          obj.node = generate_node_2D(xmax, xmin, nx);
        case 3
          obj.node = generate_node_3D(xmax, xmin, nx);          
      end
      obj.nodeno = size(obj.node, 1);
      
      % generate element list
      obj.elemOrder = 0;
      if(isTet && dof > 1)
        if(dof==3)
          tmp = delaunayTriangulation(obj.node(:, 1),obj.node(:, 2), obj.node(:, 3));
          obj.elemType = EnumElementType.Tetrahedron;
          obj.elemOrder = 0;
        else
          tmp = delaunayTriangulation(obj.node(:, 1),obj.node(:, 2));
          obj.elemType = EnumElementType.Triangle;
        end
        obj.elem = tmp.ConnectivityList;
      else
        obj.elem = hexa_element(nx, dof);
        switch(dof)
          case 1
            obj.elemType = EnumElementType.Line;
          case 2
            obj.elemType = EnumElementType.Quadrilateral;
          case 3
            obj.elemType = EnumElementType.Hexahedron;
        end
      end
      
      [obj.elemno, obj.nne] = size(obj.elem);
      if(eOrder == 1)
        [N, E] = obj.increase_element_order;
        obj.node = [obj.node; N];
        obj.elem = [obj.elem, E+obj.nodeno];
        obj.elemOrder = 1;
        obj.nodeno = size(obj.node);
        obj.nne    = size(obj.elem, 2);
      end
      obj.bndElem.ConstructBoundaryElement(obj);
    end
  end
end

function elem = hexa_element(nx, dof)

  elemno = 1;
  nne = 1;
  
  for ia = 1: dof
    elemno = elemno*nx(ia);
    nne = nne*2;
  end  
  
  elem = zeros(elemno, nne);
  
  NxNy = (nx(1)+1)*(nx(2)+1);
  Nx   = (nx(1)+1);
  
  switch(dof)
    case 1
      for ia=0: 1: nx(1)-1
        id = ia + 1;
        elem(id, 1) = ia+1;
        elem(id, 2) = ia+2;
      end
    case 2
      for ia=0: 1: nx(1)-1
        for ib=0: 1: nx(2)-1
          id = nx(1)*ib + ia + 1;
          elem(id, 1) =     ib*Nx + ia+1;
          elem(id, 2) =     ib*Nx + ia+2;
          elem(id, 3) = (ib+1)*Nx + ia+2;
          elem(id, 4) = (ib+1)*Nx + ia+1;
        end
      end
    case 3
      for ia = 0: nx(1)-1
        for ib = 0: nx(2)-1
          for ic = 0: nx(3)-1
            id = nx(1)*nx(2)*ic + nx(1)*ib + ia + 1;
            elem(id, 1) =     ic*NxNy +     ib*Nx + ia+1;
            elem(id, 2) =     ic*NxNy +     ib*Nx + ia+2;
            elem(id, 3) =     ic*NxNy + (ib+1)*Nx + ia+2;
            elem(id, 4) =     ic*NxNy + (ib+1)*Nx + ia+1;
            elem(id, 5) = (ic+1)*NxNy +     ib*Nx + ia+1;
            elem(id, 6) = (ic+1)*NxNy +     ib*Nx + ia+2;
            elem(id, 7) = (ic+1)*NxNy + (ib+1)*Nx + ia+2;
            elem(id, 8) = (ic+1)*NxNy + (ib+1)*Nx + ia+1;
          end
        end
      end
  end
end

function node = generate_node_1D(xmax, xmin, nx)
  node = zeros((nx(1)+1), 1);
  for ia = 0: 1: nx(1)
    node(ia+1) = ia*(xmax(1) - xmin(1))/nx(1) + xmin(1);
  end
end

function node = generate_node_2D(xmax, xmin, nx)
  node = zeros((nx(1)+1)*(nx(2)+1), 2);
  for ib = 0: 1: nx(2)
    y = ib*(xmax(2) - xmin(2))/nx(2) + xmin(2);
    for ia = 0: 1: nx(1)
      x = ia*(xmax(1) - xmin(1))/nx(1) + xmin(1);
      nid = ib*(nx(1)+1) + ia+1;
      node(nid, :) = [x, y];
    end
  end
end

function node = generate_node_3D(xmax, xmin, nx)
  node = zeros((nx(1)+1)*(nx(2)+1)*(nx(3)+1), 3);
  for ic = 0: nx(3)
    z = ic*(xmax(3) - xmin(3))/nx(3) + xmin(3);
    for ib = 0: nx(2)
      y = ib*(xmax(2) - xmin(2))/nx(2) + xmin(2);
      for ia = 0: nx(1)
        x = ia*(xmax(1) - xmin(1))/nx(1) + xmin(1);
        nid = ic*(nx(1)+1)*(nx(2)+1) + ib*(nx(1)+1) + ia+1;
        node(nid, :) = [x, y, z];
      end
    end
  end
end