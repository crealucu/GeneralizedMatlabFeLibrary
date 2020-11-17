classdef Mesh < handle
% Mesh - define finite element mesh
% 
%  Syntax: mesh = Mesh
% 
%   Properties:
%        elemType: type of element (Point, Line, Triangle, Quadrilateral, ...)
%       elemOrder: order of finite element
%          nodeno: number of nodes
%          elemno: number of elements
%            node: list of nodal coordinates
%            elem: list of elements
%             nsd: number of spatial dimensions
%             nne: number of nodes in an element
%         bndElem: boundary element (BoundaryElement)
% 
%  Methodes:
%     plot_mesh(b, e, disp) - plot mesh, only plot boundary elements
%                           - b: face color (0~1)
%                           - e: edge color, 'b', 'k', 'g', 'm', 'r', ...
%                           - disp: displacement in case to plot deformed shape
%     [new_node, new_elem] = increase_element_order
%                           - add more nodes and increases element order
%                           - new_node: added list of nodes
%                           - new_elem: added part of elements
%                             must shift node IDs before add this list to the element
%                             e.g.) mesh.elem = [mesh.elem, new_elem + mesh.nodeno];
%
%  Other m-files required: femlib, BoundaryElement
% 
%  See also: BoundaryElement, MeshGenerator
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 16-Oct-2020; Last revision:
%
  properties
    elemType; elemOrder; nodeno; elemno; node; elem; nsd, nne;
    bndElem;
  end
  methods
    function obj = Mesh
      obj.elemType  = [];
      obj.elemOrder = -1;
      obj.nodeno    = 0;
      obj.elemno    = 0;
      obj.node      = [];
      obj.elem      = [];
      obj.nsd       = -1;
      obj.bndElem   = BoundaryElement;
    end
    function h = plot_mesh(obj, b, e, disp)
      if(nargin<2)
        b = 0.9;
        e = 'k';
      end
      if(nargin<3)
        e = 'k';
      end
      if(nargin<4)
        disp = 0;
      end
      switch(obj.nsd)
        case 1
          h = plot_1D_mesh(obj, e, disp);
        case 2
          h = plot_2D_mesh(obj, b, e, disp);
        case 3
          h = plot_3D_mesh(obj, b, e, disp);
      end
    end
    function [new_node, new_elem] = increase_element_order(obj)
      eid = (1:obj.elemno)';
      edgeIDs = obj.bndElem.elemSideIDs.edgeNodeIDs(obj.elemType, obj.elemOrder);
      [edgeno, edge_nne] = size(edgeIDs);
      edges = zeros(obj.elemno*edgeno, edge_nne + 3);
      for ia = 0: edgeno-1
        edges(obj.elemno*ia + eid, :) = [obj.elem(:, edgeIDs(ia+1, :)), eid, ones(obj.elemno, 1) + ia, -ones(obj.elemno, 1)];
      end
      [oface, ~] = sort(edges(:, 1:end-3), 2);
      [~, ix, jx] = unique(oface, 'rows');
      
      new_node_no = numel(ix);
      new_node = zeros(new_node_no, obj.nsd);
      new_elem = zeros(obj.elemno, edgeno);
      for ia = 1: new_node_no
        fid = ix(ia);
        n = edges(fid, 1:end-3);
        x = obj.node(n, :);
        new_node(ia, :) = 0.5*sum(x);
      end
      
      for ia = 0: edgeno-1
        new_elem(:, ia+1) = jx(obj.elemno*ia + eid);
      end
    end
    function export_mesh_to_stl(obj, fn)
      [im, in] = size(obj.bndElem.elem);      
      bnd = obj.bndElem.elem(:);      
      
      e = unique(bnd);
      n = zeros(obj.nodeno, 1);
            
      n(e) = (1:numel(e))';
      bnd = reshape(n(bnd), im, in);
      
      TR = triangulation(bnd, obj.node(n>0, :));
      stlwrite(TR, fn, 'text');
    end
    function import_mesh_by_node_element(obj, N, E, kill_set_surface)
      if(nargin<4)
        kill_set_surface = false;
      end
      
      [obj.elemno, obj.nne] = size(E);
      [obj.nodeno, obj.nsd] = size(N);
      
      obj.elemOrder = 0;
      
      unknown = true;
      switch(obj.nsd)
        case 1
          if(obj.nne == 0)
            obj.elemType = EnumElementType.Point;
            unknown = false;
          else
            obj.elemType = EnumElementType.Line;
            obj.elemOrder = obj.nne - 2;
            unknown = false;
          end
        case 2
          switch(obj.nne)
            case 3
              obj.elemType = EnumElementType.Triangle;
              unknown = false;
            case 4
              obj.elemType = EnumElementType.Quadrilateral;
              unknown = false;
            case 6
              obj.elemType = EnumElementType.Triangle;
              obj.elemOrder = 1;
              unknown = false;
            case 8
              obj.elemType = EnumElementType.Quadrilateral;
              obj.elemOrder = 1;
              unknown = false;
          end
        case 3
          switch(obj.nne)
            case 4
              obj.elemType = EnumElementType.Tetrahedron;
              unknown = false;
            case 8
              obj.elemType = EnumElementType.Hexahedron;
              unknown = false;
            case 10
              obj.elemType = EnumElementType.Tetrahedron;
              obj.elemOrder = 1;
              unknown = false;
            case 20
              obj.elemType = EnumElementType.Hexahedron;
              obj.elemOrder = 1;
              unknown = false;
          end
      end
      
      if(unknown)
        fprintf('Unknown element type\n');
        obj.elemno = 0;
        obj.nne    = 0;
        obj.nodeno = 0;
        obj.nsd    = -1;
        obj.elemOrder = -1;
        return;
      end
      obj.elem = E;
      obj.node = N;
      if(kill_set_surface)
        return;
      end
      obj.bndElem.ConstructBoundaryElement(obj);
    end
    function import_mesh_by_node(obj, N, kill_set_surface)
      if(nargin<3); kill_set_surface = false; end
      E = delaunayn(N,{'Qt','Qbb','Qc','Qz'});
      obj.import_mesh_by_node_element(N, E, kill_set_surface);
    end
    function import_stl_mesh(obj, fn, eOrder)
      if(nargin<3)
        eOrder = 0;
      end
      
      if(eOrder == 0)
        o = 'linear';
      else
        o = 'quadratic';
      end
      model = createpde;
      importGeometry(model, fn);
      generateMesh(model, 'GeometricOrder', o);
      obj.import_mesh_by_node_element(model.Mesh.Nodes', model.Mesh.Elements');
    end
  end
end
function h = plot_1D_mesh(mesh, e, u)
  X = mesh.node + u;
  
  nid1 = mesh.elem(:, 1);
  nid2 = mesh.elem(:, 2);
  if(mesh.elemOrder == 1)
    nid1 = mesh.elem(:, [1,3])';
    nid2 = mesh.elem(:, [2,3])';
  end
  
  nid = mesh.elem(:, [1,2]);
  nid = unique(nid(:));
  h = figure; hold on;
  x = [X(nid1(:), 1), X(nid2(:), 1)];
  plot(x, x.*0, ['-', e]);
  plot(x, x.*0, '.k');  
  plot(X(nid, 1), nid.*0, 'b.');

  hold off
  axis equal;
  axis off;
end
function h = plot_2D_mesh(mesh, b, e, u)
  x = mesh.node + u;

  nid = 1:mesh.nne;
  if(mesh.elemOrder == 1)
    ids = 1:(mesh.nne/2);
    nid(1:2:end) = ids;
    nid(2:2:end) = (mesh.nne/2) + ids;
  end
  h = patch('Faces',mesh.elem(:, nid),'Vertices',x,'FaceColor','green');
  hold on;
  plot(mesh.node(:, 1), mesh.node(:, 2), 'b.');
  hold off;
  set(h,'facecolor', b*ones(1,3), 'edgecolor', e);
  axis equal
  axis off
end
function h = plot_3D_mesh(mesh, b, e, u)
  x = mesh.node + u;
  
  nid = 1:mesh.bndElem.nne;
  if(mesh.elemOrder == 1)
    ids = 1:(mesh.bndElem.nne/2);
    nid(1:2:end) = ids;
    nid(2:2:end) = (mesh.bndElem.nne/2) + ids;
  end

  h = patch('Faces',mesh.bndElem.elem(:, nid),'Vertices',x,'FaceColor','green');
  set(h,'facecolor', b*ones(1,3), 'edgecolor', e);
  hold on;
  n = unique(mesh.bndElem.elem(:));
  plot3(mesh.node(n, 1), mesh.node(n, 2), mesh.node(n, 3), 'b.');
  hold off;
  view(42,14);
  axis equal
  axis off
end