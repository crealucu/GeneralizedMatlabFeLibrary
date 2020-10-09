classdef BndFemLib < FemLib
% BndFemLib - class for finite element boundary integration
%             inherits from FemLib
% 
%  Syntax:  fe = BndFemLib
% 
%  Properties:
%        all properties from FemLib + 
%             bnd: struct for the boundary integration
%                 bnd.etype: boundary element type
%                 bnd.eOder: boundary element order
%                 bnd.dN_of_kez: boundary element derivative of shape functions
%                 bnd.itgPts: list of boundary element integration points 
%                 bnd.sideId: boudnary element local side ID
%                 bnd.nodeIds: list of nodes compose the boundary element
%                 bnd.nne: number of nodes for the boundary element
%      elemSideId: list of boundary element node IDs for all element types 
%          normal: normal vector of the boundary
%  Methods:
%     set_an_element(coord, eid, eType, sideId, eOrder, iOrder) - set an finite element for the boundary integration
%                         - coord : list of nodal points
%                         - eid   : element ID
%                         - eType : element type (Line, Triangle, Quadrilateral, ...)
%                         - sideID: local side ID of an element
%                         - eOrder: element order (0: linear, 1: quadratic, 2: cubic, ...), mostly support upto quadratic
%                         - iOrder: integration order (qudrature rule)
%     ElemBasis(itgID)    - update quadrature rule at the given integration number (itgID) 
%
%  Example:
%     s = 0;
%     node = [0,0; 0.5,0; 0,0.5]; % triangle 
%     sideID = 2;
%     fe = BndFemLib;
%     fe.set_an_element(node, [], EnumElementType.Triangle, sideID);
%     for ip = 1: fe.quadRule.nint
%       fe.ElemBasis(ip);
%       s = s + fe.detJxW;
%     end
%     fprintf('line length the 2nd side of the triangle = %e\n', s);
% 
%  Other m-files required: FemLib.m ElementSide.m EnumElementType.m QuadratureRules.m
%  MAT-files required: ShapeFunctions.mat
% 
%  See also: BndFemLib
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 08-Oct-2020; Last revision:
  properties
    bnd;
    elemSideId;
    normal;
  end
  methods
    function obj = BndFemLib
      obj.funcs = load('ShapeFunctions.mat');
      eSide = ElementSide;
      obj.elemSideId = eSide.Side;
    end
    function set_an_element(obj, coord, eid, eType, sideId, eOrder, iOrder)
      if(nargin < 6); eOrder = 0;  end % will set default
      if(nargin < 7); iOrder = -1; end % will set default quadrature rule

      obj.node = coord;
      [obj.nne, obj.nsd] = size(obj.node);
      obj.curr_elmid = eid;  
      
      obj.sf = ShapeFunction(eType, eOrder);
      obj.N_of_kez = obj.funcs.N.(sprintf('%s', eType)){eOrder+1};
      obj.dN_of_kez = obj.funcs.dN.(sprintf('%s', eType)){eOrder+1};

      [obj.bnd.etype, obj.bnd.eOder] = VolToBnd(eType, eOrder);      
      obj.quadRule = QuadratureRules(obj.bnd.etype, iOrder);
      obj.bnd.dN_of_kez = obj.funcs.dN.(sprintf('%s', obj.bnd.etype)){obj.bnd.eOder+1};
      obj.bnd.itgPts = intPtsVolToBnd(eType, sideId, obj.quadRule.points);
      obj.bnd.sideId = sideId;
      [obj.bnd.nodeIds, obj.bnd.nne] = nodeIdVolToBnd(obj.elemSideId, eType, eOrder, sideId);
    end
    function ElemBasis(obj, itgID)
      %
      obj.curr_itg_id = itgID;
      if(obj.compute_N_dN(itgID) < 0)
        return;
      end

      if(obj.sf.nsd>1)
        iPt = obj.quadRule.points(itgID, :);
        dN_kez = obj.bnd.dN_of_kez(iPt(1), iPt(2), iPt(3));
        X = zeros(obj.bnd.nne, 3);
        X(:, 1:obj.nsd) = obj.node(obj.bnd.nodeIds, :);
        dX_kez = X'*dN_kez;
        tmp = ones(3, 2);
        tmp(1:obj.nsd, :) = dX_kez(1:obj.nsd, 1:2);
        n = cross(tmp(:, 1), tmp(:, 2));
        J = norm(n);
        obj.normal = n./J;
      else
        J = 1;
        n = [-1, 1];
        obj.normal = n(obj.bnd.sideId);
      end      
      obj.detJxW = J*obj.quadRule.weights(itgID);
    end
  end
end

function [bndEtype, bndEoder] = VolToBnd(eType, eOrder)
  bndEoder = eOrder;
  switch(eType)
    case EnumElementType.Line
      bndEtype = EnumElementType.Point;
      bndEoder = 0;
    case EnumElementType.Triangle
      bndEtype = EnumElementType.Line;
    case EnumElementType.Quadrilateral
      bndEtype = EnumElementType.Line;
    case EnumElementType.Tetrahedron
      bndEtype = EnumElementType.Triangle;
    case EnumElementType.Hexahedron
      bndEtype = EnumElementType.Quadrilateral;
  end  
end      
function itgPts = intPtsVolToBnd(eType, sideId, pts)
  ptsno = size(pts, 1);
  z = zeros(ptsno, 1);
  switch(eType)
    case EnumElementType.Line
      if(sideId == 1)
        itgPts = [-1, 0.0, 0.0];
      else
        itgPts = [1, 0.0, 0.0];
      end
    case EnumElementType.Triangle
      pts(:, 1) = 0.5*(pts(:, 1) + 1.0);
      switch(sideId)
        case 1
          itgPts = pts;
        case 2
          itgPts = [pts(:,1), 1-pts(:,1), z];
        case 3
          itgPts = [z, pts(:, 1), z];
      end
    case EnumElementType.Quadrilateral
      switch(sideId)
        case 1
          itgPts = pts;
        case 2
          itgPts = [z, pts(:, 1), z];
        case 3
          itgPts = pts;
        case 4
          itgPts = [z, pts(:, 1), z];
      end
    case EnumElementType.Tetrahedron
      switch(sideId)      
        case 1
          itgPts = pts;
        case 2
          itgPts = [pts(:, 1), z, pts(:, 2)];
        case 3
          itgPts = [pts(:, 1:2), 1-pts(:, 1)-pts(:, 2)];
        case 4
          itgPts = [z, pts(:, 1:2)];
      end
    case EnumElementType.Hexahedron
      switch(sideId)
        case 1
          itgPts = [pts(:, 1:2), z-1];
        case 2
          itgPts = [pts(:, 1:2), z+1];
        case 3
          itgPts = [pts(:, 1), z-1, pts(:, 2)];
        case 4
          itgPts = [z+1, pts(:, 1:2)];
        case 5
          itgPts = [pts(:, 1), z+1, pts(:, 2)];
        case 6
          itgPts = [z-1, pts(:, 1:2)];
      end
  end
end
function [nodeId, bnd_nne] = nodeIdVolToBnd(elmSideIds, eType, eOrder, sideId)
  faces = elmSideIds.(sprintf('%s', eType)){eOrder+1};
  nodeId = faces{sideId};
  bnd_nne = numel(nodeId);
end