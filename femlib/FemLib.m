classdef FemLib < handle
% FemLib - class for finite element volume integration
% 
%  Syntax:  fe = FemLib
% 
%  Properties:
%           funcs: all shape functions and their derivative as a function of
%                  paramatric coordinate (ksi, eta, zeta)
%        quadRule: quadrature rule for given element
%              sf: shape function info
%            node: list of nodes for the element
%        N_of_kez: shape functions for a given element as a function of
%                  paramatric coordinate (ksi, eta, zeta)
%       dN_of_kez: shape derivative for a given element as a function of
%                  paramatric coordinate (ksi, eta, zeta)
%               N: shap functions at an integration point. (N = [N1; N2; N3:..])
%              dN: dNdX at an integration point. (dN = [dN1dX, dN1dY, dN1dZ; dN2dX, dN2dY, ...])
%            detJ: det(dXdkez)
%          detJxW: det(dXdkez) x weight
%             nne: number of nodes
%             nsd: number of spatial dimesions
%      curr_elmid: current element ID
%     curr_itg_id: current integration point ID
%              ST: list of 
% 
%  Methods:
%     set_an_element(coord, eid, eType, eOrder, iOrder) - set an finite element
%                         - coord : list of nodal points
%                         - eid   : element ID
%                         - eType : element type (Line, Triangle, Quadrilateral, ...)
%                         - eOrder: element order (0: linear, 1: quadratic, 2: cubic, ...), mostly support upto quadratic
%                         - iOrder: integration order (qudrature rule)
%     ElemBasis(itgID)    - update quadrature rule at the given integration
%                           point (itgID), it call compute_N_dN function
%     compute_N_dN(itgID) - compute shape functions and their derviative
%                           shape function at the integration point
%     x = getLocalEvalPt  - return a current integration point on local coordinates
%     x = getGlobalEvalPt - Return current integration point on global coordinates
%     update_shape_tensor - compute B matrix (obj.ST)
%
%  Example:
%     fe = FemLib;
%      v = 0;
%     node = [1; 3.2]; % line element 
%     fe.set_an_element(node, [], EnumElementType.Line);
%     for ip = 1: feV.quadRule.nint
%       fe.ElemBasis(ip);
%       v = v + feV.detJxW;
%     end
%     fprintf('line length = %e\n', v);
% 
%  Other m-files required: EnumElementType.m QuadratureRules.m
%  MAT-files required: ShapeFunctions.mat
% 
%  See also: BndFemLib
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 08-Oct-2020; Last revision:
  properties
    funcs; quadRule; sf;
    node;
    N_of_kez;
    dN_of_kez; 
    N; dN; detJ; detJxW
    nne; nsd;
    curr_elmid; curr_itg_id;
    ST;
  end
  methods
    function obj = FemLib
      obj.funcs = load('ShapeFunctions.mat');
    end
    function set_an_element(obj, coord, eid, eType, eOrder, iOrder)
      if(nargin < 5); eOrder = 0;  end % will set default
      if(nargin < 6); iOrder = -1; end % will set default quadrature rule
      
      obj.node = coord;
      [obj.nne, obj.nsd] = size(obj.node);
      obj.curr_elmid = eid;
     
      % compute dNdX
      obj.sf = ShapeFunction(eType, eOrder);
      obj.N_of_kez = obj.funcs.N.(sprintf('%s', eType)){eOrder+1};
      obj.dN_of_kez = obj.funcs.dN.(sprintf('%s', eType)){eOrder+1};      

      % assign quadrature rule
      obj.quadRule = QuadratureRules(eType, iOrder);
    end
    function ElemBasis(obj, itgID)
      %
      obj.curr_itg_id = itgID;
      J = obj.compute_N_dN(itgID);

      if(J < 0)
        return;
      end
      
      if(obj.nsd > obj.sf.nsd)
        vol2bnd.nne = obj.nne;
        vol2bnd.dN_of_kez = obj.dN_of_kez;
        vol2bnd.nodeIds = 1:obj.nne;
        [J, ~] = obj.compute_J_on_higher_dim(itgID, vol2bnd);
      end
      obj.detJ = J;
      obj.detJxW = J*obj.quadRule.weights(itgID);
    end
    function J = compute_N_dN(obj, itgID)
      iPt = obj.quadRule.points(itgID, :);
      obj.N = obj.N_of_kez(iPt(1), iPt(2), iPt(3));
      dNdkez = obj.dN_of_kez(iPt(1), iPt(2), iPt(3));

      X = zeros(obj.nne, 3);
      X(:, 1:obj.nsd) = obj.node;
      dXdkez = eye(3) + (X - obj.sf.node)'*dNdkez;
      J = det(dXdkez);
      if(J < 0)
        fprintf("ERROR: isoparametric J is %e (<=0)\n", J);
        return;
      end
      
      dXdkezI = dXdkez\eye(3);
      obj.dN = 1.0/J*dNdkez*dXdkezI;
    end
    function x = getLocalEvalPt(obj)
      % Return current integration point in local coordinates
      x = obj.quadRule.points(obj.curr_itg_id, :)';
    end
    function x = getGlobalEvalPt(obj)
      % Get current integration point in global coordinates
      x = obj.node'*obj.N;
    end
    function update_shape_tensor(obj)
      obj.ST = cell(obj.nne, obj.nsd);
      for a = 1: obj.nne
        for b = 1: obj.nsd
          st = zeros(obj.nsd, obj.nsd);
          for I = 1: obj.nsd
            if(I~=b)
              continue;
            end
            for J = 1: obj.nsd
              st(I,J) = obj.dN(a,J);
            end
          end
          obj.ST{a,b} = st;
        end
      end
    end
  end
end
