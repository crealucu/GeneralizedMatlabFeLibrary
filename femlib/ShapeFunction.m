classdef ShapeFunction < handle
% ShapeFunction - class for finite element shape functions
% 
%  Syntax:  sf = ShapeFunction(elementType, elementOrder)
% 
%  properties:
%       u_of_X: polynominal equation represents a finite element shape
%               function w.r.t X = {ksi, eta, zet} coordinate
%         node: list of nodal coordinate where shape functions are 1
%            N: shape function as a function of X = {ksi, eta, zeta}
%           dN: derivative of shape function w.r.t X = {ksi, eta, zeta}
%            X: {ksi, eta, zeta}
%     elemType: type of element (Point, Line, Triangle, Quadrilateral, ...)
%        order: order of finite element
%          nsd: number of spatial dimensions
% 
%  methods:
%     ShapeFunction(eType, eOrder) - class constructor. set element type and order of the element
%     compute_shape_function       - compute shape functions
%     plot_element                 - plot isoparametric element
%
%  Example:
%     sf = ShapeFunction(EnumElementType.Triangle, 0);
%     sf.compute_shape_function;
%     sf = 
%     ShapeFunction with properties:
%       u_of_X: @(c,ksi,eta,zeta)c(1)+c(2)*ksi+c(3)*eta
%         node: [3×3 double]
%            N: [3×1 sym]
%           dN: [3×3 sym]
%            X: [1×3 sym]
%     elemType: Triangle
%        order: 0
%          nsd: 2
%
%  Other m-files required: EnumElementType.m
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 08-Oct-2020; Last revision:  
  properties
    u_of_X; node; N; dN;
    X;
    elemType;
    order;
    nsd;
  end
  methods
    function obj = ShapeFunction(eType, eOrder)
      syms ksi
      syms eta
      syms zeta
      obj.X = [ksi, eta, zeta];
      obj.elemType = eType;
      if(nargin<2)
        eOrder = 0;
      end
      obj.order = eOrder;
      switch(eType)
        case EnumElementType.Point
          [obj.u_of_X, obj.node, obj.nsd] = Point();
        case EnumElementType.Line
          [obj.u_of_X, obj.node, obj.nsd] = Line(eOrder);
        case EnumElementType.Triangle
          [obj.u_of_X, obj.node, obj.nsd] = Triangle(eOrder);
        case EnumElementType.Quadrilateral
          [obj.u_of_X, obj.node, obj.nsd] = Quadrilateral(eOrder);
        case EnumElementType.Tetrahedron
          [obj.u_of_X, obj.node, obj.nsd] = Tetrahedron(eOrder);
        case EnumElementType.Hexahedron
          [obj.u_of_X, obj.node, obj.nsd] = Hexahedron(eOrder);
      end
    end
    function compute_shape_function(obj)
      nodeno = size(obj.node, 1);
      
      % compute A in equation u = A*c
      c = sym('c', [nodeno, 1]);
      A = zeros(nodeno, nodeno);
      for ia = 1: nodeno
        x = zeros(3,1);
        for ib = 1: 3
          x(ib) = obj.node(ia, ib);
        end
        tmp = obj.u_of_X(c, x(1), x(2), x(3));
        for ib = 1: nodeno
          A(ia, ib) = diff(tmp, c(ib));
        end
      end
      
      % compute c = M*u; u = 1 on its node otherwize 0
      M = A\eye(nodeno);
      
      % compute shape function
      % N_i = c(:, i)'*ksi_eta_zeta
      % where ksi_eta_zeta = [1; ksi; eta; zeta; ...]
      % and computed by ksi_eta_zeta(j) = u_of_X([0, 0, ..., 1(jth), ...], ksi, eta, zeta)
      obj.N = c(1)*zeros(nodeno, 1); % simbolic zero
      for ia = 1: nodeno
        d = c(1)*0;
        for ib = 1: nodeno
          x = zeros(nodeno, 1);
          x(ib) = 1;
          d = d + M(ib, ia)*obj.u_of_X(x, obj.X(1), obj.X(2), obj.X(3));
        end
        obj.N(ia) = simplify(d);
      end
      
      obj.dN = c(1)*zeros(nodeno, 3);
      for ia = 1: nodeno
        for ib = 1: 3
          obj.dN(ia, ib) = diff(obj.N(ia), obj.X(ib));
        end
      end
    end
    function plot_element(obj)
      n = obj.nsd;
      if(n < 2)
        n = 2;
      end
      e = [];
      switch(obj.elemType)
        case EnumElementType.Point
          e = [1,1];
        case EnumElementType.Line
          e = [1,2];
        case EnumElementType.Triangle
          e = [1,2; 2,3; 3,1];
        case EnumElementType.Quadrilateral
          e = [1,2; 2,3; 3,4; 4,1];
        case EnumElementType.Tetrahedron
          e = [1,2; 2,3; 3,1; 1,4; 2,4; 3,4];
        case EnumElementType.Hexahedron
          e = [1,2; 2,3; 3,4; 4,1; 5,6; 6,7; 7,8; 8,5;
               1,5; 2,6; 3,7; 4,8];
      end
      
      fd = FIGDATA;
      for ia = 1: size(e, 1)
        fd.add_xy([obj.node(e(ia, 1), 1:n); obj.node(e(ia, 2), 1:n)], '-k')
      end
      fd.add_xy(obj.node(:, 1:n), 'ko');
      
      fd.myplots;
      axis off;
    end
  end
end

function [u_of_X, node, nsd] = Point()
  nsd = 0;  
  u_of_X = @(c, ksi, eta, zeta) c(1);
  node  = [0,0,0];
end

function [u_of_X, node, nsd] = Line(order)
  nsd = 1;
  switch(order)
    case 0
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi;
      node  = [-1,0,0; 1,0,0];
    case 1
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*ksi*ksi;
      node  = [-1,0,0; 1,0,0; 0,0,0];
    case 2
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*ksi*ksi + c(4)*ksi*ksi*ksi;
      node  = [-1,0,0; 1,0,0; -1/3,0,0; 1/3,0,0];
    otherwise
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi;
      node  = [-1,0,0; 1,0,0];
  end
end

function [u_of_X, node, nsd] = Triangle(order)
  nsd = 2;
  switch(order)
    case 0
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta;
      node  = [0,0,0; 1,0,0; 0,1,0];
    case 1
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*ksi*ksi + c(5)*ksi*eta + c(6)*eta*eta;
      node  = [0,0,0; 1,0,0; 0,1,0; 1/2,0,0; 1/2,1/2,0; 0,1/2,0];
    otherwise
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta;
      node  = [0,0,0; 1,0,0; 0,1,0];
  end
end

function [u_of_X, node, nsd] = Quadrilateral(order)
  nsd = 2;
  switch(order)
    case 0
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*ksi*eta;
      node  = [-1,-1,0; 1,-1,0; 1,1,0; -1,1,0];
    case 1
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*ksi*eta + ...
                                    c(5)*ksi*ksi + c(6)*eta*eta + c(7)*ksi*ksi*eta + c(8)*eta*eta*ksi;
      node  = [-1,-1,0; 1,-1,0; 1,1,0; -1,1,0;
                0,-1,0; 1, 0,0; 0,1,0; -1,0,0];
    otherwise
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*ksi*eta;
      node  = [-1,-1,0; 1,-1,0; 1,1,0; -1,1,0];
  end
end

function [u_of_X, node, nsd] = Tetrahedron(order)
  nsd = 3;
  switch(order)
    case 0
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta;
      node  = [0,0,0; 1,0,0; 0,1,0; 0, 0, 1];
    case 1
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta + ...
                                    c(5)*ksi*eta + c(6)*eta*zeta + c(7)*zeta*ksi + ...
                                    c(8)*ksi*ksi + c(9)*eta*eta  + c(10)*zeta*zeta;
      node  = [0,0,0; 1,0,0; 0,1,0; 0,0,1; 1/2,0,0; 1/2,1/2,0; 0,1/2,0; 0,0,1/2; 1/2,0,1/2; 0,1/2,1/2];
    otherwise
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta;
      node  = [0,0,0; 1,0,0; 0,1,0; 0, 0, 1];
  end
end

function [u_of_X, node, nsd] = Hexahedron(order)
  nsd = 3;
  switch(order)
    case 0
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta + ...
                                    c(5)*ksi*eta + c(6)*eta*zeta + c(7)*zeta*ksi + c(8)*ksi*eta*zeta;
      node  = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];
    case 1
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta + ...
                                    c(5)*ksi*eta + c(6)*eta*zeta + c(7)*zeta*ksi + c(8)*ksi*eta*zeta + ...
                                    c(9)*ksi*ksi + c(10)*eta*eta + c(11)*zeta*zeta + ...
                                    c(12)*ksi*ksi*eta + c(13)*ksi*ksi*zeta + c(14)*ksi*ksi*eta*zeta + ...
                                    c(15)*eta*eta*ksi + c(16)*eta*eta*zeta + c(17)*eta*eta*zeta*ksi + ...
                                    c(18)*zeta*zeta*ksi + c(19)*zeta*zeta*eta + c(20)*zeta*zeta*ksi*eta;
                                    
      node  = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1;
               0,-1,-1; 1,0,-1; 0,1,-1; -1,0,-1; 0,-1, 1; 1,0,1; 0,1,1; -1,0,1;
               -1,-1,0; 1,-1,0; 1,1,0; -1,1,0];
    otherwise
      u_of_X = @(c, ksi, eta, zeta) c(1) + c(2)*ksi + c(3)*eta + c(4)*zeta + ...
                                    c(5)*ksi*eta + c(6)*eta*zeta + c(7)*zeta*ksi + c(8)*ksi*eta*zeta;
      node  = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1];
  end
end