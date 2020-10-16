classdef ElementSide
% ElementSide - Define list of nodes composing element side and edge
%  Syntax:  elemSide = ElementSide;
%% 
%  Properties:
%     Side
%     Edge
%%
%  Methods:
%     argument for member functions denoted below:
%       - eType  : element type (Line, Triangle, Quadrilateral, ...)
%       - eOrder : element order (0: linear, 1: quadratic, 2: cubic, ...), mostly support upto quadratic
%       - sideId : element local side ID composing a volume element 
%     [nodeId, bnd_nne] = nodeIdVolToBnd(obj, eType, eOrder, sideId) 
%                         - get side boundary element at a given side (sideId)
%     node = sideNodeIDs(eType, eOrder) - get list of local element nodes composing boundary elements(face)
%     node = edgeNodeIDs(eType, eOrder) - get list of local element nodes composing element edges
%%
%  Example:
%     elemSide = ElementSide;
%     elemOrder = 0;
%     sideId = 1;
%     elemSide.Side.Triangle{elemOrder + 1}{sideId}
%     ans =
%          1     2
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 08-Oct-2020; Last revision:
  properties
    Side; Edge;
  end
  methods
    function obj = ElementSide
      % set side
      L0 = {1, 2};
      obj.Side.(sprintf('%s', EnumElementType.Line)) = {L0, L0, L0, L0};
      T0 = {[1,2], [2,3], [3,1]};
      T1 = {[1,2,4], [2,3,5], [3,1,6]};
      obj.Side.(sprintf('%s', EnumElementType.Triangle)) = {T0, T1};
      Q0 = {[1,2], [2,3], [3,4], [4,1]};
      Q1 = {[1,2,5], [2,3,6], [3,4,7], [4,1,8]};
      obj.Side.(sprintf('%s', EnumElementType.Quadrilateral)) = {Q0, Q1};
      Tet0 = {[1,3,2], [1,2,4], [2,3,4], [1,4,3]};
      Tet1 = {[1,3,2,7,6,5], [1,2,4,5,9,8], [2,3,4,6,10,9], [1,4,3,8,10,7]};
      obj.Side.(sprintf('%s', EnumElementType.Tetrahedron)) = {Tet0, Tet1};
      H0 = {[1,4,3,2], [5,6,7,8], [1,2,6,5], [8,7,3,4], [2,3,7,6], [5,8,4,1]};
      H1 = {[1,4,3,2,12,11,10,9], [5,6,7,8,13,14,15,16], ...
        [1,2,6,5,9,18,13,17], [8,7,3,4,15,19,11,20], ...
        [2,3,7,6,10,19,14,18], [5,8,4,1,16,20,12,17]};
      obj.Side.(sprintf('%s', EnumElementType.Hexahedron)) = {H0, H1};
      
      % set edge
      L0 = {1, 2}; L1 = {[1,3], [3,2]}; L2 = {[1,3], [3,4], [4,2]}; L3 = {[1,3], [3,4], [4,5], [5,2]};
      obj.Edge.(sprintf('%s', EnumElementType.Line)) = {L0, L1, L2, L3};
      T0 = {[1,2], [2,3], [3,1]};
      T1 = {[1,2,4], [2,3,5], [3,1,6]};
      obj.Edge.(sprintf('%s', EnumElementType.Triangle)) = {T0, T1};
      Q0 = {[1,2], [2,3], [3,4], [4,1]};
      Q1 = {[1,2,5], [2,3,6], [3,4,7], [4,1,8]};
      obj.Edge.(sprintf('%s', EnumElementType.Quadrilateral)) = {Q0, Q1};
      Tet0 = {[1,2], [2,3], [3,1], [1,4], [2,4], [3,4]};
      Tet1 = cell(1, numel(Tet0));
      for ia = 1: numel(Tet0)
        Tet1{(ia-1)*2+1} = [Tet0{ia}(1), 4+ia];
        Tet1{(ia-1)*2+2} = [4+ia, Tet0{ia}(2)];
      end
      obj.Edge.(sprintf('%s', EnumElementType.Tetrahedron)) = {Tet0, Tet1};
      H0 = {[1,2],[2,3],[3,4],[4,1],[5,6],[6,7],[7,8],[8,5],[1,5],[2,6],[3,7],[4,8]};
      H1 = cell(1, numel(H0));
      for ia = 1: numel(H0)
        H1{(ia-1)*2+1} = [H0{ia}(1), 8+ia];
        H1{(ia-1)*2+2} = [8+ia, H0{ia}(2)];
      end
      obj.Edge.(sprintf('%s', EnumElementType.Hexahedron)) = {H0, H1};
    end
    function [nodeId, bnd_nne] = nodeIdVolToBnd(obj, eType, eOrder, sideId)
      faces = obj.Side.(sprintf('%s', eType)){eOrder+1};
      nodeId = faces{sideId};
      bnd_nne = numel(nodeId);
    end
    function node = sideNodeIDs(obj, eType, eOrder)
      if(eType == EnumElementType.Line && eOrder == 0)
        node = cell2mat(obj.Side.(sprintf('%s', eType)){eOrder+1});
      else      
        node = cell2mat(obj.Side.(sprintf('%s', eType)){eOrder+1}(:));
      end
    end
    function node = edgeNodeIDs(obj, eType, eOrder)
      if(eType == EnumElementType.Line && eOrder == 0)
        node = cell2mat(obj.Edge.(sprintf('%s', eType)){eOrder+1});
      else
        node = cell2mat(obj.Edge.(sprintf('%s', eType)){eOrder+1}(:));
      end
    end    
  end
end