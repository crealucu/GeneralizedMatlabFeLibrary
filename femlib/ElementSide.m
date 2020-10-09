classdef ElementSide
% ElementSide - Define list of nodes composing element side
%  Syntax:  elemSide = ElementSide;
% 
%  Methods:
%  Properties:
%     Side
%  Example:
%     elemSide = ElementSide;
%     elemOrder = 0;
%     sideId = 1;
%     elemSide.Triangle{elemOrder + 1}{sideId}
%     ans =
%          1     2
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 08-Oct-2020; Last revision:
  properties
    Side;
  end
  methods
    function obj = ElementSide
      L0 = {1, 2};
      obj.Side.(sprintf('%s', EnumElementType.Line)) = {L0, L0, L0, L0};
      T0 = {[1,2], [2,3], [3,1]};
      T1 = {[1,2,4], [2,3,5], [3,1,6]};
      obj.Side.(sprintf('%s', EnumElementType.Triangle)) = {T0, T1};
      Q0 = {[1,2], [2,3], [3,4], [4,1]};
      Q1 = {[1,2,5], [2,3,6], [3,4,7], [4,1,8]};
      obj.Side.(sprintf('%s', EnumElementType.Quadrilateral)) = {Q0, Q1};
      Tet0 = {[1,3,2], [1,2,4], [2,3,4], [1,4,3]};
      Tet1 = {[1,2,3,7,6,5], [1,2,4,5,9,8], [2,3,4,6,10,9], [1,4,3,8,10,7]};      
      obj.Side.(sprintf('%s', EnumElementType.Tetrahedron)) = {Tet0, Tet1};
      H0 = {[1,4,3,2], [5,6,7,8], [1,2,6,5], [8,7,3,4], [2,3,7,6], [5,8,4,1]};
      H1 = {[1,4,3,2,12,11,10,9], [5,6,7,8,13,14,15,16], ...
            [1,2,6,5,9,18,13,17], [8,7,3,4,15,19,11,20], ...
            [2,3,7,6,10,19,14,18], [5,8,4,1,16,20,12,17]};
      obj.Side.(sprintf('%s', EnumElementType.Hexahedron)) = {H0, H1};
    end
  end
end