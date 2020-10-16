classdef BoundaryElement < handle
% BoundaryElement - class for constructing mesh boundary
%
%  Syntax: bndElem   = BoundaryElement
%
%  Properties:
%                elem: list of boundary element, currently assumed homogeneous element
%     LocalElemSideId: list of element side IDs for elem
%            BndToVol: list of volume element IDs for elem
%       bndElemSideId: list of user defined element side IDs
%       BndNodeSideId: list of node's user defined side IDs, each node can
%                      have multiple IDs or no IDs
%      MarkedNodeList: list of user define side IDs on nodes
%                      = {[list of nodes for side ID 1], [for side ID 2], [for side ID 3], ...}
%              elemno: number of boundary elements
%                 nne: number of nodes in each boundary elements
%         elemSideIDs: element local side IDs are defined
% 
%  Methods:
%     ConstructBoundaryElement(vol) - construct boudnary elements from volume elements
%                                   - vol: volume mesh class (Mesh.m)
%     MarkBndWithXmaxXmin(node)     - assign user defined side IDs automatically based on nodal position
%                                   - this function will construct 
%                                   - MarkedNodeList = {[node on min. X], [node on max. X], ...
%                                                       [node on min. Y], [node on max. Y], ...
%                                                       [node on min. Z], [node on max. Z]};
%     SetBoundaryNodeSideId(nodeno) - construct list of node's user defined side IDs (MarkedNodeList)
%     SetBoundaryElemSideId()       - construct list of user defined element side IDs (MarkedNodeList)
% 
%  Other m-files required: femlib
%  See also: Mesh, MeshGenerator
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 16-Oct-2020; Last revision:
  properties
    elem;
    LocalElemSideId; BndToVol;
    bndElemSideId; BndNodeSideId; MarkedNodeList;
    elemno; nne;
    elemSideIDs;
  end
  methods
    function obj = BoundaryElement
      obj.elemno = 0;
      obj.LocalElemSideId = [];
      obj.BndToVol = [];
      obj.elem    = [];
      obj.bndElemSideId = [];
      obj.BndNodeSideId = {};
      obj.MarkedNodeList = {};
      obj.elemSideIDs = ElementSide;
    end
    function ConstructBoundaryElement(obj, vol)
      %
      eid = (1:vol.elemno)';
      faceIDs = obj.elemSideIDs.sideNodeIDs(vol.elemType, vol.elemOrder);
      [faceno, face_nne] = size(faceIDs);
      
      face = zeros(vol.elemno*faceno, face_nne + 2);
      
      % list of face (boundary) = [list of 1st faces, element IDs, 1st local side IDs
      %                            list of 2nd faces, element IDs, 2nd local side IDs
      %                            :,                 :,           :                 ]
      % a face = [node 1, node 2, ..., node n];
      for ia = 0: faceno-1
        face(vol.elemno*ia + eid, :) = [vol.elem(:, faceIDs(ia+1, :)), eid, ones(vol.elemno, 1) + ia];
      end
      
      % identify duplications of the faces
      [oface, ~] = sort(face(:, 1:end-2), 2);
      [~, ix, jx] = unique(oface, 'rows');
       
      % unique faces are boundary
      N = histcounts(jx, 1:(size(oface, 1)+1));
      
      qx = (N==1);
      obj.LocalElemSideId = face(ix(qx), end);
      obj.BndToVol        = face(ix(qx), end-1);
      obj.elem = face(ix(qx), 1:end-2);
      [obj.elemno, obj.nne] = size(obj.elem);

    end    
    function MarkBndWithXmaxXmin(obj, node)
      % compute maxima, minima
      xmin = min(node);
      xmax = max(node);
      dim = numel(xmin);
      nodeno = size(node, 1);
      
      % initialize list
      obj.MarkedNodeList = cell(1,dim*2);
      id = (1:nodeno)';
      
      % mark nodes for each coordinate: X, Y, and Z
      for a = 1: dim
        obj.MarkedNodeList{2*(a-1) + 1} = id(node(:, a)==xmin(a));
        obj.MarkedNodeList{2*(a-1) + 2} = id(node(:, a)==xmax(a));
      end
      
      % construct assign side IDs for boundary elements and nodes
      obj.SetBoundaryNodeSideId(nodeno);
      obj.SetBoundaryElemSideId;
    end
    function SetBoundaryNodeSideId(obj, nodeno)
      obj.BndNodeSideId = cell(nodeno, 1);
      m = numel(obj.MarkedNodeList);
      for ia = 1: m
        N = obj.MarkedNodeList{ia};
        n = numel(N);
        for ib = 1: n
          obj.BndNodeSideId{N(ib)}(end+1) = ia;
        end
      end
    end    
    function SetBoundaryElemSideId(obj)
      %
      % find maximum number of boundary that can touch one node
      max_face_no = 0;
      for ia = 1: numel(obj.BndNodeSideId)
        max_face_no = max(max_face_no, numel(obj.BndNodeSideId{ia}));
      end
      
      % construct boundary element side IDs (user defined) list
      obj.bndElemSideId = zeros(obj.elemno, 1) - 1;

      for ia = 1: obj.elemno
        tmp_face_ids = zeros(1, obj.nne*max_face_no) - 1; % the size is initially unkown
        % collect all node user defined side IDs for this boundary element
        for ib = 1: obj.nne
          iDs = obj.BndNodeSideId{obj.elem(ia, ib), 1};
          idno = numel(iDs);
          if(idno == 0); continue; end
          tmp_face_ids(end+(1:idno)) = iDs;
        end
        
        tmp_face_ids(tmp_face_ids<0) = [];
        if(numel(tmp_face_ids)<0); continue; end
        
        % count number of side IDs and if one of side IDs is counted same
        % as nne, the side IDs is for the boundary element
        [~, ix, jx] = unique(tmp_face_ids);        
        N = histcounts(jx, 1:(numel(tmp_face_ids)+1));
        qx = (N==obj.nne);
        obj.bndElemSideId(ia) = tmp_face_ids(ix(qx));
      end
    end
  end
end