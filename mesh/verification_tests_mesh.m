function verification_tests_mesh
% verification_tests_mesh - perform verification tests for mesh
% a box mesh for each element type is generated and integrated over volume
% and a face.
% 
%  Syntax: verification_tests_mesh;
% 
%  Outputs:
%     test results will be printed
%
%  Other m-files required: femlib, Mesh.m, MeshGenerator.m, BoundaryElement.m
% 
%  See also: verification_tests_boundary_integration, verification_tests_volume_integration

% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 16-Oct-2020; Last revision:
%
  tol = 1.0e-12;
  xmin = [0.0, 0.0, 0.0];
  xmax = [2.0, 1.0, 1.0];
  nx   = [4,4,4];

  v0 = 2.0;
  s0 = 2.0;
  % problem set
  %     dof, isTet, elemType
  info = {1, false, 0
          1, false, 1
          2, false, 0
          2, false, 1
          2, true,  0
          2, true,  1
          3, false, 0
          3, false, 1
          3, true,  0
          3, true,  1          
          };
  for ia = 1: size(info, 1)
    mesh = MeshGenerator;
  
    dof    = info{ia, 1};
    isTet  = info{ia, 2};
    eOrder = info{ia, 3};

    mesh.generate_box(xmin, xmax, nx, dof, isTet, eOrder);
    mesh.bndElem.MarkBndWithXmaxXmin(mesh.node);
    
    fprintf('test for %s(order = %d): ', mesh.elemType, mesh.elemOrder);    
    s = compute_area(mesh);
    v = compute_volume(mesh);
    if(dof == 1)
      err_s = abs(s);
    else
      err_s = sqrt((s - s0)^2);
    end
    err_v = sqrt((v - v0)^2);

    if((err_s + err_v) < tol)
      cprintf('blue', 'Pass\n');
    else
      cprintf('red', ': Failed\n');
    end
    fprintf('volume = %.17f, err = %e\n', v, err_v);    
    fprintf('area   = %.17f, err = %e\n', s, err_s);
  end
end

function v = compute_volume(mesh)
  v = 0;
  for eid = 1: mesh.elemno
    fe = FemLib;
    node = mesh.node(mesh.elem(eid, :), :);
    fe.set_an_element(node, eid, mesh.elemType, mesh.elemOrder);
    for ip = 1: fe.quadRule.nint
      fe.ElemBasis(ip);
      v = v + fe.detJxW;
    end
  end 
end

function s = compute_area(mesh)
  s = 0;
  for fid = 1: mesh.bndElem.elemno
    eid    = mesh.bndElem.BndToVol(fid);
    sideId = mesh.bndElem.LocalElemSideId(fid);
    if(mesh.bndElem.bndElemSideId(fid) ~= 3); continue; end
    fe = BndFemLib;
    node = mesh.node(mesh.elem(eid, :), :);
    fe.set_an_element(node, eid, mesh.elemType, sideId, mesh.elemOrder);
    for ip = 1: fe.quadRule.nint
      fe.ElemBasis(ip);
      s = s + fe.detJxW;
    end
  end 
end