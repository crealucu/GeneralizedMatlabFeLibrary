function setup_GMFeL(uninstall)
% setup_GMFeL - install GFEM library
%
%  Syntax:  setup_GMFeL(uninstall)
%  inputs:  \param[in] uninstall - if yes, uninstall the library
%  Other m-files required: matlab_functions library
%
% Author: Sangmin Lee, Ph.D.
% email: crealucu@gmail.com
% 04-Oct-2020; Last revision:
  if(nargin<1)
    uninstall = false;
  end
  
  cdir = pwd;
  subDirsTasks = {% dirname {install task1, task2, ...} {uninstall task1, task2, ...}
                  'femlib', {'generate_shape_functions', 'set_lib_path'} {'delete ShapeFunctions.mat', 'set_lib_path(true)'}
                  };
  for ia = 1: size(subDirsTasks)
    if(uninstall)
      fprintf('Uninstall for %s ...\n', subDirsTasks{ia, 1});
    else
      fprintf('Install for %s ...\n', subDirsTasks{ia, 1});
    end
    cd(subDirsTasks{ia, 1});
    for ib = 1: numel(subDirsTasks{ia, 2})
      if(uninstall)
        eval(subDirsTasks{ia, 3}{ib});
      else
        eval(subDirsTasks{ia, 2}{ib});
      end
    end
    cd(cdir);
  end
end
