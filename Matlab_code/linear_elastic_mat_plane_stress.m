% -----------------------------------------------------
% Author: Stefan Prueger, TU BAF, Nonlinear FEM, SS2023
% created: 2023-05-08
% -----------------------------------------------------
function [stress,C_t] = linear_elastic_mat_plane_stress (strain,props)
  Emod=props(1); nu=props(2);
  C_t=Emod/(1-nu^2)*[[1,nu,0];[nu,1,0];[0,0,(1-nu)/2]];
  stress=C_t*strain;
end
