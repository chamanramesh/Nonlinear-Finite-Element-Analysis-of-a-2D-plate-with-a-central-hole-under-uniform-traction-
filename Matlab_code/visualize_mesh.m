% -----------------------------------------------------
% Author: Stefan Prueger, TU BAF, Nonlinear FEM, SS2023
% created: 2023-05-22
% -----------------------------------------------------
function [] = visualize_mesh (node_list,elem_list,n_elem,elem_type,figID)
%
figure(figID)
clf;
hold on
%
switch elem_type
  case 'Quad'
    ind=[1,5,2,6,3,7,4,8,1];
  case 'Tri'
    ind=[1,4,2,5,3,6,1];
end
for i=1:n_elem
  nodes_c_elem=elem_list(i,2:end);
  for k=2:length(ind)
    x_coord_segment=[node_list(nodes_c_elem(ind(k-1)),2),node_list(nodes_c_elem(ind(k)),2)];
    y_coord_segment=[node_list(nodes_c_elem(ind(k-1)),3),node_list(nodes_c_elem(ind(k)),3)];
    plot(x_coord_segment,y_coord_segment,'ko-')
  end
end
axis equal
xlabel('x'); ylabel('y');
end
