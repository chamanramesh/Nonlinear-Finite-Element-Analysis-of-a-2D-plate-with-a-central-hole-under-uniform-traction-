%%%%%%NLFEM Assignment 2023%%%%%%
%***************Code developed by********************
%         1.Chaman Ramesh [Matr.No:67771]           %
%       for       Variant 1                         %
%***************************************************%

%********main routine for the problem***************%

clear all; close all;

% material properties
Emod=210;         % Young's modulus in GPa
nu=0.3;          % Poisson's ratio
sigmaY0=0.2;    % initial yield stress in GPa
%sigmaY0=210;      % Linear case, or Emod=simgaY0
hard_iso=1.0;     % linear, isotropic hardening modulus in GPa
deltaY=0.0;      % increase in yield stress stemming from exponential hardening in GPa
h_sh=15;          % shape parameter for exponential hardening
hard_kin=0.0;     % linear kinematic hardening modulus in GPa
hard_funct='lin';
%
props=[Emod,nu,sigmaY0,hard_iso,deltaY,h_sh,hard_kin];
n_sdv_p_GP=7;
%
visualize_mesh_data=true;
%

% input file
foi='Mesh_Tri_Fine'; % 'Mesh_Tri_Coarse.mat' 'Mesh_Tri_Fine.mat' 'Mesh_Tri_vFine.mat'
elem_type='Tri';           % 'Tri'
%
% read mesh information
load(foi);
[n_nodes,n_dim]=size(node_list); n_dim=n_dim-1;
[n_elem,n_node_p_elem]=size(elem_list); n_node_p_elem=n_node_p_elem-1;
%
% visualization of initial mesh data
if visualize_mesh_data
    figID=1;
    visualize_mesh (node_list,elem_list,n_elem,elem_type,figID);
end

a=2.5; % a=0 uniform mesh,  radius of the circle 
node_coordinates = node_list(:, 2:end);
element_ids = elem_list(:, 2:end);

% Integration points and weights for triangular element 
L1 = 0.816847572980459;
L2 = 0.091576213509771;
L3 = 0.108103018168070;
L4 = 0.445948490915965;
W1 = 0.109951743655322 / 2.0;
W2 = 0.223381589678011 / 2.0;

i = [1,2,3,4,5,6];
si = [L1, L2, L2, L3, L4, L4];
ti = [L2, L1, L2, L4, L3, L4];
wi = [W1, W1, W1, W2, W2, W2];

% gauss point and weights for edge elements
se = [0.211324865405187, 0.788675134594813];
we = [0.5, 0.5]; 

% For disolaying displacement vs load step graph, at 50, 0 node to stroe
% the displacmenet
Edge_node_stored = [0];
Load_tau_stored = [0];

n_e = size(elem_list,1);
sdv_global = zeros(n_e*6,n_sdv_p_GP); % storing sdv at each gauss point for non linear calculation
sdv=zeros(1,n_sdv_p_GP); % sdv of current gauss point

sigma_infi = 0.195;
tau_n =  0; %current load 
delta_tau = sigma_infi/10; %Time step, to get more data points for plotting
tau_n =  tau_n + delta_tau;

u_e = zeros(12,1); % Elemental displacement
U_k = zeros(size(node_list,1)*2,1); % Global displacment
U_k_bck = zeros(size(node_list,1)*2,1); % Backup of gloabal displacment to store previous dispplacment values

new_matrix = zeros(size(elem_list,1)*6, 5);
node_elem = size(node_list,1)*2;
tolerance = 1e-6; 
k_max = 18; % Set time increments

%load increment loop
while(tau_n <= sigma_infi + 0.0001) % added tolerance to avoid termination of while loop like 0.195001 is considerd greater and is terminated

    delta_U_k = zeros(node_elem,1);
    k=1;
    U_k_bck = U_k;
    sdv_bck = sdv_global; % for each gauss point we have an indivdual arrary of sdv 
    new_matrix_bck = new_matrix;

    % Newton-Raphson iterations
    while(true)
        K_g = zeros(size(node_list,1)*2,size(node_list,1)*2);
        R_g = zeros(size(node_list,1)*2,1);
        sdv_counter = 1;

        % Essential boundary conditons
        rowsToRemove = node_sets{1, 2} * 2;
        colsToRemove = ((node_sets{4, 2} * 2) - 1);
        removed_rows_columns = [rowsToRemove,colsToRemove];

        % Element loop
        for e = 1:size(element_ids, 1)
            element = element_ids(e, :);
            element_coordinates = node_coordinates(element, :);

            A_e = zeros(12, node_elem);
            %Call function to use the Assembley matrix
            A_e = Assembly(element,A_e);
            A_e = sparse(A_e);

            u_e = A_e * U_k;

            Fe_int = zeros(12,1);
            Fe_ext = zeros(12,1);
            K = zeros(12, 12);

            %Call element routine
            [K,Fe_int,Fe_ext,matrix,sdv_global, sdv_counter] = element_routine(i,si,ti,wi,element_coordinates,u_e,node_coordinates,e,elem_sets,node_sets,props, sdv, hard_funct, element,elem_list,se,we,tau_n,sdv_global,sdv_counter, K, Fe_int, Fe_ext);

            %Compute K_global
            K = A_e' * K * A_e;
            K_g = K_g + K;

            %Compute R_global
            R = A_e' * (Fe_int - Fe_ext);
            R_g = R_g + R;

            % Store element stress values in the new matrix
            row_index = (e-1)*6 + 1;
            new_matrix(row_index:e*6, :) = matrix;
        end

        % Applying essential boundary conditions to remove respective rows and columns
        K_g(removed_rows_columns, :) = [];
        R_g(removed_rows_columns, :) = [];
        K_g(:, removed_rows_columns) = [];

        delta_U_k = -((K_g)\R_g);

        I = 1;
        % Assigning nodal values from delta_u to U_k(Global_U)
        for b = 1:node_elem 
            if any(removed_rows_columns == b) %at boundry, dont update, but at other point update u global, becuse at boundary alreday displacmenrt is 0
                continue
            else
                U_k(b) = U_k(b) + delta_U_k(I);
                I = I + 1; % incrementing counter for delta_u
            end
        end

        % Check convergence criteria

        %fprintf("CurrLoad: %f | delta_tau: %f | k: %d | norm(delta_u): %f\n", tau_n,delta_tau, k, norm(delta_U_k));
        if (norm(R_g)<=tolerance && norm(delta_U_k)<=tolerance)
            fprintf("Newton Raphson converged...\n")
            node_at_corner = node_list(((node_list(:,2)==50) & (node_list(:,3)==0)),1);
            Edge_node_stored = [Edge_node_stored; U_k(node_at_corner*2-1)]; % extract global displacment for the corresponding node number
            Load_tau_stored = [Load_tau_stored; (tau_n/sigma_infi)];
            tau_n = tau_n + delta_tau;
            break;
        elseif k>=k_max
            fprintf("Newton Raphson not converged...")
            tau_n = tau_n - delta_tau;
            delta_tau = delta_tau/10;
            tau_n = tau_n + delta_tau;
            U_k = U_k_bck;
            sdv_global = sdv_bck;
            new_matrix = new_matrix_bck;
            break;
        else
            k = k+1;
            continue;
        end
    end
end


%------------------------------------Plotting-------------------------------------------------------%

%Extracting the stored stress values
tolerance = 1;
filtered_indices_s = (new_matrix(:, end-1) - min(new_matrix(:,end-1))) < tolerance;
matrix_90 = [new_matrix(filtered_indices_s, 1:3), new_matrix(filtered_indices_s, end-1:end)]; 

filtered_indices_t = (new_matrix(:, end) - min(new_matrix(:,end))) < tolerance;
matrix_0 = [new_matrix(filtered_indices_t, 1:3), new_matrix(filtered_indices_t, end-1:end)];

fs = new_matrix(:, end-1);
ft = new_matrix(:, end);
filtered_indices_st = abs(fs - ft) < tolerance;
matrix_45 = [new_matrix(filtered_indices_st, 1:3), new_matrix(filtered_indices_st, end-1:end)];


%Converting from cartesian to polar coordinates
% Convert s values to polar coordinates (r_s, theta_s)
[r_s, theta_s] = cart2pol(matrix_90(:, 4), matrix_90(:, 5));
matrix_s0_polar = [matrix_90(:, 1:3), r_s, theta_s];

% Convert t values to polar coordinates (r_t, theta_t)
[r_t, theta_t] = cart2pol(matrix_0(:, 4), matrix_0(:, 5));
filtered_matrix_t_polar = [matrix_0(:, 1:3), r_t, theta_t];

% Convert difference values (s-t) to polar coordinates (r_st, theta_st)
[r_st, theta_st] = cart2pol(matrix_45(:, 4), matrix_45(:, 5));
filtered_matrix_st_polar = [matrix_45(:, 1:3), r_st, theta_st];


for i = 1:size(matrix_s0_polar,1)

    theta = matrix_s0_polar(i, 4);
    r = matrix_s0_polar(i, 5);

    sigma_11 = matrix_s0_polar(i, 1);
    sigma_12 = matrix_s0_polar(i, 2);
    sigma_22 = matrix_s0_polar(i, 3);

    T = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    sigma_s = [sigma_11, sigma_12; sigma_12, sigma_22];

    sigma_result_s = T * sigma_s * T';

    transformed_stress_components_s_f(i, :) = [sigma_result_s(1, 1), sigma_result_s(1, 2), sigma_result_s(2, 2), r]; % 90 Degree fem plot


    % Anaytical stress calculations
    term1 = (1 - (a / r)^2);
    term2 = (1 - 4*(a / r)^2 + 3*(a / r)^4) * cos(2*theta);
    sigma_rr = sigma_infi / 2 * (term1 + term2);

    term1 = (1 + (a / r)^2);
    term2 = (1 + 3*(a / r)^4) * cos(2*theta);
    sigma_thetatheta = sigma_infi / 2 * (term1 - term2);

    term1 = (1 + 2*(a / r)^2);
    term2 = (3*(a / r)^4);
    sigma_rtheta = - sigma_infi / 2 *((term1 - term2))*sin(2*theta);

    transformed_stress_components_s_a(i, :) = [sigma_rr, sigma_rtheta, sigma_thetatheta, r]; % 90 Degree analytical plot
end

for j = 1:size(filtered_matrix_t_polar,1)

    theta = filtered_matrix_t_polar(j, 4);
    r = filtered_matrix_t_polar(j, 5);

    sigma_11 = filtered_matrix_t_polar(j, 1);
    sigma_12 = filtered_matrix_t_polar(j, 2);
    sigma_22 = filtered_matrix_t_polar(j, 3);

    T = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    sigma_s = [sigma_11, sigma_12; sigma_12, sigma_22];

    sigma_result_t = T * sigma_s * T';
    transformed_stress_components_t_f(j, :) = [sigma_result_t(1, 1), sigma_result_t(1, 2), sigma_result_t(2, 2), r]; % 0 Degree fem plot


    % Anaytical stress calculations
    term1 = (1 - (a / r)^2);
    term2 = (1 - 4*(a / r)^2 + 3*(a / r)^4) * cos(2*theta);
    sigma_rr = sigma_infi / 2 * (term1 + term2);

    term1 = (1 + (a / r)^2);
    term2 = (1 + 3*(a / r)^4) * cos(2*theta);
    sigma_thetatheta = sigma_infi / 2 * (term1 - term2);

    term1 = (1 + 2*(a / r)^2);
    term2 = (3*(a / r)^4);
    sigma_rtheta = - sigma_infi / 2 *((term1 - term2))*sin(2*theta);

    transformed_stress_components_t_a(j, :) = [sigma_rr, sigma_rtheta, sigma_thetatheta, r]; % 0 Degree analytical plot
end

for k = 1:size(filtered_matrix_st_polar,1)

    theta = filtered_matrix_st_polar(k, 4);
    r = filtered_matrix_st_polar(k, 5);

    sigma_11 = filtered_matrix_st_polar(k, 1);
    sigma_12 = filtered_matrix_st_polar(k, 2);
    sigma_22 = filtered_matrix_st_polar(k, 3);

    T = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    sigma_s = [sigma_11, sigma_12; sigma_12, sigma_22];

    sigma_result_st = T * sigma_s * T';
    transformed_stress_components_st_f(k, :) = [sigma_result_st(1, 1), sigma_result_st(1, 2), sigma_result_st(2, 2), r]; % 45 Degree fem plot


    % Anaytical stress calculations
    term1 = (1 - (a / r)^2);
    term2 = (1 - 4*(a / r)^2 + 3*(a / r)^4) * cos(2*theta);
    sigma_rr = sigma_infi / 2 * (term1 + term2);

    term1 = (1 + (a / r)^2);
    term2 = (1 + 3*(a / r)^4) * cos(2*theta);
    sigma_thetatheta = sigma_infi / 2 * (term1 - term2);

    term1 = (1 + 2*(a / r)^2);
    term2 = (3*(a / r)^4);
    sigma_rtheta = - sigma_infi / 2 *((term1 - term2))*sin(2*theta);

    transformed_stress_components_st_a(k, :) = [sigma_rr, sigma_rtheta, sigma_thetatheta, r]; % 45 Degree fem plot
end
%--------------------------------------------------------------------------------------------------------------------------%
%Plotting graph


%Theta = 90 Degree
% Extract and sort fem stresses and analytcal stresses in ascending order based on the last column (r)
sorted_data = sortrows(transformed_stress_components_s_f, size(transformed_stress_components_s_f, 2));
sigma11 = sorted_data(:, 1);
sigma12 = sorted_data(:, 2);
sigma22 = sorted_data(:, 3);
r = sorted_data(:, 4);

sorted_analytical_data = sortrows(transformed_stress_components_s_a, size(transformed_stress_components_s_a, 2));
analytical_sigma_rr = sorted_analytical_data(:, 1);
analytical_sigma_rtheta = sorted_analytical_data(:, 2);
analytical_sigma_thetatheta = sorted_analytical_data(:, 3);
analytical_r = sorted_analytical_data(:, 4);

% Plot for Stress(Gpa) vs Radial distance(mm) for theta = 90 Degree
figure;
hold on;
plot(r, sigma11, 'Color', 'r', 'LineWidth', 1);
plot(r, sigma12, 'Color', 'g', 'LineWidth', 1);
plot(r, sigma22, 'Color', 'b', 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rr, 'Color', 'm','marker', 'square', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rtheta, 'Color', 'y','marker', '*', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_thetatheta, 'Color', 'k','marker', 'diamond', 'MarkerSize', 6, 'LineWidth', 1);
grid on;
hold off;
title('90 Degree');

xlabel('Radial Distance (mm)');
ylabel('Stress (Gpa)');
legend('σrr FEM', 'σrθ FEM', 'σθθ FEM', 'σrr Analytical', 'σrθ Analytical', 'σθθ Analytical');


%Theta = 0 Degree
% Extract and sort fem stresses and analytcal stresses in ascending order based on the last column (r)
sorted_data = sortrows(transformed_stress_components_t_f, size(transformed_stress_components_t_f, 2));
sigma11 = sorted_data(:, 1);
sigma12 = sorted_data(:, 2);
sigma22 = sorted_data(:, 3);
r = sorted_data(:, 4);

sorted_analytical_data = sortrows(transformed_stress_components_t_a, size(transformed_stress_components_t_a, 2));
analytical_sigma_rr = sorted_analytical_data(:, 1);
analytical_sigma_rtheta = sorted_analytical_data(:, 2);
analytical_sigma_thetatheta = sorted_analytical_data(:, 3);
analytical_r = sorted_analytical_data(:, 4);

% Plot for Stress(Gpa) vs Radial distance(mm) for theta = 0 Degree.
figure;
hold on;
plot(r, sigma11, 'Color', 'r', 'LineWidth', 1);
plot(r, sigma12, 'Color', 'g', 'LineWidth', 1);
plot(r, sigma22, 'Color', 'b', 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rr, 'Color', 'm','marker', 'square', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rtheta, 'Color', 'y','marker', '*', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_thetatheta, 'Color', 'k','marker', 'diamond', 'MarkerSize', 6, 'LineWidth', 1);
grid on;
hold off;
title('0 Degree');
xlabel('Radial Distance (mm)');
ylabel('Stress (Gpa)');
legend('σrr FEM', 'σrθ FEM', 'σθθ FEM', 'σrr Analytical', 'σrθ Analytical', 'σθθ Analytical', 'Location', 'southeast');


%Theta = 45 Degree
% Extract and sort fem stresses and analytcal stresses in ascending order based on the last column (r)
sorted_data = sortrows(transformed_stress_components_st_f, size(transformed_stress_components_st_f, 2));
sigma11 = sorted_data(:, 1);
sigma12 = sorted_data(:, 2);
sigma22 = sorted_data(:, 3);
r = sorted_data(:, 4);

sorted_analytical_data = sortrows(transformed_stress_components_st_a, size(transformed_stress_components_st_a, 2));
analytical_sigma_rr = sorted_analytical_data(:, 1);
analytical_sigma_rtheta = sorted_analytical_data(:, 2);
analytical_sigma_thetatheta = sorted_analytical_data(:, 3);
analytical_r = sorted_analytical_data(:, 4);

% Plot for Stress(Gpa) vs Radial distance(mm) for theta = 45 Degree.
figure;
hold on;
plot(r, sigma11, 'Color', 'r', 'LineWidth', 1);
plot(r, sigma12, 'Color', 'g', 'LineWidth', 1);
plot(r, sigma22, 'Color', 'b', 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rr, 'Color', 'm','marker', 'square', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_rtheta, 'Color', 'y', 'marker', '*', 'MarkerSize', 6, 'LineWidth', 1);
plot(analytical_r, analytical_sigma_thetatheta, 'Color', 'k','marker', 'diamond', 'MarkerSize', 6, 'LineWidth', 1);
grid on;
hold off;
title('45 Degree');
xlabel('Radial Distance (mm)');
ylabel('Stress (Gpa)');
legend('σrr FEM', 'σrθ FEM', 'σθθ FEM', 'σrr Analytical', 'σrθ Analytical', 'σθθ Analytical');


% Plot for Displacmenet vs Time.
figure;
plot(Load_tau_stored, Edge_node_stored, 'Color', '#0072BD','Marker','pentagram');
grid on;
title('Displacement vs Time');
xlabel('Loadstep Time');
ylabel('Displacement (mm)');
legend('Coarse Mesh', 'Location', 'northwest'); % Change the label according to the type of meshing that is being used to run the program.