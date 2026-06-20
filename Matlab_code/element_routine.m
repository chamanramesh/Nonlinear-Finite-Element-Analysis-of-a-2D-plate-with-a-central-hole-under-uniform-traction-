%%%%%%NLFEM Assignment 2023%%%%%%
%***************Code developed by********************
%         1.Chaman Ramesh [Matr.No:67771]           %
%       for       Variant 1                         %
%***************************************************%

%********Element Routine for the problem************%


function [K, Fe_int, Fe_ext, matrix, sdv_global, sdv_counter] = element_routine(i, si, ti, wi, element_coordinates, u_e, node_coordinates, e, elem_sets, node_sets, props, sdv, hard_funct, element, elem_list, se, we, tau_n, sdv_global, sdv_counter, K, Fe_int, Fe_ext)

% Store stress values and s and t values
matrix = zeros(6, 5);

for j = i
    s = si(j);
    t = ti(j);
    w = wi(j);

    % Calculate shape functions N1 to N6
    N1 = 2*s^2 + 2*t^2 + 4*s*t - 3*t - 3*s + 1;
    N2 = 2*s^2 - s;
    N3 = 2*t^2 - t;
    N4 = 4*s - 4*s^2 - 4*s*t;
    N5 = 4*s*t;
    N6 = 4*t - 4*s*t - 4*t^2;

    N = [N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0;
        0, N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6];

    % Calculate derivatives of shape functions with respect to s and t
    dN1_ds = 4*s + 4*t - 3;
    dN2_ds = 4*s - 1;
    dN3_ds = 0;
    dN4_ds = 4 - 8*s - 4*t;
    dN5_ds = 4*t;
    dN6_ds = -4*t;

    dN1_dt = 4*s + 4*t - 3;
    dN2_dt = 0;
    dN3_dt = 4*t - 1;
    dN4_dt = -4*s;
    dN5_dt = 4*s;
    dN6_dt = 4 - 4*s - 8*t;

    % Compute the dN_dE matrix
    dN_dE = [dN1_ds, dN2_ds, dN3_ds, dN4_ds, dN5_ds, dN6_ds;
        dN1_dt, dN2_dt, dN3_dt, dN4_dt, dN5_dt, dN6_dt];

    % Compute the Jacobian matrix J
    J = dN_dE * element_coordinates;

    % Compute the B matrix using the inverse of J
    B = inv(J) * dN_dE;

    % Compute B_matrix
    B_MAT = [B(1, 1), 0, B(1, 2), 0, B(1, 3), 0, B(1, 4), 0, B(1, 5), 0, B(1, 6), 0;
        0, B(2, 1), 0, B(2, 2), 0, B(2, 3), 0, B(2, 4), 0, B(2, 5), 0, B(2, 6);
        B(2, 1), B(1, 1), B(2, 2), B(1, 2), B(2, 3), B(1, 3), B(2, 4), B(1, 4), B(2, 5), B(1, 5), B(2, 6), B(1, 6)];

    % Compute the strain
    strain = B_MAT * u_e;

    % Retrieve the current state variables (sdv)
    sdv = sdv_global(sdv_counter, :);

    % Elastic-plastic material behavior
    [stress, C_t, sdv, eps33] = elastic_plastic_von_mises_model_plane_stress(strain, sdv, props, hard_funct);

    % Update the global state variables
    sdv_global(sdv_counter, :) = sdv;
    sdv_counter = sdv_counter + 1;

    % Get the element coordinates and reshape to 12x1
    element_coord = reshape(node_coordinates(element, :)', [], 1);

    % Calculate s and t values
    x = N * element_coord;
    s_value = x(1);
    t_value = x(2);

    % Store stress, s, and t values in the matrix
    matrix(j, 1:3) = stress;
    matrix(j, 4) = s_value;
    matrix(j, 5) = t_value;

    % Compute the stiffness matrix
    Kt = B_MAT' * C_t * B_MAT * det(J) * w; % only b mat is in unit domain, so we multiply it with jacobian to transform to elemnt domain
    K = K + Kt;

    % Compute the internal forces
    Fint = B_MAT' * stress * det(J) * w;
    Fe_int = Fe_int + Fint;
end

% Calculate the external force Fe_ext
if any(cell2mat(elem_sets(2, 2)) == e)
    for j = 1:2
        elements = node_sets{2, 2};
        c = elem_list(e, 2:4);
        is_present = ismember(c, elements);

        if is_present(1) && is_present(3)
            t = se(j);
            s = 0;
        elseif is_present(1) && is_present(2)
            t = 0;
            s = se(j);
        else
            s = 1 - se(j);
            t = se(j);
        end

        trac = [tau_n; 0];
        w = we(j);

        N1 = 2*s^2 + 2*t^2 + 4*s*t - 3*t - 3*s + 1;
        N2 = 2*s^2 - s;
        N3 = 2*t^2 - t;
        N4 = 4*s - 4*s^2 - 4*s*t;
        N5 = 4*s*t;
        N6 = 4*t - 4*s*t - 4*t^2;

        N = [N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6, 0;
            0, N1, 0, N2, 0, N3, 0, N4, 0, N5, 0, N6];

        dN1_ds = 4*s + 4*t - 3;
        dN2_ds = 4*s - 1;
        dN3_ds = 0;
        dN4_ds = 4 - 8*s - 4*t;
        dN5_ds = 4*t;
        dN6_ds = -4*t;

        dN1_dt = 4*s + 4*t - 3;
        dN2_dt = 0;
        dN3_dt = 4*t - 1;
        dN4_dt = -4*s;
        dN5_dt = 4*s;
        dN6_dt = 4 - 4*s - 8*t;

        % Compute the dN_dE matrix
        dN_dE = [dN1_ds, dN2_ds, dN3_ds, dN4_ds, dN5_ds, dN6_ds;
            dN1_dt, dN2_dt, dN3_dt, dN4_dt, dN5_dt, dN6_dt];

        % Compute the Jacobian matrix J
        J = dN_dE * element_coordinates;
        J = J';

        if is_present(1) && is_present(3)
            detJ = norm(J(:, 2));
        elseif is_present(1) && is_present(2)
            detJ = norm(J(:, 1));
        else
            detJ = norm(J * [-1; 1]);
        end

        % Compute the external force Fe_ext using shape functions and traction
        Fe_ext = Fe_ext + N' * trac * detJ * w;
    end
end
end