M_sym = [M11, M12, M13;
         M12, M22, M23;
         M13, M23, M33];

% Compute Christoffel symbols:
% C_ijk = 1/2 * (dM_ij/dtheta_k + dM_ik/dtheta_j - dM_jk/dtheta_i)
% We'll do this in a triple nested loop or systematically.

C_sym = sym('C',[3,3,3]);

theta_vec = [theta1; theta2; theta3];
for i=1:3
    for j=1:3
        for k=1:3
            C_sym(i,j,k) = ( diff(M_sym(i,j),theta_vec(k)) ...
                            + diff(M_sym(i,k),theta_vec(j)) ...
                            - diff(M_sym(j,k),theta_vec(i)) )/2;
        end
    end
end

% Now C(θ,θ̇) is computed by summing over k: C_ij * dtheta_k
% C(θ, θ̇) matrix:
C_matrix = sym('Cmat',[3 3]);
for i=1:3
    for j=1:3
        C_ij = 0;
        for k=1:3
            C_ij = C_ij + C_sym(i,j,k)*[dtheta1; dtheta2; dtheta3](k);
        end
        C_matrix(i,j) = expand(C_ij); % expand to get full expression
    end
end

disp('C(θ,θ̇) = ');
disp(C_matrix);
