clear; clc;
syms theta1 theta2 theta3 real
syms dtheta1 dtheta2 dtheta3 real
syms m1 m2 m3 L1 L2 lc1 lc2 lc3 I1 I2 I3 g real

% Define the absolute angles (as given)
alpha1 = theta1;
alpha2 = theta1 + theta2;
alpha3 = theta1 + theta2 + theta3;

% Positions of COMs:
x_c1 = lc1*cos(alpha1); y_c1 = lc1*sin(alpha1);
x_c2 = L1*cos(alpha1) + lc2*cos(alpha2); y_c2 = L1*sin(alpha1) + lc2*sin(alpha2);
x_c3 = L1*cos(alpha1) + L2*cos(alpha2) + lc3*cos(alpha3);
y_c3 = L1*sin(alpha1) + L2*sin(alpha2) + lc3*sin(alpha3);

% Velocities:
dx_c1 = diff(x_c1,theta1)*dtheta1 + diff(x_c1,theta2)*dtheta2 + diff(x_c1,theta3)*dtheta3;
dy_c1 = diff(y_c1,theta1)*dtheta1 + diff(y_c1,theta2)*dtheta2 + diff(y_c1,theta3)*dtheta3;
dx_c2 = diff(x_c2,theta1)*dtheta1 + diff(x_c2,theta2)*dtheta2 + diff(x_c2,theta3)*dtheta3;
dy_c2 = diff(y_c2,theta1)*dtheta1 + diff(y_c2,theta2)*dtheta2 + diff(y_c2,theta3)*dtheta3;
dx_c3 = diff(x_c3,theta1)*dtheta1 + diff(x_c3,theta2)*dtheta2 + diff(x_c3,theta3)*dtheta3;
dy_c3 = diff(y_c3,theta1)*dtheta1 + diff(y_c3,theta2)*dtheta2 + diff(y_c3,theta3)*dtheta3;

% Angular velocities:
omega1 = dtheta1;
omega2 = dtheta1 + dtheta2;
omega3 = dtheta1 + dtheta2 + dtheta3;

% Kinetic Energy:
T = (1/2)*m1*(dx_c1^2+dy_c1^2) + (1/2)*I1*omega1^2 ...
  + (1/2)*m2*(dx_c2^2+dy_c2^2) + (1/2)*I2*omega2^2 ...
  + (1/2)*m3*(dx_c3^2+dy_c3^2) + (1/2)*I3*omega3^2;

% Extract Mass Matrix M(θ):
T_expanded = expand(T);
M11 = diff(diff(T_expanded,dtheta1),dtheta1);
M12 = diff(diff(T_expanded,dtheta1),dtheta2);
M13 = diff(diff(T_expanded,dtheta1),dtheta3);
M22 = diff(diff(T_expanded,dtheta2),dtheta2);
M23 = diff(diff(T_expanded,dtheta2),dtheta3);
M33 = diff(diff(T_expanded,dtheta3),dtheta3);

% Due to symmetry:
M_full = [M11, M12, M13;
          M12, M22, M23;
          M13, M23, M33];

disp('Mass Matrix (M_full):');
disp(M_full);

% Potential Energy V:
V = m1*g*y_c1 + m2*g*y_c2 + m3*g*y_c3;

% Gravity vector N(θ) = [∂V/∂θ1; ∂V/∂θ2; ∂V/∂θ3]
N1 = diff(V,theta1);
N2 = diff(V,theta2);
N3 = diff(V,theta3);
N_full = [N1; N2; N3];

disp('Gravity vector (N_full):');
disp(N_full);

% Coriolis/centrifugal terms from Christoffel symbols:
% C_ijk = 1/2(dM_ij/dθ_k + dM_ik/dθ_j - dM_jk/dθ_i)
theta_vec = [theta1; theta2; theta3];
C_sym = sym('C',[3,3,3]);
for i=1:3
    for j=1:3
        for k=1:3
            C_sym(i,j,k) = (diff(M_full(i,j),theta_vec(k)) + diff(M_full(i,k),theta_vec(j)) - diff(M_full(j,k),theta_vec(i)))/2;
        end
    end
end

% C(θ,θ̇): C_ij = sum_k C_ijk dθ_k
C_matrix = sym('Cmat',[3 3]);
for i=1:3
    for j=1:3
        tmp = 0;
        for k=1:3
            tmp = tmp + C_sym(i,j,k)*[dtheta1; dtheta2; dtheta3]*(k);
        end
        C_matrix(i,j) = expand(tmp);
    end
end

disp('Coriolis matrix (C_matrix):');
disp(C_matrix);
