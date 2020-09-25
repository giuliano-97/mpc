function param = compute_controller_base_parameters
    % load truck parameters
    load('system/parameters_truck');
    
    % CT parameters    
    Bcd = diag([1/truck.m1, 1/truck.m2, 1/truck.m3]);

    A_ = [-(truck.a1o + truck.a12),                            truck.a12,                         0;
                         truck.a12, -(truck.a12 + truck.a23 + truck.a2o),                 truck.a23;
                                 0,                            truck.a23, -(truck.a23 + truck.a3o)];

    Bc_ = [1, 0; 
           0, 1;
           0, 0];
    
    Ac = Bcd * A_;
    Bc = Bcd * Bc_;

    dc = truck.w + [truck.a1o * truck.To; 
                    truck.a2o * truck.To; 
                    truck.a3o * truck.To];

    
    % (2) Discretization (Euler)
    Ts = 60;    
    A = (eye(3) + Ts*Ac);
    B = Ts*Bc;
    Bd = Ts*Bcd;
    
    % (3) Set point computation
    d = Bd*dc;
    T1ss = truck.b_ref(1);
    T2ss = truck.b_ref(2);
    
    S = eye(3) - A;
    K = [B, -S(:,3)];
    H = [S(:,1:2), -d];
    b = [T1ss;T2ss;1];

    ss = K \ (H*b);

    T_sp = [T1ss;T2ss;ss(3)];
    p_sp = [ss(1);ss(2)];
    
    % (4) system constraints
    Pcons = truck.InputConstraints;
    Tcons = truck.StateConstraints;
    
    % (4) constraints for delta formulation
    Ucons = Pcons - p_sp;
    Xcons = Tcons - T_sp;
    
    % (5) LQR cost function
    Q = 3850*eye(3);
    R = 0.1*eye(2);
    
    % put everything together
    param.A = A;
    param.B = B;
    param.C = eye(3);
    param.Bd = Bd;
    param.dc = dc;
    param.Q = Q;
    param.R = R;
    param.Ts = Ts;
    param.T_sp = T_sp;
    param.p_sp = p_sp;
    param.Ucons = Ucons;
    param.Xcons = Xcons;
    param.Tcons = Tcons;
    param.Pcons = Pcons;
    
    % First set of initial conditions    
    param.x0_1 = [3;1;0];
    param.T0_1 = T_sp + param.x0_1;
    
    % Second set of initial conditions
    param.x0_2 = [-1; -0.3; -4.5];
    param.T0_2 = T_sp + [-1; -0.3; -4.5];
    
    % Third set of initial conditions
    param.x0_3 = [12, 12, 12] - T_sp;
    param.T0_3 = [12, 12, 12];
    
    % (5) Precompute LQR controller parameters
    [F,P] = dlqr(param.A,param.B,param.Q,param.R);
    param.F = -F;
    param.P = P;
    
    % (18) Soft constrained MPC parameters
    param.ni = 500000;
    param.S = 100;

    param.Ax = [-1,0,0;  1,0,0;  0,-1,0;  0,1,0;  0,0,-1;  0,0,1];
    param.Au = [-1,0;1,0;0,-1;0,1];
    param.bx = [-Xcons(1,1); Xcons(1,2); -Xcons(2,1); Xcons(2,2); -Xcons(3,1); Xcons(3,2)];
    param.bu = [-Ucons(1,1); Ucons(1,2); -Ucons(2,1); Ucons(2,2)];
    
    % (22) Offset-free MPC
    param.A_aug = [A, Bd; zeros(3,3), eye(3)];
    param.B_aug = [B; zeros(3,2)];
    param.C = eye(3);
    param.Cd = zeros(3);
    param.C_aug = [param.C, param.Cd];

    poles = [0.4,0.3,0.1,0.2,0.1,0.2];
    param.L = (place(param.A_aug', param.C_aug', poles))';

    param.r = [T1ss; T2ss];
end

