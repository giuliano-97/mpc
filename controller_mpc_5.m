% BRIEF:
%   Controller function template. This function can be freely modified but
%   input and output dimension MUST NOT be changed.
% INPUT:
%   T: Measured system temperatures, dimension (3,1)
% OUTPUT:
%   p: Cooling power, dimension (2,1)
function p = controller_mpc_5(T)
% controller variables
persistent param yalmip_optimizer curr_T_hat curr_d_hat

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init();
    curr_T_hat = T;
    curr_d_hat = param.dc;
end

% Get steady state params with current state and disturbance estimates
[Ts,ps] = getSS(param, curr_d_hat);

% Evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(...
    {T, curr_d_hat, Ts, ps});
if (errorcode ~= 0)
      warning('MPC infeasible');
end
p = u_mpc;

% Re-estimate
[curr_T_hat, curr_d_hat] = estimate(param, curr_T_hat, curr_d_hat, T, p);

end

%% Compute ss state and control input
function [Ts, ps] = getSS(param, curr_d_hat)

    A = [param.A - eye(3), param.B;
         param.H*param.C,       zeros(2,2)];
     
    b = [-param.Bd * curr_d_hat;
         param.r-param.H*param.Cd*curr_d_hat];
    
    x = A \ b;
    
    Ts = x(1:3);
    ps = x(4:5);
end

%% Step state estimator
function [curr_T_hat, curr_d_hat] = estimate(param, curr_T_hat, curr_d_hat, T_out, u_curr)
    Td_hat = [curr_T_hat; curr_d_hat];
    Td_hat = param.A_aug * Td_hat + param.B_aug * u_curr + ...
        param.L * (T_out - param.C_aug * Td_hat);
   
    curr_T_hat = Td_hat(1:3);
    curr_d_hat = Td_hat(4:6);
    % disp(curr_d_hat);
end


function [param, yalmip_optimizer] = init()
% initializes the controller on first call and returns parameters and
% Yalmip optimizer object

param = compute_controller_base_parameters; % get basic controller parameters



param.H = [1,0,0;
           0,1,0];
    
% Get LQR polytopic invariant set
[A_x, b_x] = compute_X_LQR;

% Dims
N = 30;
nx = size(param.A,1);
nu = size(param.B,2);

% State and control input
P = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
T = sdpvar(repmat(nx,1,N),ones(1,N),'full');
D = sdpvar(repmat(nx,1,N),ones(1,N),'full');

Ts = sdpvar(nx,1,'full');
Ps = sdpvar(nu,1,'full');

% Slack variables
objective = 0;
constraints = [];
for k = 1:N-1
    % Add term to the objective
    objective = objective +  (T{:,k} - Ts)'*param.Q*(T{:,k} - Ts) ...
        + (P{:,k} - Ps)'*param.R*(P{:,k} - Ps);
    % Add dynamics equality constraints
    constraints = [constraints, ...
        T{:,k+1} == param.A * T{:,k} + param.B * P{:,k} + param.Bd * D{:,k}];
    constraints = [constraints, ...
        D{k+1} == D{k}];
    % Add state box constraint
    constraints = [constraints, param.Tcons(:,1) <= T{k+1} <= param.Tcons(:,2)];
    % Add input box constraint
    constraints = [constraints, param.Pcons(:,1) <= P{k} <= param.Pcons(:,2)];
end

% Terminal cost
objective = objective + (T{N} - Ts)'*param.P*(T{N} - Ts);

% Relaxed terminal state constraints
constraints = [constraints, A_x * (T{N} - Ts) <= b_x];

ops = sdpsettings('verbose',0,'solver','quadprog');
fprintf('JMPC_dummy = %f',value(objective));
yalmip_optimizer = optimizer(constraints,objective,ops, ...
    {T{1}, D{1}, Ts, Ps}, P{1});
end
