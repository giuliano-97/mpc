% BRIEF:
%   Controller function template. This function can be freely modified but
%   input and output dimension MUST NOT be changed.
% INPUT:
%   T: Measured system temperatures, dimension (3,1)
% OUTPUT:
%   p: Cooling power, dimension (2,1)
function p = controller_mpc_4(T)
% controller variables
persistent param yalmip_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init();
end

%% Evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T - param.T_sp);
if (errorcode ~= 0)
      warning('MPC infeasible');
end
p = u_mpc + param.p_sp;
end

function [param, yalmip_optimizer] = init()
% initializes the controller on first call and returns parameters and
% Yalmip optimizer object

param = compute_controller_base_parameters; % get basic controller parameters

% Get LQR polytopic invariant set
[A_x, b_x] = compute_X_LQR;

% Dims
N = 30;
nx = size(param.A,1);
nu = size(param.B,2);
nt = size(b_x,1);

% State and control input
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');

% Slack variables
E = sdpvar(repmat(nx*2,1,N),ones(1,N),'full'); % Box constraints
Et = sdpvar(nt,1);  % Terminal set constraints

objective = 0;
% This this we enforce the state constraints (w/ slack) also on the
% initial conditions
constraints = E{1} >= zeros(nx*2,1);
% constraints = [constraints, param.Xcons(:,1) - E{1}(1:nx) <= X{1}];
% constraints = [constraints, X{1} <= param.Xcons(:,2) + E{1}((nx+1):end)];
for k = 1:N-1
    % Add term to the objective
    objective = objective +  X{k}'*param.Q*X{k} + U{k}'*param.R*U{k} + ...
        E{k}'*param.S*E{k} + ... % Quadratic slack penalty
        param.ni * norm(E{k},1); % Linear slack penalty
    % Add dynamics equality constraints
    constraints = [constraints, X{k+1} == param.A * X{k} + param.B * U{k}];
    % Add state box constraint with slack
  constraints = [constraints, param.Ax * X{k+1} <= param.bx + E{k+1}];
    % Add slackness positiveness constraint
    constraints = [constraints, E{k+1} >= zeros(nx*2,1)];
    % Add input box constraint
    constraints = [constraints, param.Ucons(:,1) <= U{k} <= param.Ucons(:,2)];
end

% Terminal cost
objective = objective + X{N}'*param.P*X{N} + ...
    E{N}'*param.S*E{N} + ...    % Quadratic slack penalty
    param.ni * norm(E{N},1);    % Linear slack penalty;

% Terminal state constraints relaxation and slack variable penalty
constraints = [constraints, A_x * X{N} <= b_x + Et, Et >= zeros(nt,1)];
objective = objective + param.ni*norm(Et,1) + Et'*param.S*Et;


ops = sdpsettings('verbose',0,'solver','quadprog');
fprintf('JMPC_dummy = %f',value(objective));
yalmip_optimizer = optimizer(constraints,objective,ops,X{1},U{1});
end

