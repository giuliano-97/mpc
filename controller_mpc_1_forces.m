% BRIEF:
%   Controller function template. This function can be freely modified but
%   input and output dimension MUST NOT be changed.
% INPUT:
%   T: Measured system temperatures, dimension (3,1)
% OUTPUT:
%   p: Cooling power, dimension (2,1)
function p = controller_mpc_1_forces(T)
% controller variables
persistent param forces_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, forces_optimizer] = init();
end

%% evaluate control action by solving MPC problem, e.g.
[u_mpc,errorcode] = forces_optimizer(T-param.T_sp);
if (errorcode ~= 1)
      warning('MPC infeasible');
end
p = u_mpc + param.p_sp;

end

function [param, forces_optimizer] = init()
% initializes the controller on first call and returns parameters and
% Yalmip optimizer object

param = compute_controller_base_parameters; % get basic controller parameters

%% implement your MPC using Yalmip2Forces interface here, e.g.
N = 30;
nx = size(param.A,1);
nu = size(param.B,2);

U = sdpvar(repmat(nu,1,N-1),repmat(1,1,N-1),'full');
X = sdpvar(repmat(nx,1,N),repmat(1,1,N),'full');

objective = 0;
constraints = [];
for k = 1:N-1
  % system dynamics constraints
  constraints = [constraints, X{k+1} == param.A*X{k} + param.B*U{k}];
  % state constraints
  constraints = [constraints, param.Xcons(:,1) <= X{k+1} <= param.Xcons(:,2)];
  % input constraints
  constraints = [constraints, param.Ucons(:,1) <= U{k} <= param.Ucons(:,2)];
  objective = objective + X{k}'*param.Q*X{k} + U{k}'*param.R*U{k};
end

%terminal cost
objective = objective + X{end}'*param.P*X{end};

% initial state 
x0 = sdpvar(3,1);
constraints = [constraints, X{1} == x0];

codeoptions = getOptions('solvername');
codeoptions.optlevel = 3;
fprintf('JMPC_dummy = %f',value(objective));
forces_optimizer = optimizerFORCES(constraints,objective,codeoptions,x0,U{1},{'xinit'}, {'u0'});
end

