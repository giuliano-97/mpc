% BRRIEF:
%   Template for explicit invariant set computation. You MUST NOT change
%   the output.
% OUTPUT:
%   A_x, b_x: Describes polytopic X_LQR = {x| A_x * x <= b_x}
function [A_x, b_x] = compute_X_LQR
    % get basic controller parameters
    param = compute_controller_base_parameters;
     % Here you need to implement the X_LQR computation and assign the result.
    
    system = LTISystem('A', param.A+param.B*param.F);
    Xp = Polyhedron('A',[eye(3); -eye(3); param.F; -param.F], 'b', [param.Xcons(:,2);-param.Xcons(:,1); param.Ucons(:,2);-param.Ucons(:,1)]);
    system.x.with('setConstraint');
    system.x.setConstraint = Xp;
    
    InvSetLQR = system.invariantSet()
%     figure()
%     InvSetLQR.plot()
    
    A_x = InvSetLQR.A;
    b_x = InvSetLQR.b;
end

