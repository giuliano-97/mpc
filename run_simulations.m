% Init
clear all
close all
addpath(genpath(cd));
load('system/parameters_scenarios.mat');
figIdx = 1;

% clear persisten variables
clear controller_lqr; 
clear controller_mpc_1;
clear controller_mpc_2;
clear controller_mpc_3;
clear controller_mpc_4;
clear controller_mpc_5;
param = compute_controller_base_parameters;
sim_LQR = [true, true];
sim_MPC1 = [true, true];
sim_MPC2 = [true, true];
sim_MPC3 = [true, true, true];
sim_MPC4 = [true, true, true];
sim_MPC5 = [true, true, true;  % scen1
            true, true, true]; % scen2
sim_MPC_FORCES = true;

%% point 5
% Execute simulations starting from T_01 with scen1
% Using the LQR
if sim_LQR(1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_1, @controller_lqr, scen1);
    sgtitle('LQR - T_01 - scen1');
%     savefig("LQR_plot_T_01");
    if  norm(T(:,30) - param.T_sp) > 0.2 * norm(param.x0_1)
        warning('Requirement 5.2 not satisfied by LQR when starting from T_01!\n');
    end
end


%% point 6
J_inf = (T(:,1)-param.T_sp)'*param.P*(T(:,1)-param.T_sp);


%% point 7
% Using the LQR
if sim_LQR(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_lqr, scen1);
    sgtitle('LQR - T_02 - scen1');
end


%% point 8
% [A_x, b_x] = compute_X_LQR();


%% point 9
% Using the MPC controller MPC1
if sim_MPC1(1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, p] = simulate_truck(param.T0_1, @controller_mpc_1, scen1);
    sgtitle('MPC1 - T_01 - scen1');
    % Compute JMPC2(x0_1) - finite horizon cost
    JMPC1 = 0;
    for i= 1:30
        JMPC1 = JMPC1 + (T(:,i)-param.T_sp)'*param.Q*(T(:,i)-param.T_sp) + ...
            (p(:,i)-param.p_sp)'*param.R*(p(:,i)-param.p_sp);
    end
    % Add LQR terminal cost
    JMPC1 = JMPC1 + (T(:,i+1)-param.T_sp)'*param.P*(T(:,i+1)-param.T_sp);
end
% Using the MPC controller MPC1
if sim_MPC1(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_mpc_1, scen1);
    sgtitle('MPC1 - T_02 - scen1');
end


%% point 12
% Using the MPC controller MPC2
if sim_MPC2(1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, p] = simulate_truck(param.T0_1, @controller_mpc_2, scen1);
    sgtitle('MPC2 - T_01 - scen1');
    % Compute JMPC2(x0_1)
    JMPC2 = 0;
    for i= 1:30
        JMPC2 = JMPC2 + (T(:,i)-param.T_sp)'*param.Q*(T(:,i)-param.T_sp) + ...
            (p(:,i)-param.p_sp)'*param.R*(p(:,i)-param.p_sp);
    end
end


%% point 13
% Using the MPC controller MPC2
if sim_MPC2(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_mpc_2, scen1);
    sgtitle('MPC2 - T_02 - scen1');
end


%% point 15
% Using the MPC controller MPC3
if sim_MPC3(1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_1, @controller_mpc_3, scen1);
    sgtitle('MPC3 - T_01 - scen1');
end

% Using the MPC controller MPC3
if sim_MPC3(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_mpc_3, scen1);
    sgtitle('MPC3 - T_02 - scen1');
end


%% point 17
% Using the MPC controller MPC3
if sim_MPC3(3)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_3, @controller_mpc_3, scen1);
    sgtitle('MPC3 - T_03 - scen1');
end


%% point 18
% Using the MPC controller MPC4
if sim_MPC4(3)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_3, @controller_mpc_4, scen1);
    sgtitle('MPC4 - T_03 - scen1');
end


%% point 19
% Using the MPC controller MPC3
if sim_MPC3(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_mpc_3, scen1);
    sgtitle('MPC3 - T_02 - scen1');
end


% Using the MPC controller MPC4
if sim_MPC4(2)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_2, @controller_mpc_4, scen1);
    sgtitle('MPC4 - T_02 - scen1');
end


%% point 22
% Using the MPC controller MPC5
clear controller_mpc_5
if sim_MPC5(2,1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_1, @controller_mpc_5, scen2);
    sgtitle('MPC5 - T_01 - scen2');
end

% Using the MPC controller MPC3
if sim_MPC3(1)
    figure(figIdx); figIdx = figIdx + 1; 
    [T, ~] = simulate_truck(param.T0_1, @controller_mpc_3, scen2);
    sgtitle('MPC3 - T_01 - scen2');
end


%% point 23
% execute simulation starting from T0_2 using FORCES Pro MPC controller with scenario 1 and get the average
%  running time
if sim_MPC_FORCES
    %initialization
    figure(figIdx); figIdx = figIdx + 1; 
    [~,~,t_sim] = simulate_truck(param.T0_2, @controller_mpc_1, scen1);
    sgtitle('MPC1 - T_02 - scen1');
    
    figure(figIdx); figIdx = figIdx + 1; 
    [T_force,p_force,t_sim_forces] = simulate_truck(param.T0_2, @controller_mpc_1_forces, scen1);
    sgtitle('MPC1 FORCES - T_02 - scen1');
%     t = [];
%     t_forces = [];
%     for i = 1:10
%     [~,~,t_sim] = simulate_truck(param.T0_2, @controller_mpc_1, scen1);
%     [T_force,p_force,t_sim_forces] = simulate_truck(param.T0_2, @controller_mpc_1_forces, scen1);
%     % [~,~,t_sim_forces] = simulate_truck(T0_2, @controller_mpc_1_forces, scen1);
%     t = [t; t_sim];
%     t_forces = [t_forces; t_sim_forces];
%     end
%     fprintf('average running times for controller_mpc_1_forces is %s\n',num2str(mean(t_forces)))
%     fprintf('average running times for controller_mpc_1 is %s\n',num2str(mean(t)))
end



