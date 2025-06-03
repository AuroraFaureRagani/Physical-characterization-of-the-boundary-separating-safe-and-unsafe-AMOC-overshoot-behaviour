%% --- OVERSHOOT TRAJECTORIES WITH INCREASING LINEAR FORCING ---
% Creates panels a-b in Figure 8

% Clear workspace
clearvars
close all

% Define simulation
sim = [3, 4, 5, 6];

% Load simulation parameters from Excel file
table_sim = xlsread('data/sim_PL.xlsx'); 
sim_PL = table_sim(:,1);                  % Simulation indices (used for indexing rows)
t1 = table_sim(:,2);                      % Time when ramp starts
t2 = table_sim(:,3);                      % Time when ramp ends
gamma_max = table_sim(:,4);               % Maximum forcing amplitude
label = {'a', 'b', 'c', 'd'};             % Labels for plotting

% Load bifurcation diagram data
bif_diagr = dlmread("data/p.AUR1");

% Physical constants and scaling factors
factatl = -0.913 * 8.308e-02;
threshold = [0.1855, 5.1793];             % [gamma_threshold, Psi_threshold]
T0 = 15;                                  % Reference temperature [°C]
S0 = 35;                                  % Reference salinity [psu]

% Define gamma (freshwater forcing) function with a piecewise linear ramp
dgamma_PL = @(t1, t2, ymax, t) -(...
    heaviside(t1 - t) .* ((ymax / t1) * t) + ...
    heaviside(t - t1) .* heaviside(t2 - t) .* (((threshold(1)/2 - ymax) / (t2 - t1)) .* t + ymax - t1 * ((threshold(1)/2 - ymax) / (t2 - t1))) + ...
    heaviside(t - t2) .* threshold(1)/2 ...
    ) / factatl;

% Set color scheme for plotting
color = lines(numel(sim) + 1);

%% Loop through simulations and plot results

k = 0;  % Index for color and label tracking
for i = sim
    k = k + 1;
    
    % Load simulation results and compute freshwater forcing
    anal1 = dlmread("data/p.TPL" + i);
    time_i = anal1(:,3);   % Time in years
    gamma_i = -dgamma_PL(t1(sim_PL(i)), t2(sim_PL(i)), gamma_max(sim_PL(i)), time_i) * factatl;
    

    % Plot freshwater forcing after year 355
    figure(1)
    kk = find(time_i == 355);  % Starting point for plotting (transition time)
    plot(time_i(kk:end), gamma_i(kk:end), 'LineWidth', 2, 'Color', color(k, :), 'DisplayName', "sim " + label{k});
    hold on

    % Plot state vs forcing (bifurcation diagram)
    figure(2)
    plot(gamma_i(kk:end), anal1(kk:end,5), 'LineWidth', 2, 'Color', color(k, :), 'DisplayName', "sim " + label{k});
    hold on

    % Plot earlier trajectory in black (before transition)
    plot(gamma_i(1:kk), anal1(1:kk,5), 'LineWidth', 1.2, 'Color', 'k', 'HandleVisibility', 'off');
    hold on
end

%% Finalize freshwater forcing plot
figure(1)
plot(time_i(1:kk), gamma_i(1:kk), 'LineWidth', 1.2, 'Color', 'k', 'HandleVisibility', 'off');  % Pre-transition trajectory
yline(threshold(1), 'LineWidth', 1.2, 'HandleVisibility', 'off');  % Critical gamma threshold
xlabel('time (years)', 'FontSize', 16);
ylabel('\gamma_A  (Sv)', 'FontSize', 16);
xlim([0, 600]);
grid on
title('Freshwater forcing', 'FontSize', 16);

%% Finalize bifurcation diagram plot
figure(2)
plot_bif_diagr_function();  % Plot background bifurcation diagram (custom function)
grid on
xlim([0, 0.3]);
xlabel('\gamma_A  (Sv)', 'FontSize', 16);
ylabel('\Psi_A  (Sv)', 'FontSize', 16);
title('AMOC strength', 'FontSize', 16);
legend('FontSize', 16);  % Add legend with simulation labels
