%% --- OVERSHOOT TRAJECTORIES WITH INCREASING LINEAR FORCING ---
% Creates panels in Figure 1

% Clear workspace
clearvars
close all

% Define simulations and corresponding rates of forcing change
sim = [2, 4];
rate = -[1/2, 1/8];

% Build file names for loading simulation data
fname = "data/p.BI2_T" + string(sim);
% --- DEFINE CONSTANTS ---
udim  = 0.1;                   % [m/s]   Characteristic velocity scale
r0dim = 6.4e6;                 % [m]     Earth's radius
factatl     = -0.913 * 8.308e-02;               % Scaling factor relative to Atlantic Ocean
factrate    = 2.64e-03;                         % [10^{-4} Sv/yr] Conversion factor for rate
facttime    = r0dim / udim / (360*24*60*60);    % [yr] Time scale (from velocity and Earth's radius)



% Define critical threshold values (for gamma and streamfunction)
threshold = [0.1855, 5.1793];
gamma_0 = 0;
gamma_lim = threshold(1) / 2; 

% --- PLOT FRESHWATER FORCING vs TIME ---
figure(1)
yline(threshold(1), 'LineWidth', 2, 'HandleVisibility','off'); 

for i = 1:numel(sim)

    % Load file
    xx = dlmread(fname(i));
    
    % Define time and freshwater forging
    time = xx(:,3);
    gamma = time * rate(i) * factrate * factatl;
    
    % Plot freshwater forcing vs time
    figure(1)
    plot(time, gamma, 'LineWidth', 2);    
    hold on

    % --- PLOT BIFURCATION PROGRESSION (STATE vs FORCING) ---
    figure(2)
    plot(gamma, xx(:,5), 'LineWidth', 2); 
    hold on    
end

% Finalize figure 1 (Freshwater forcing vs time)
figure(1)
ylabel('\gamma_A (Sv)', 'FontSize', 16);
xlabel('time (years)', 'FontSize', 16);
title('Freshwater forcing', 'FontSize', 16);
grid on

% Finalize figure 2 (Bifurcation diagram)
figure(2)
xlim([0, 0.3]);
bif_diagr = dlmread("data/p.bifOG");
hold on
plot_bif_diagr_function()
xlabel('\gamma_A (Sv)', 'FontSize', 16);
ylabel('\Psi_A (Sv)', 'FontSize', 16);
title('Bifurcation diagram', 'FontSize', 16);
grid on
