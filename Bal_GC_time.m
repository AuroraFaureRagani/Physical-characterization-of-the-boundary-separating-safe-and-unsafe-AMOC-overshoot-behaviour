%% --- SALT BALANCES STUDY ---
% Creates 
% panels a-b of Figure 2
% panles a-b of Figure 3
% panels a-d of Figure 4
% panels a-b of Figure 5


% Clear workspace
clear all
close all

% Load simulation parameters
table_sim = xlsread('data/sim_PL.xlsx');

% Extract relevant columns from the table
sim_PL = table_sim(:,1);       % Index of simulations
t1 = table_sim(:,2);           % First time of forcing
t2 = table_sim(:,3);           % Second time of forcing
gamma_max = table_sim(:,4);    % Maximum value of freshwater forcing

%% - READ BIFURCATION DIAGRAM/TRANSIENTS -

% Load bifurcation and transient simulation data

% Simulation cases:
% 1 - 2 have same peak and growth, but different relaxation rates
%     1: recovery (safe)
%     2: collapse (unsafe, shown as dashed line)

i = 1; j = 2;

% Load transient simulations
anal1 = dlmread("data/p.TPL" + i);
anal2 = dlmread("data/p.TPL" + j);
bif_diagr = dlmread("data/p.bifOG");

% Load bifurcation diagram
threshold = [0.1855 , 5.1793];
factatl = -0.913*8.308e-02;
ft      = 1000*0.035*0.1*4000*6.3e06/1.0e09;
T0      = 15;                  %[deg C]  Reference temperature
S0      = 35;                  %[psu]    Reference salinity
alphaS  = 0.76;                %[g/kg]   Haline contraction coefficient
alphaT  = 0.1;                 %[kg/m^3] Thermal expansion coefficient

% Piecewise linear forcing function
dgamma_PL = @(t1, t2, ymax, t) -(heaviside(t1-t).*((ymax/t1)*t) + ...
    heaviside(t-t1).*heaviside(t2-t).*(((threshold(1)/2-ymax)/(t2-t1)).*t + ymax-t1.*((threshold(1)/2-ymax)/(t2-t1))) + ...
    heaviside(t-t2).*threshold(1)/2)/factatl;

% Coordinates of bifurcation tipping point for reference in plots
x_pos1 = -0.24550E+01*factatl;
y_pos1 = 0.53833E+01; 

%% Compute freshwater forcing time series
time_i = anal1(:,3);
gamma_i = -dgamma_PL(t1(sim_PL(i)), t2(sim_PL(i)), gamma_max(sim_PL(i)), time_i)*factatl;

time_j = anal2(:,3);
gamma_j = -dgamma_PL(t1(sim_PL(j)), t2(sim_PL(j)), gamma_max(sim_PL(j)), time_j)*factatl;

%% Compute meridional differences (north-south)

% Ocean basin volumes
vols = 0.24931;
voln = 0.13207;

% Salt components for both simulations
salt1 = anal1(:, 24);
salt2 = anal2(:, 24);

% Density, salinity and temperature differences (north-south) for both simulations
dS1 = anal1(:, end-2)/voln - anal1(:, end-1)/vols;
dT1 = anal1(:, end-5)/voln - anal1(:, end-4)/vols;
drho1 = anal1(:, end-8)/voln - anal1(:, end-7)/vols;
dS2 = anal2(:, end-2)/voln - anal2(:, end-1)/vols;
dT2 = anal2(:, end-5)/voln - anal2(:, end-4)/vols;
drho2 = anal2(:, end-8)/voln - anal2(:, end-7)/vols;

% Advective, diffusive and net freshwater flux components
advn1 = anal1(:, 12);
advs1 = anal1(:, 13);
difn1 = anal1(:, 14);
difs1 = anal1(:, 15);
emp1 = anal1(:, 16);
check1 = anal1(:, 17);

advn2 = anal2(:, 12);
advs2 = anal2(:, 13);
difn2 = anal2(:, 14);
difs2 = anal2(:, 15);
emp2 = anal2(:, 16);
check2 = anal2(:, 17);

% North/South salinity for later use
Snorth1 = anal1(:, end-2);
Ssouth1 = anal1(:, end-1);
Snorth2 = anal2(:, end-2);
Ssouth2 = anal2(:, end-1);

%% Plot freshwater forcing
figure()
hold on
yline(x_pos1, 'HandleVisibility','off', 'LineWidth', 1.2)
plot(time_i, gamma_i, 'LineWidth', 2, 'color', [0 0.4470 0.7410])
plot(time_j, gamma_j, '--', 'LineWidth', 2, 'Color',  [0.4660 0.6740 0.1880]), hold on %'color', [0.4660 0.6740 0.1880]

xlabel('time (years)','FontSize', 16)
ylabel('\gamma_A  (Sv)', 'FontSize', 16)

xlim([0,1000])
ylim([0 0.3])
grid on
title('Freshwater forcing', 'FontSize',16)


%% Plot bifurcation diagram and transients
figure()
plot_bif_diagr_function() % Custom function to plot bifurcation background

% Add transient trajectories on top
plot(gamma_i, anal1(:,5), 'LineWidth', 2, 'color', [0 0.4470 0.7410], 'DisplayName', 'recovery')
plot(gamma_j, anal2(:,5), '--', 'LineWidth', 2, 'Color',  [0.4660 0.6740 0.1880], 'DisplayName', 'collapse') %'color', [0.4660 0.6740 0.1880]
grid on
xlim([0, 0.3])
xlabel('\gamma_A  (Sv)', 'FontSize', 16)
ylabel('\Psi_A  (Sv)', 'FontSize', 16)
title('AMOC strength', 'FontSize',16)

lgd = legend('FontSize',16);


%% Plot time evolution of meridional density difference and AMOC

figure()
colororder({'k','k'})

% Left y-axis: density difference
yyaxis left
plot(time_i, drho1, '-','LineWidth',2, 'Color',"#77AC30", 'DisplayName', '$\Delta \rho$'), hold on
plot(time_j, drho2, '--','LineWidth',2, 'Color',"#77AC30", 'HandleVisibility', 'off')
yline(0, 'HandleVisibility', 'off', 'LineWidth', 1.2)
ylabel('density difference (kg m^{-3})', 'FontSize',16)

yyaxis right

% Right y-axis: AMOC strength
plot(time_i, anal1(:,5),'-','LineWidth',1.2, 'Color',"k", 'DisplayName', ...
    'AMOC'), hold on
plot(time_j, anal2(:,5),'--','LineWidth',1.2, 'Color',"k", ...
    'HandleVisibility', 'off'), hold on
ylabel('AMOC strenght (Sv)', 'FontSize',16)

% Mark forcing intervals
xline(t1(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')

grid on
xlabel('time (years)', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 2;
title('Meridional density difference', 'FontSize',16)


%% Plot AMOC vs density difference
figure()
plot(drho1(1:400), anal1(1:400,5), 'LineWidth',2, 'color', [0 0.4470 0.7410], 'DisplayName', 'recovery'), hold on
plot(drho2, anal2(:,5), '--', 'LineWidth',2, 'color', [0.4660 0.6740 0.1880], 'DisplayName', 'collapse'); hold on;

grid on
xlabel('density difference (kg m^{-3})', 'FontSize', 16)
ylabel('AMOC strenght (Sv)', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 2;
title('AMOC versus density difference', 'FontSize', 16)

%% Plot meridional temperature and salinity contributions to density

% Temperature contribution
figure()
colororder({'k','k'})
plot(time_i, 10*dT1*alphaT, '-','LineWidth',2, 'Color',[0 0.4470 0.7410], 'DisplayName', '$\Delta T$'), hold on
plot(time_j, 10*dT2*alphaT, '--','LineWidth',2, 'Color',[0 0.4470 0.7410], 'HandleVisibility', 'off')

xline(t1(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')

grid on
ylabel('\DeltaT density contribution (kg/m^3) ' , 'FontSize',16)
xlabel('time (years)', 'FontSize',16)
title('Meridional temperature difference', 'FontSize',16)

% Salinity contribution
figure()
plot(time_i, 10*dS1*alphaS, '-','LineWidth',2, 'Color',"#A2142F", 'DisplayName', '$\Delta S$'), hold on
plot(time_j, 10*dS2*alphaS, '--','LineWidth',2, 'Color',"#A2142F", 'HandleVisibility', 'off')

xline(t1(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',1.2,'HandleVisibility', 'off')

grid on
ylabel('\DeltaS density contribution (kg/m^3)', 'FontSize',16)
xlabel('time (years)', 'FontSize',16)
title('Meridional salinity difference', 'FontSize',16)

%% salt balances - whole atlantic
figure()

% Plot time derivative of the total Atlantic salt content
plot(time_i, ft*gradient(salt1),'-','LineWidth',2, 'Color',"#77AC30", 'DisplayName', ...
     '$\frac{d}{dt} \int_{Atl} S \, dV $'), hold on
plot(time_j, ft*gradient(salt2),'--','LineWidth',2, 'Color',"#77AC30", ...
     'HandleVisibility', 'off'), hold on

% Plot E-P (evaporation minus precipitation) term
plot(time_i, ft*anal1(:,16),'-','LineWidth',2,'Color','#7E2F8E', 'DisplayName', ...
    '$\Phi^s$')
    %'$\Phi^s$') 
hold on
plot(time_j, ft*anal2(:,16),'--','LineWidth',2,'Color','#7E2F8E', 'HandleVisibility', 'off') 
hold on

% Plot net lateral salt fluxes into the Atlantic (sum over four boundaries)
plot(time_i, -ft*(anal1(:,12)-anal1(:,14)-(anal1(:,13)-anal1(:, 15))),'-','LineWidth',2,'Color', 'b', 'DisplayName', ...
    '$\Phi^{lat}$') 
hold on
plot(time_j, -ft*(anal2(:,12)-anal2(:,14)-(anal2(:,13)-anal2(:, 15))),'--','LineWidth',2,'Color', 'b', 'Handlevisibility', 'off') 
hold on

% Plot residual in the salt budget
plot(time_i, (ft*(anal1(:,17)+gradient(salt1))),'-','LineWidth',2,'Color',[0 0 0], 'DisplayName', '$\Phi^b$') 
hold on
plot(time_j, (ft*(anal2(:,17)+gradient(salt2))),'--','LineWidth',2,'Color',[0 0 0], 'Handlevisibility', 'off') 
hold on

% Add reference lines and formatting
yline(0, 'HandleVisibility', 'off')
xline(t1(i), 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',0.5,'HandleVisibility', 'off')
box('on');
xlim([0 max([time_i; time_j])])


% Label and format the plot
title('Salt balance vs time', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 4;
xlabel('time (years)', 'FontSize',16);
ylabel('salt transport (10^9 kg/s)', 'FontSize',16)


%% Detailed Breakdown of Lateral Salt Flux Terms (Whole Atlantic)
figure()

% Plot individual advective salt fluxes at northern and southern boundaries
plot(time_i, ft*advn1,'-','LineWidth',1,'Color','#4DBEEE', 'DisplayName','$\Phi^a(\theta_n)$'), hold on
plot(time_j, ft*advn2,'--','LineWidth',1,'Color','#4DBEEE', 'HandleVisibility', 'off') 

plot(time_i, ft*advs1,'-','LineWidth',1,'Color',[1 0 0], 'DisplayName','$\Phi^a(\theta_{s})$' )
plot(time_j, ft*advs2,'--','LineWidth',1,'Color',[1 0 0], 'HandleVisibility', 'off' )

% Plot individual diffusive salt fluxes at northern and southern boundaries
plot(time_i, ft*difn1,'-','LineWidth',1,'Color','#77AC30', 'DisplayName', '$\Phi^d(\theta_n)$')
plot(time_j, ft*difn2,'--','LineWidth',1,'Color','#77AC30', 'HandleVisibility', 'off') 

plot(time_i, ft*difs1,'-','LineWidth',1,'Color',[0 0 1], 'DisplayName', '$\Phi^d(\theta_s)$') 
plot(time_j, ft*difs2,'--','LineWidth',1,'Color',[0 0 1], 'HandleVisibility', 'off') 
hold on

% Plot total lateral flux for reference
plot(time_i, -ft*(anal1(:,12)-anal1(:,14)-(anal1(:,13)-anal1(:, 15))),'-','LineWidth',2,'Color', 'b', 'DisplayName', '$\Phi^{lat}$') 
plot(time_j, -ft*(anal2(:,12)-anal2(:,14)-(anal2(:,13)-anal2(:, 15))),'--','LineWidth',2,'Color', 'b', 'Handlevisibility', 'off') 

% Add reference lines and formatting
yline(0, 'HandleVisibility', 'off')
xline(t1(i), 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',0.5,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',0.5,'HandleVisibility', 'off')
box('on');
xlim([0 max([time_i; time_j])])

title('Lateral salt fluxes vs time', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 3;
xlabel('time (years)', 'FontSize',16)
ylabel('salt transport (10^9 kg/s)', 'FontSize',16)

%% Salt Transport in Northern Box

% Extract individual contributions in the northern box
advn1_1 = anal1(:, 18);
difn1_1 = anal1(:, 19);
emp1_1 = anal1(:, 20);

advn1_2 = anal2(:, 18);
difn1_2 = anal2(:, 19);
emp1_2 = anal2(:, 20);


figure()

% Plot time derivative of salt content in northern box
plot(time_i, ft*gradient(Snorth1),'-','LineWidth',2, 'Color',"#77AC30", 'DisplayName', ...
     ['$\frac{d}{dt} \int_{north box} S \, dV $']), hold on
plot(time_j, ft*gradient(Snorth2),'--','LineWidth',2, 'Color',"#77AC30", ...
     'HandleVisibility', 'off')

% Plot E-P contribution
plot(time_i, ft*emp1_1,'-','LineWidth',2,'Color','#7E2F8E', 'DisplayName', '$\Phi^s$') 
plot(time_j, ft*emp1_2,'--','LineWidth',2,'Color','#7E2F8E', 'HandleVisibility', 'off') 

% Compute and plot net lateral flux into the north box
latfluxnorth_1 = -ft*(advn1 - difn1 - (advn1_1 - difn1_1));
latfluxnorth_2 = -ft*(anal2(:,12) - anal2(:,14) - (advn1_2 - difn1_2));

plot(time_i, latfluxnorth_1, '-','LineWidth',2,'Color',[0 0 1], 'DisplayName','$\Phi^{lat}$')
plot(time_j, latfluxnorth_2,'--','LineWidth',2,'Color',[0 0 1], 'Handlevisibility', 'off') 
hold on

% Plot residual in the salt balance
plot(time_i, -latfluxnorth_1 - ft*emp1_1 + ft*gradient(Snorth1),'-','LineWidth',2,'Color',[0 0 0], 'Displayname', '$\Phi^b$') 
plot(time_j, -latfluxnorth_2 - ft*emp1_2 + ft*gradient(Snorth2),'--','LineWidth',2,'Color',[0 0 0], 'Handlevisibility', 'off') 
hold on

% Add reference lines and formatting
yline(0, 'HandleVisibility', 'off')
xline(t1(i), 'LineWidth',1,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',1,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',1,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',1,'HandleVisibility', 'off')
xlim([0 max([time_i; time_j])])

title('Salt balance vs time - north', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 4;
xlabel('time (years)', 'FontSize',16)
ylabel('salt transport (10^9 kg/s)', 'FontSize',16)

%% Advective salt inflow at northern boundary

figure()

% Advective salt inflow at northern boundary
plot(time_i, ft*anal1(:,12),'-','LineWidth',1,'Color','#4DBEEE', 'DisplayName','$\Phi^a(\theta_n)$'), hold on
plot(time_j, ft*anal2(:,12),'--','LineWidth',1,'Color','#4DBEEE', 'HandleVisibility', 'off') 

% Advective salt inflow at southern boundary of the box
plot(time_i, ft*advn1_1,'-','LineWidth',1,'Color',[1 0 0], 'DisplayName','$\Phi^a(\theta_{s})$' )
plot(time_j, ft*advn1_2,'--','LineWidth',1,'Color',[1 0 0], 'HandleVisibility', 'off' )

% Diffusive salt flux at northern boundary
plot(time_i, ft*anal1(:,14),'-','LineWidth',1,'Color','#77AC30', 'DisplayName', '$\Phi^d(\theta_n)$')
plot(time_j, ft*anal2(:,14),'--','LineWidth',1,'Color','#77AC30', 'HandleVisibility', 'off') 
hold on

% Diffusive salt flux at southern boundary of the box
plot(time_i, ft*difn1_1,'-','LineWidth',1,'Color',[0 0 1], 'DisplayName', '$\Phi^d(\theta_s)$') 
plot(time_j, ft*difn1_2,'--','LineWidth',1,'Color',[0 0 1], 'HandleVisibility', 'off') 
hold on

% Plot total lateral flux for reference
latfluxnorth_1 = -ft*(advn1 - difn1 - (advn1_1 - difn1_1));
latfluxnorth_2 = -ft*(anal2(:,12) - anal2(:,14) - (advn1_2 - difn1_2));

plot(time_i, latfluxnorth_1, '-','LineWidth',2,'Color',[0 0 1], 'DisplayName', ...
    '$\Phi^{lat}$') 
plot(time_j, latfluxnorth_2,'--','LineWidth',2,'Color',[0 0 1], 'Handlevisibility', 'off') 

% Add reference lines and formatting
yline(0, 'HandleVisibility', 'off')
xline(t1(i), 'LineWidth',1,'HandleVisibility', 'off')
xline(t2(i), 'LineWidth',1,'HandleVisibility', 'off')
xline(t1(j), '--', 'LineWidth',1,'HandleVisibility', 'off')
xline(t2(j), '--', 'LineWidth',1,'HandleVisibility', 'off')
xlim([0 max([time_i; time_j])])

title('Lateral salt fluxes vs time - north', 'FontSize',16)
lgd = legend('FontSize',16,'Location','southoutside', 'Interpreter', 'latex', 'Orientation','horizontal');
lgd.NumColumns = 3;
xlabel('time (years)','FontSize',16)
ylabel('salt transport (10^9 kg/s)','FontSize',16)
