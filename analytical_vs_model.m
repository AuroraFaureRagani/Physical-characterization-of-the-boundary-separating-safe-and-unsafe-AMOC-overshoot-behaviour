%% --- COMPARISON of analytical results and global ocean model simulations ---
% Creates panels c-f of Figure 8
%% Comparison of analytical results and global ocean model simulations

% Clear workspace
clearvars
close all

% Select simulation case
simulation = 'd';
switch simulation
    case 'a'
        sim_PL = 3;
    case 'b'
        sim_PL = 4;
    case 'c'
        sim_PL = 5;
    case 'd'
        sim_PL = 6;
end

% Load data: simulation time series and bifurcation diagram
fname_PL = ["data/p.TPL" + string(sim_PL)];
xx = dlmread(fname_PL);            % Time series from simulation
bif_diagr = dlmread("data/p.AUR1");% Bifurcation diagram

color = [0.0000    0.4470    0.7410]; % Plot color

%% Define constants

udim  = 0.1;                 % [m/s] velocity scale
r0dim = 6.4e6;               % [m] radius of Earth
T0    = 15;                  % [°C] reference temperature
S0    = 35;                  % [psu] reference salinity
RtD   = 180/pi;              % radians to degrees
factatl     = -0.913*8.308e-02;             % scaling factor for AMOC strength
factrate    = 2.64e-03;                     % [10^-4 Sv/yr] scaling factor
facttime    = r0dim/udim/(360*24*60*60);    % [yr] dimensional time scale

% Extract saddle-node bifurcation point from bifurcation diagram
j = min(find(bif_diagr(:,1)==2));
threshold = [bif_diagr(j,3)*factatl, bif_diagr(j,5)]; % [freshwater flux, AMOC_strength]

%% Read tipping time and amplitude from Excel summary of simulations

table_sim = xlsread('data/sim_PL.xlsx');

t1 = table_sim(sim_PL,2);           % Time when forcing switch (from increasing to decreasing)
t2 = table_sim(sim_PL,3);           % Time when forcing ends
gamma_max = table_sim(sim_PL,4);    % Max forcing amplitude

% Define forcing profile gamma(t) with three segments: ramp-up, ramp-down, then constant
gamma_lim = threshold(1)/2;
dgamma_PL = @(t1, t2, ymax, t) -( ...
    heaviside(t1-t).*((ymax/t1)*t) + ...
    heaviside(t-t1).*heaviside(t2-t).*(((threshold(1)/2-ymax)/(t2-t1)).*t + ymax-t1.*((threshold(1)/2-ymax)/(t2-t1))) + ...
    heaviside(t-t2).*threshold(1)/2 ...
) / factatl;

%% Read domain and grid mask from fort.44 (model grid)

[n, m, l, la, nun, xmin, xmax, ymin, ymax, hdim, x, y, z, xu, yv, zw, landm] = readfort44('fort.44');
surfm = landm(2:n+1,2:m+1,l+1);  % Surface land mask (interior points only)
dx = (xu(n+1)-xu(1))/n;
dy = (yv(m+1)-yv(1))/m;
dz = (zw(l+1)-zw(1))/l;
[qz,dfzt,dfzw] = gridstretch(zw);  % Vertical grid stretching

%% Fit the saddle-node region of the bifurcation diagram and define rescaled parameters

y_sc = bif_diagr(13:15, 3)*factatl;     % AMOC strength
x_sc = bif_diagr(13:15, 5);             % Freshwater flux
p = polyfit(x_sc, y_sc, 2);             % Fit parabola (SN bifurcation)
a = p(1); b = p(2); c = p(3);


% Define tipping time scale tau_sn based on parabolic fit
m1 = gamma_max/t1;
m2 = (gamma_max - threshold(1)/2)/(t2 - t1);
tau_sn = c/m1 - b^2/(4*a*m1);

% Define forcing parameters in scaled time
t0 = t1 - 100;
t_end = xx(end, 3);
time1 = [t0:0.1:t1]./tau_sn;
time2 = [t1:0.1:t_end]./tau_sn;
time = [time1, time2];

% Define coefficients of piecewise forcing function f(τ)
alpha = -m1*a*tau_sn^3;
beta = -m2*a*tau_sn^3*1;
h = (m1 + m2)*t1;
k = (alpha + beta)*(t1/tau_sn - 1); 

% Define the piecewise forcing function
forcing1 = @(t) -alpha.*(1-t);
forcing2 = @(t) +k+beta.*(1-t);
forcing = @(t) ...
    (heaviside(t1/tau_sn - t).*forcing1(t) + ...
    heaviside(t - t1/tau_sn).*heaviside(t2/tau_sn - t).*forcing2(t));

% Determine initial condition x0
jj = max(find(xx(:,3)<t0));
if (xx(1,3)==1)
    if isempty(jj); jj = 1; end
    x0 = -a*(xx(jj,5)+b/(2*a))*(tau_sn);
else
    x0 = sqrt(-forcing1(t0/tau_sn));
end

% Mask invalid values (where sqrt is imaginary)
valid1 = forcing1(time1) <= 0.05;
valid2 = forcing2(time2) <= 0.05;

% Plot forcing envelope (horizontal parabolas) and simulation
figure()
plot(time1(valid1), sqrt(-forcing1(time1(valid1))), 'k', ...
     time1(valid1), -sqrt(-forcing1(time1(valid1))), 'k', ...
     'HandleVisibility', 'off', 'LineWidth', 1.2); hold on
plot(time2(valid2), sqrt(-forcing2(time2(valid2))), 'k', ...
     time2(valid2), -sqrt(-forcing2(time2(valid2))), 'k', ...
     'HandleVisibility', 'off', 'LineWidth', 1.2);
plot(xx(:,3)./tau_sn, -a*(xx(:,5)+b/(2*a)).*tau_sn, '--', ...
     'DisplayName', "sim " + simulation, 'LineWidth', 2, 'color', color);
grid on

%% Airy function and solution in first forcing regime (sol 1)

z = -forcing1(time1)./alpha^(2/3);
z_bar = -forcing1(t0/tau_sn)./alpha^(2/3);

% Define Airy functions and derivatives
A   = @(y) airy(0,y);
Ap  = @(y) airy(1,y);
B   = @(y) airy(2,y);
Bp  = @(y) airy(3,y);

% Coefficients of solution (matching initial conditions)
C1 = (Bp(z_bar) + x0/(alpha^(1/3))*B(z_bar))/(A(z_bar)*Bp(z_bar) - B(z_bar)*Ap(z_bar));
C2 = (-x0/alpha^(1/3)*A(z_bar)-Ap(z_bar))/(A(z_bar)*Bp(z_bar) - B(z_bar)*Ap(z_bar));

% Solution 1
sol1 = @(y) -(alpha^(1/3)).*(C1*Ap(y) + C2*Bp(y))./(C1*A(y) + C2*B(y));

% Plot analytical solution (region 1)
plot(time1, sol1(z), 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880], 'DisplayName', "sol A");

%% Solution in second forcing regime (sol 2)

% Compute value at transition point
t_switch = t1/tau_sn;
x1 = sol1(-forcing1(t_switch)./alpha^(2/3));
plot(t_switch, x1, 'o', 'HandleVisibility', 'off')
ylim([-10,10])

% Define time vector and scaled forcing for second regime
time2_update = t_switch:0.01:t_end;
w       = -forcing2(time2_update)./beta^(2/3);
w_bar   = -forcing2(t_switch)/beta^(2/3);

% Coefficients for solution 2
K1 = (Bp(w_bar) - x1/beta^(1/3)*B(w_bar))/(A(w_bar)*Bp(w_bar) - B(w_bar)*Ap(w_bar));
K2 = (+x1/beta^(1/3)*A(w_bar)-Ap(w_bar))/(A(w_bar)*Bp(w_bar) - B(w_bar)*Ap(w_bar));

% Solution 2
sol2 = @(y) (beta^(1/3)).*(K1*Ap(y) + K2*Bp(y))./(K1*A(y) + K2*B(y));

% Plot analytical solution (region 2)
plot(time2_update, sol2(w),'LineWidth', 2, 'Color', [0.6350 0.0780 0.1840], 'DisplayName', "sol B");

grid on
xline(t1/tau_sn, 'HandleVisibility', 'off')
xlim([t0/tau_sn, t_end/tau_sn])
xlabel('\tau', 'FontSize', 16)
ylabel('x', 'FontSize', 16)
if (simulation == 'd')
    legend('Location', 'northeast', 'FontSize', 16)
else
    legend('Location', 'southeast', 'FontSize', 16)
end
title("Simulation " + simulation, 'FontSize', 16)
subtitle("\alpha = " + alpha + ", \tau_1 = " + t1/tau_sn + ", \beta = " + beta, 'FontSize', 16 )
