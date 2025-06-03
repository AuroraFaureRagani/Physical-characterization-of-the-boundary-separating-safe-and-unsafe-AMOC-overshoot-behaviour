%% --- ANALYTICAL RESULTS ---
% Creates panels a-b in Figure 6

% This script plots the forcing and analytical solutions for different values
% of the parameters alpha and beta. The solution is constructed in two parts
% corresponding to different regimes of the forcing function.

% Clear workspace
clearvars
close all

% Number of solution cases to plot
n = 2;

% Parameters for each case
beta_vec = [10, 8];         % Slope of forcing after t1
alpha_vec = [4, 4];         % Slope of forcing before t1
t1_vec = [1.8, 1.8];        % Switching time from alpha to beta slope

% Global parameters for forcing
t0 = -2;                    % Start time of integration
t_tip = 1;                  % Tipping point time (vertex of first parabola)
t_end = 4;                  % End time of integration

% Colors and styles for plotting
colors = lines(n);
linestyles = {'-', '--'};

% Loop through each parameter set
for i = 1:n

    % Assign alpha, beta, t1 for current simulation
    alpha = alpha_vec(i);
    beta = beta_vec(i);
    t1 = t1_vec(i);

    color = colors(i, :);             % Line color
    linestyle = linestyles{i};        % Line style

    % Define time intervals for forcing
    time1 = t0:0.01:t1;               % Before switching forcing
    time2 = t1:0.01:t_end;            % After switching forcing
    time = [time1, time2];            % Full time span

    %% Define piecewise linear forcing function:
    % First phase: linear ramp up to t1
    forcing1 = @(t) -alpha.*(t_tip - t);
    % Value of forcing at transition
    max_forcing = forcing1(t1);
    % Intercept value to make forcing continuous at t1
    k = alpha*(t1 - t_tip) + beta*(t1 - 1);
    % Second phase: linear decay after t1
    forcing2 = @(t) k - beta.*(t - 1);
    % Full piecewise forcing function
    forcing = @(t) (heaviside(t1 - t) .* forcing1(t) + heaviside(t - t1) .* forcing2(t));

    %% Plot forcing
    figure(1)
    plot(time1, forcing(time1), 'Color', color, 'LineWidth', 2, 'LineStyle', linestyle, 'HandleVisibility', 'off'); hold on
    plot(time2, forcing(time2), 'Color', color, 'LineWidth', 2, ...
        'DisplayName', "$\alpha = $" + alpha + ", $\beta = $" + beta + ", $ max(f) = $" + max_forcing);
    xlabel('time', 'FontSize', 16)
    ylabel('f', 'FontSize', 16)
    yline(0, 'LineWidth', 1.2, 'HandleVisibility', 'off')
    title('Forcing', 'FontSize', 16)
    grid on
    ylim([forcing1(t0), forcing1(t1) + 3])
    lgd = legend('FontSize', 16, 'Location', 'southoutside', 'Interpreter', 'latex', 'Orientation', 'horizontal');
    lgd.NumColumns = 1;

    %% Plot fixed point curves (±√(-forcing))
    figure(2)
    valid1 = forcing1(time) <= 0;    % Real domain for sqrt
    valid2 = forcing2(time) <= 0.01; % Real domain for sqrt

    % Plot parabola from forcing1
    plot(time(valid1),  sqrt(-forcing1(time(valid1))), ...
         time(valid1), -sqrt(-forcing1(time(valid1))), ...
         'Color', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off'); hold on

    % Plot parabola from forcing2
    plot(time(valid2),  sqrt(-forcing2(time(valid2))), ...
         time(valid2), -sqrt(-forcing2(time(valid2))), ...
         'Color', color, 'LineWidth', 1.2, 'HandleVisibility', 'off'); hold on

    %% Airy functions and solution in first forcing regime (sol 1)
    % First interval (t0 to t1) solution (from Appendix)
    z = -forcing1(time1) ./ alpha^(2/3);           % Rescaled argument for Airy
    z_bar = -forcing1(t0) ./ alpha^(2/3);          % Initial rescaled forcing
    x0 = sqrt(-forcing1(t0));                      % Initial condition x(t0)

    % Define Airy functions and derivatives
    A   = @(y) airy(0, y); Ap = @(y) airy(1, y);
    B   = @(y) airy(2, y); Bp = @(y) airy(3, y);

    % Coefficients of solution (from initial conditions)
    C1 = (Bp(z_bar) + x0 / alpha^(1/3) * B(z_bar)) / ...
         (A(z_bar) * Bp(z_bar) - B(z_bar) * Ap(z_bar));
    C2 = (-x0 / alpha^(1/3) * A(z_bar) - Ap(z_bar)) / ...
         (A(z_bar) * Bp(z_bar) - B(z_bar) * Ap(z_bar));

    % Solution 1
    sol1 = @(y) -(alpha^(1/3)) .* (C1 * Ap(y) + C2 * Bp(y)) ./ (C1 * A(y) + C2 * B(y));
    % Define tipping time
    t_star = 1 + 2.388 * alpha^(-1/3);

    % Plot first part of the solution
    plot(time1, sol1(z), 'LineWidth', 2, 'Color', color, 'LineStyle', linestyle, 'HandleVisibility', 'off'); hold on

    %% Solution in first forcing regime (sol 2)
    x1 = sol1(-forcing1(t1) / alpha^(2/3));        % Match with end of sol1
    w = -forcing2(time2) / beta^(2/3);             % Rescaled argument
    w_bar = -forcing2(t1) / beta^(2/3);            % Rescaled transition value

    % Coefficients for solution 2 (from continuity at t1)
    K1 = (Bp(w_bar) - x1 / beta^(1/3) * B(w_bar)) / ...
         (A(w_bar) * Bp(w_bar) - B(w_bar) * Ap(w_bar));
    K2 = (x1 / beta^(1/3) * A(w_bar) - Ap(w_bar)) / ...
         (A(w_bar) * Bp(w_bar) - B(w_bar) * Ap(w_bar));

    % Solution 2
    sol2 = @(y) (beta^(1/3)) .* (K1 * Ap(y) + K2 * Bp(y)) ./ (K1 * A(y) + K2 * B(y));

    % Plot second part of the solution
    plot(time2, sol2(w), 'LineWidth', 2, 'Color', color, ...
        'DisplayName', "$\alpha = $" + alpha + ", $\beta = $" + beta + ", $max(f) = $" + max_forcing);

    % Format plot
    grid on
    ylim([-10 10])
    xlim([t0, time2(end)])
    xlabel('time', 'FontSize', 16)
    ylabel('x', 'FontSize', 16)
    title('Analytical solution', 'FontSize', 16)
    lgd = legend('FontSize', 16, 'Location', 'southoutside', 'Interpreter', 'latex', 'Orientation', 'horizontal');
    lgd.NumColumns = 1;
end
