%% --- ANALYTICAL RESULTS ---
%% critical curve for collapse conditions
% Creates panels a-d of Figure 7

% Clear workspace
clearvars
close all

% Initial condition parameters
t0 = 0.8476;      % Initial time
x0 = 3.1518;      % Initial state value at time t0
t_tip = 1;        % Time at which tipping starts
alpha = 84.0607;  % Slope of the forcing before the switching point t1

fig = 'a';        % Panel from the paper to reproduce ('a', 'b', 'c', or 'd')
                  % If 'a' or 'b' is selected, both figures will be generated

% Choose parameter ranges based on selected figure panel
if (fig == 'a' || fig == 'b')
    beta_vec = [92, 88, 84, 80, 76];     % Slopes for the forcing after t1
    t1_vec = 1.22:0.01:1.23;             % Switching times
elseif(fig == 'c')
    beta_vec = 110:-0.05:20;
    t1_vec = 1.14:0.02:1.24;
elseif(fig == 'd')
    beta_vec = 250:-0.05:10;
    t1_vec = 1.10:0.01:1.3;
end

% Variables for tracking collapse behavior
beta_tilde  = [];      % Store maximum beta values that lead to collapse
collapse    = 0;       % Boolean to flag if collapse has occurred
k_col       = 0;       % Index to select color from colormap
colors = lines(numel(t1_vec));  % Generate distinguishable colors for plotting

% Initialize progress bar
f = waitbar(0, '0.0 %', 'Name', 'Creating figure', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');

% Main loop over t1 values
for t1 = t1_vec
    tstar = [];       % Collapse times
    t1star = [];      % t1 values leading to collapse
    betastar = [];    % beta values that lead to collapse

    % Update progress bar
    waitbar((t1 - t1_vec(1)) / (t1_vec(end) - t1_vec(1)), f, ...
        sprintf('%.1f %%', (t1 - t1_vec(1)) / (t1_vec(end) - t1_vec(1)) * 100))

    % Stop if user cancels progress bar
    if getappdata(f, 'canceling')
        break
    end

    k_col = k_col + 1;  
    if (fig == 'a' || fig == 'b')
        color = colors(k_col, :);  % Select color for current t1
    end

    % Loop over beta values
    for beta = beta_vec
        if getappdata(f, 'canceling')
            break
        end
        collapse = 0;

        % Define time intervals before and after switching time
        time1 = t0:0.001:t1;
        time2 = t1:0.001:3;
        time = [time1, time2];

        % Define piecewise linear forcing
        forcing1 = @(t) -alpha * (t_tip - t);                     % Before t1
        k = alpha * (t1 - t_tip) + beta * (t1 - 1);               % Ensure continuity
        forcing2 = @(t) k - beta * (t - 1);                       % After t1
        forcing = @(t) (heaviside(t1 - t).*forcing1(t) + ...
                        heaviside(t - t1).*forcing2(t));         % Combined forcing

        % Plot forcing for figures a and b
        if (fig == 'a' || fig == 'b')
            figure(1)
            plot(time, forcing(time), 'color', color, 'LineStyle', '-', 'LineWidth', 1.2)
            hold on
            grid on
            title('Forcing', 'FontSize', 16)
            xlabel('\tau', 'FontSize', 16)
            ylabel('f(\tau)', 'FontSize', 16)
            xlim([1.14, 1.3])
        end

        % Define Airy functions and their derivatives
        A   = @(y) airy(0, y);
        Ap  = @(y) airy(1, y);
        B   = @(y) airy(2, y);
        Bp  = @(y) airy(3, y);

        % --- Solution before t1 (from initial condition) ---
        z = -forcing1(time1) / alpha^(2/3);
        z_bar = -forcing1(t0) / alpha^(2/3);  % Value at initial time

        % Match initial condition using Airy function solution
        C1 = (Bp(z_bar) + x0 / alpha^(1/3) * B(z_bar)) / (A(z_bar) * Bp(z_bar) - B(z_bar) * Ap(z_bar));
        C2 = (-x0 / alpha^(1/3) * A(z_bar) - Ap(z_bar)) / (A(z_bar) * Bp(z_bar) - B(z_bar) * Ap(z_bar));
        sol1 = @(y) -(alpha^(1/3)) * (C1 * Ap(y) + C2 * Bp(y)) ./ (C1 * A(y) + C2 * B(y));

        % --- Solution after t1 ---
        x1 = sol1(-forcing1(t1) / alpha^(2/3));  % Value at switching point
        w = -forcing2(time2) / beta^(2/3);
        w_bar = -forcing2(t1) / beta^(2/3);

        % Match continuity condition using Airy functions
        K1 = (Bp(w_bar) - x1 / beta^(1/3) * B(w_bar)) / (A(w_bar) * Bp(w_bar) - B(w_bar) * Ap(w_bar));
        K2 = (x1 / beta^(1/3) * A(w_bar) - Ap(w_bar)) / (A(w_bar) * Bp(w_bar) - B(w_bar) * Ap(w_bar));
        sol2 = @(y) (beta^(1/3)) * (K1 * Ap(y) + K2 * Bp(y)) ./ (K1 * A(y) + K2 * B(y));

        % Plot the solution if figure is 'a' or 'b'
        if (fig == 'a' || fig == 'b')
            figure(2)
            plot(time1, sol1(z), 'LineWidth', 1.2, 'Color', [0.4660 0.6740 0.1880])  % Before t1
            hold on
            plot(time2, sol2(w), 'LineWidth', 1.2, 'Color', color, 'LineStyle', '-') % After t1
            plot(t1, x1, 'o')  % Mark matching point
            xlabel('\tau', 'FontSize', 16)
            ylabel('x', 'FontSize', 16)
            ylim([-10, 5])
            xlim([t0, 2])
            title('Solution for different \beta and \tau_1', 'FontSize', 16)
            grid on
        end

        % --- Collapse check: if solution dips below threshold ---
        [m, i_min] = min([sol1(z), sol2(w)]);
        if (m < -8)
            t1star = [t1star, t1];
            tstar = [tstar, time(i_min)];
            betastar = [betastar, beta];
            collapse = 1;
        end
    end

    % Plot tipping time vs beta (for panel c)
    if(fig == 'c')
        figure(3)
        plot(betastar, tstar, 'LineWidth', 2, 'DisplayName', "\tau_1 = " + t1)
        hold on
    end

    % Store max beta that leads to collapse (for panel d)
    beta_tilde = [beta_tilde, max(betastar)];
end

%% Final plots for critical boundaries (if user didn’t cancel)

if isempty(getappdata(f, 'canceling'))
    if(fig == 'c')
        figure(3)
        grid on
        xlim([beta_vec(end)-1, beta_vec(1)+1])
        ylim([1.5, 2.1])
        xlabel('\beta', 'FontSize', 16)
        ylabel('\tau^*', 'FontSize', 16)
        legend('FontSize', 16, 'Location','eastoutside')
        title('\tau^* vs \beta', 'FontSize', 16)
    end

    if(fig == 'd')
        figure(4)
        plot(t1_vec, beta_tilde, '-o', ...
            'LineWidth', 2, ...
            'MarkerSize', 3, ...
            'MarkerEdgeColor', 'b', ...
            'MarkerFaceColor', 'b', ...
            'DisplayName', "\alpha = " + alpha)
        hold on
        grid on
        xlabel('\tau_1', 'FontSize', 16)
        ylabel('$\tilde{\beta}$', 'Interpreter', 'Latex', 'FontSize', 16)
        title('$\tilde{\beta}$ vs $\tau_1$', 'Interpreter', 'Latex', 'FontSize', 16)
        text(1.22, 40, 'unsafe', 'FontSize', 16)
        text(1.17, 100, 'safe', 'FontSize', 16)
    end
end

% Close progress bar
delete(f)
