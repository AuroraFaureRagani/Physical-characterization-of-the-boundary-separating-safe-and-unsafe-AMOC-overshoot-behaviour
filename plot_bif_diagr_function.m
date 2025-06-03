hold on

% First saddle node
x_pos1 = -0.24550E+01*factatl;
y_pos1 = 0.53833E+01; 
scatter(x_pos1, y_pos1, "filled", 'o', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')

% Second saddle node
x_pos2 = -0.68801E+00*factatl;
y_pos2 = 0.79277E+00; 
scatter(x_pos2, y_pos2, "filled", 'o', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off')

% Upper stable branch (before the first saddle node) - solid line
stable_upper = bif_diagr(bif_diagr(:,5) >= y_pos1, :);
plot([stable_upper(:, 3)*factatl; x_pos1], [stable_upper(:, 5); y_pos1], '-', 'Color', 'k', 'HandleVisibility', 'off', 'LineWidth', 1.2)

% Unstable middle branch - dashed line
unstable_middle = bif_diagr((bif_diagr(:,5) < y_pos1) & (bif_diagr(:,5) >= y_pos2), :);
plot([unstable_middle(:, 3)*factatl; x_pos2], [unstable_middle(:, 5); y_pos2], '--', 'Color', 'k', 'HandleVisibility', 'off', 'LineWidth', 1.2)

% Lower stable branch (at zero) - solid line
stable_lower = bif_diagr(bif_diagr(:,5) < y_pos2, :);
plot([stable_lower(:, 3)*factatl; linspace(bif_diagr(end, 3)*factatl, 0.3, 100)'], [stable_lower(:, 5); zeros(100,1)], '-', 'Color', 'k', 'HandleVisibility', 'off', 'LineWidth', 1.2)
