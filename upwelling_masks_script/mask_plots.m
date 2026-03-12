clear; clc; close all;

% ========================================================================
% PLOT THREE OFFSHORE MASKS ON GEBCO BACKGROUND - CORRECTED VERSION
% ========================================================================
% The mask regions will now be clearly visible as solid grey areas

save_path_plot = 'Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\plot\plot_mask';
save_mask = 'Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\upwelling_masks_script';

if ~isfolder(save_path_plot)
    mkdir(save_path_plot);
end

% Load GEBCO data
fprintf('Loading GEBCO data...\n');
[lon_gebco, lat_gebco, elev_gebco] = open_gebco;

fprintf('GEBCO data loaded:\n');
fprintf('  Lon range: %.2f to %.2f\n', min(lon_gebco), max(lon_gebco));
fprintf('  Lat range: %.2f to %.2f\n', min(lat_gebco), max(lat_gebco));
fprintf('  Elevation grid size: %d x %d\n', size(elev_gebco, 1), size(elev_gebco, 2));

% ========================================================================
% LOAD PRE-COMPUTED MASKS AND DISTANCE FIELD
% ========================================================================

fprintf('\nLoading masks...\n');
mask_file = fullfile(save_mask, 'ccmp_offshore_masks.mat');

if isfile(mask_file)
    load(mask_file, 'mask_50km', 'mask_150km', 'mask_300km', ...
        'lon_sub', 'lat_sub', 'distance_offshore_km');
    fprintf('Masks loaded successfully\n');
    fprintf('Mask grid size: %d lon x %d lat\n', length(lon_sub), length(lat_sub));
else
    error('Mask file not found at: %s\n', mask_file);
end

% Create meshgrid for plotting
[LON_MASK, LAT_MASK] = meshgrid(lon_sub, lat_sub);

% ========================================================================
% FIGURE 1: 50 KM MASK ON GEBCO
% ========================================================================

fprintf('\n========================================\n');
fprintf('Creating mask visualizations...\n');
fprintf('========================================\n\n');

fprintf('Creating 50 km mask figure...\n');

figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as greyscale background
contourf(lon_gebco, lat_gebco, elev_gebco, 12, 'LineColor', 'none', 'HandleVisibility', 'off');
colormap(gca, flipud(gray));
hold on

% Plot coastline (0m contour) and 200m isobath from GEBCO
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 2, ...
    'DisplayName','Coastline (0 m)');

contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5, ...
    'DisplayName','Isobath (200 m)');

legend

% Convert mask to display format
mask_display = double(mask_50km);

% Create RGB image where mask=1 becomes grey and mask=0 becomes transparent
[ny, nx] = size(mask_display);

mask_rgb = ones(ny, nx, 3);

% Use surf to plot the mask with transparency
[X, Y] = meshgrid(lon_sub, lat_sub);
h_surf = surf(X, Y, zeros(size(X)), mask_rgb, 'EdgeColor', 'none', 'HandleVisibility','off');

% Set transparency based on mask values
alpha_matrix = zeros(ny, nx);
alpha_matrix(mask_display == 1) = 0.5;  % 50% opacity for masked cells
alpha_matrix(mask_display == 0) = 0;    % Fully transparent for non-masked cells
set(h_surf, 'AlphaData', alpha_matrix);

view(2);  % Top-down 2D view
set(h_surf, 'FaceAlpha', 'flat');

% Overlay distance contours
[C0, h0] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 0], 'k', 'LineWidth', 2.5, ...
       'DisplayName', 'Parallel line to coast (0 km offshore)');
clabel(C0, h0, 'FontSize', 13, 'Color', 'black', 'FontWeight', 'bold');

[C50, h50] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [50 50], 'Color', [0.6 0 0], 'LineWidth', 2.5, ...
    'DisplayName', 'Parallel line to coast (50 km offshore)');
clabel(C50, h50, 'FontSize', 13, 'Color', [0.6 0 0], 'FontWeight', 'bold');


xlabel('Longitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Latitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
title('Offshore Mask: 0 to 50 km (mask = red shaded)', ...
    'FontSize', 22, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on
set(gca,'GridColor',[0.7 0.7 0.7], ...   % light grey grid
        'GridAlpha',0.3, ...             % transparency
        'MinorGridAlpha',0.15)

xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 2;
ax.Layer = 'top';
box on
legend('Location','northwest')
saveas(gcf, fullfile(save_path_plot, 'mask_50km_on_gebco.png'), 'png');
savefig(fullfile(save_path_plot,'mask_50km_on_gebco.fig'));
fprintf('  ✓ Figure saved: mask_50km_on_gebco.png\n');
close(gcf);

% ========================================================================
% FIGURE 2: 150 KM MASK ON GEBCO
% ========================================================================

fprintf('Creating 150 km mask figure...\n');

figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as greyscale background
contourf(lon_gebco, lat_gebco, elev_gebco, 12, 'LineColor','none', 'HandleVisibility', 'off');
colormap(gca, flipud(gray));
hold on

% Plot coastline and 200m isobath
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 2, ...
    'DisplayName','Coastline (0 m)');

contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5, ...
    'DisplayName','Isobath (200 m)');

legend

% Convert mask to display format
mask_display = double(mask_150km);

% Create RGB image
[ny, nx] = size(mask_display);
mask_rgb = ones(ny, nx, 3);


% Use surf to plot the mask with transparency
[X, Y] = meshgrid(lon_sub, lat_sub);
h_surf = surf(X, Y, zeros(size(X)), mask_rgb, 'EdgeColor', 'none', 'HandleVisibility','off');

% Set transparency based on mask values
alpha_matrix = zeros(ny, nx);
alpha_matrix(mask_display == 1) = 0.5;  % 50% opacity for masked cells
alpha_matrix(mask_display == 0) = 0;    % Fully transparent for non-masked cells
set(h_surf, 'AlphaData', alpha_matrix);

view(2);  % Top-down 2D view
set(h_surf, 'FaceAlpha', 'flat');

% Overlay distance contours
[C0, h0] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 0], 'k', 'LineWidth', 2.5, ...
     'DisplayName', 'Parallel line to coast (0 km offshore)');
clabel(C0, h0, 'FontSize', 13, 'Color', 'black', 'FontWeight', 'bold');

[C150, h150] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [150 150], 'Color', [0.6 0 0], 'LineWidth', 2.5, ...
     'DisplayName', 'Parallel line to coast (150 km offshore)');
clabel(C150, h150, 'FontSize', 13, 'Color', [0.6 0 0], 'FontWeight', 'bold');


xlabel('Longitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Latitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
title('Offshore Mask: 0 to 150 km (mask = red shaded)', ...
    'FontSize', 22, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on
set(gca,'GridColor',[0.7 0.7 0.7], ...   % light grey grid
        'GridAlpha',0.3, ...             % transparency
        'MinorGridAlpha',0.15)

xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 2;
ax.Layer = 'top';
box on
legend('Location','northwest')
saveas(gcf, fullfile(save_path_plot, 'mask_150km_on_gebco.png'), 'png');
savefig(fullfile(save_path_plot,'mask_150km_on_gebco.fig'));
fprintf('  ✓ Figure saved: mask_150km_on_gebco.png\n');
close(gcf);

% ========================================================================
% FIGURE 3: 300 KM MASK ON GEBCO
% ========================================================================

fprintf('Creating 300 km mask figure...\n');

figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as greyscale background
contourf(lon_gebco, lat_gebco, elev_gebco, 12, 'LineColor','none',  'HandleVisibility', 'off');
colormap(gca, flipud(gray));
hold on

% Plot coastline and 200m isobat
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 2, ...
    'DisplayName','Coastline (0 m)');

contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5, ...
    'DisplayName','Isobath (200 m)');

legend
% Convert mask to display format
mask_display = double(mask_300km);

% Create RGB image
[ny, nx] = size(mask_display);
mask_rgb = ones(ny, nx, 3);


% Use surf to plot the mask with transparency
[X, Y] = meshgrid(lon_sub, lat_sub);
h_surf = surf(X, Y, zeros(size(X)), mask_rgb, 'EdgeColor', 'none', 'HandleVisibility','off');

% Set transparency based on mask values
alpha_matrix = zeros(ny, nx);
alpha_matrix(mask_display == 1) = 0.5;  % 50% opacity for masked cells
alpha_matrix(mask_display == 0) = 0;    % Fully transparent for non-masked cells
set(h_surf, 'AlphaData', alpha_matrix);

view(2);  % Top-down 2D view
set(h_surf, 'FaceAlpha', 'flat');

% Overlay distance contours
[C0, h0] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 0], 'k', 'LineWidth', 2.5, ...
     'DisplayName', 'Parallel line to coast (0 km offshore)');
clabel(C0, h0, 'FontSize', 13, 'Color', 'black', 'FontWeight', 'bold');

[C300, h300] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [300 300], 'Color', [0.6 0 0], 'LineWidth', 2.5, ...
     'DisplayName', 'Parallel line to coast (300 km offshore)');
clabel(C300, h300, 'FontSize', 13, 'Color', [0.6 0 0], 'FontWeight', 'bold');

xlabel('Longitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Latitude (°)', 'FontSize', 20, 'FontWeight', 'bold');
title('Offshore Mask: 0 to 300 km (mask = red shaded)', ...
    'FontSize', 22, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on
set(gca,'GridColor',[0.7 0.7 0.7], ...   % light grey grid
        'GridAlpha',0.3, ...             % transparency
        'MinorGridAlpha',0.15)

xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

ax = gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1;
ax.Layer = 'top';
box on
legend('Location','northwest')
saveas(gcf, fullfile(save_path_plot, 'mask_300km_on_gebco.png'), 'png');
savefig(fullfile(save_path_plot,'mask_300km_on_gebco.fig'));
fprintf('  ✓ Figure saved: mask_300km_on_gebco.png\n');
close(gcf);

fprintf('\n========================================\n');
fprintf('✓ All three mask figures created!\n');
fprintf('Saved to: %s\n', save_path_plot);
fprintf('========================================\n');