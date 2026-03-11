clear; clc; close all;

% ========================================================================
% PLOT THREE OFFSHORE MASKS ON GEBCO BACKGROUND
% ========================================================================

save_path_plot = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/plot/plot_mask';
save_mask = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/upwelling_masks_script';

if ~isfolder(save_path_plot)
    mkdir(save_path_plot);
end

% Load GEBCO data using your function
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

fprintf('\nCreating 50 km mask figure...\n');

fig1 = figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as background (inverted gray)
contourf(lon_gebco, lat_gebco, elev_gebco, 8, 'LineColor', 'none');
colormap(gca, flipud(gray));
hold on;

% Plot coastline (0m contour) and 200m isobath from GEBCO
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 1.5);      % 0 m = coast
contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5);  % 200 m isobath

% Overlay mask as semi-transparent grey patch
mask_50_display = double(mask_50km);

% Use patch objects for better transparency control
for i = 1:size(mask_50_display, 1)
    for j = 1:size(mask_50_display, 2)
        if mask_50_display(i, j) == 1
            % Draw small rectangle for each masked grid cell
            rect_lon = [lon_sub(j) - (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) - (lon_sub(2)-lon_sub(1))/2];
            rect_lat = [lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2];
            patch(rect_lon, rect_lat, [0.6 0.6 0.6], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4);
        end
    end
end

% Overlay distance contours (boundaries)
[C, h_contour] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 50], 'r', 'LineWidth', 2.5);
clabel(C, h_contour, 'FontSize', 11, 'Color', 'red', 'FontWeight', 'bold');

xlabel('Longitude (°)', 'FontSize', 12);
ylabel('Latitude (°)', 'FontSize', 12);
title('Offshore Mask: 0 to 50 km (grey shaded region)', ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on;
grid minor;

% Set axis limits to GEBCO extent
xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

saveas(fig1, fullfile(save_path_plot, 'mask_50km_on_gebco.png'));
fprintf('  ✓ Figure saved: mask_50km_on_gebco.png\n');
close(fig1);

% ========================================================================
% FIGURE 2: 150 KM MASK ON GEBCO
% ========================================================================

fprintf('Creating 150 km mask figure...\n');

fig2 = figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as background
contourf(lon_gebco, lat_gebco, elev_gebco, 8, 'LineColor', 'none');
colormap(gca, flipud(gray));
hold on;

% Plot coastline and 200m isobath
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 1.5);      % 0 m = coast
contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5);  % 200 m isobath

% Overlay mask as semi-transparent grey patch
mask_150_display = double(mask_150km);

for i = 1:size(mask_150_display, 1)
    for j = 1:size(mask_150_display, 2)
        if mask_150_display(i, j) == 1
            rect_lon = [lon_sub(j) - (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) - (lon_sub(2)-lon_sub(1))/2];
            rect_lat = [lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2];
            patch(rect_lon, rect_lat, [0.6 0.6 0.6], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4);
        end
    end
end

% Overlay distance contours
[C, h_contour] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 150], 'r', 'LineWidth', 2.5);
clabel(C, h_contour, 'FontSize', 11, 'Color', 'red', 'FontWeight', 'bold');

xlabel('Longitude (°)', 'FontSize', 12);
ylabel('Latitude (°)', 'FontSize', 12);
title('Offshore Mask: 0 to 150 km (grey shaded region)', ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on;
grid minor;

xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

saveas(fig2, fullfile(save_path_plot, 'mask_150km_on_gebco.png'));
fprintf('  ✓ Figure saved: mask_150km_on_gebco.png\n');
close(fig2);

% ========================================================================
% FIGURE 3: 300 KM MASK ON GEBCO
% ========================================================================

fprintf('Creating 300 km mask figure...\n');

fig3 = figure('Position', [100 100 1200 900], 'Color', 'w');

% Plot GEBCO elevation as background
contourf(lon_gebco, lat_gebco, elev_gebco, 8, 'LineColor', 'none');
colormap(gca, flipud(gray));
hold on;

% Plot coastline and 200m isobath
contour(lon_gebco, lat_gebco, elev_gebco, [0 0], 'k', 'LineWidth', 1.5);      % 0 m = coast
contour(lon_gebco, lat_gebco, elev_gebco, [-200 -200], 'b', 'LineWidth', 1.5);  % 200 m isobath

% Overlay mask as semi-transparent grey patch
mask_300_display = double(mask_300km);

for i = 1:size(mask_300_display, 1)
    for j = 1:size(mask_300_display, 2)
        if mask_300_display(i, j) == 1
            rect_lon = [lon_sub(j) - (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) + (lon_sub(2)-lon_sub(1))/2, ...
                       lon_sub(j) - (lon_sub(2)-lon_sub(1))/2];
            rect_lat = [lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) - (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2, ...
                       lat_sub(i) + (lat_sub(2)-lat_sub(1))/2];
            patch(rect_lon, rect_lat, [0.6 0.6 0.6], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4);
        end
    end
end

% Overlay distance contours
[C, h_contour] = contour(LON_MASK, LAT_MASK, distance_offshore_km, ...
    [0 300], 'r', 'LineWidth', 2.5);
clabel(C, h_contour, 'FontSize', 11, 'Color', 'red', 'FontWeight', 'bold');

xlabel('Longitude (°)', 'FontSize', 12);
ylabel('Latitude (°)', 'FontSize', 12);
title('Offshore Mask: 0 to 300 km (grey shaded region)', ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'YDir', 'normal');
grid on;
grid minor;

xlim([min(lon_gebco) max(lon_gebco)]);
ylim([min(lat_gebco) max(lat_gebco)]);

saveas(fig3, fullfile(save_path_plot, 'mask_300km_on_gebco.png'));
fprintf('  ✓ Figure saved: mask_300km_on_gebco.png\n');
close(fig3);

fprintf('\n========================================\n');
fprintf('✓ All three mask figures created!\n');
fprintf('Saved to: %s\n', save_path_plot);
fprintf('========================================\n');