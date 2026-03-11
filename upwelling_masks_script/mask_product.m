
clear; clc; close all;

% ========================================================================
% CREATE THREE OFFSHORE DISTANCE MASKS (50, 150, 300 km)
% ========================================================================

save_path_plot = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/plot/plot_mask';
save_mask= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/upwelling_masks_script';

if ~isfolder(save_path_plot)
    mkdir(save_path_plot);
end

% Full bounding box
LONMIN = -40;
LONMAX = 8;
LATMIN = 8;
LATMAX = 40;

fprintf('========================================\n');
fprintf('STEP 1: DEFINE THE TWO CAPES\n');
fprintf('========================================\n');

% Cape Bojador (northernmost)
cb_lat = 26.11;
cb_lon = -14.33;

% Cape Juby (southernmost)
cj_lat = 27.96;
cj_lon = -13.20;

fprintf('Cape Bojador: Lat = %.2f°N, Lon = %.2f°W\n', cb_lat, abs(cb_lon));
fprintf('Cape Juby:    Lat = %.2f°N, Lon = %.2f°W\n', cj_lat, abs(cj_lon));

% Convert to radians
rad = pi / 180;
cb_lat_rad = cb_lat * rad;
cb_lon_rad = cb_lon * rad;
cj_lat_rad = cj_lat * rad;
cj_lon_rad = cj_lon * rad;

fprintf('\n========================================\n');
fprintf('STEP 2: CALCULATE COASTLINE DIRECTION\n');
fprintf('========================================\n');

% Vector FROM Cape Bojador TO Cape Juby (this is the direction along the coast)
dlat_coast = cj_lat_rad - cb_lat_rad;
dlon_coast = (cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad + cb_lat_rad) / 2);

fprintf('Vector from Bojador to Juby:\n');
fprintf('  Δlat = %.6f rad = %.4f°\n', dlat_coast, dlat_coast/rad);
fprintf('  Δlon = %.6f rad = %.4f°\n', dlon_coast, dlon_coast/rad);

% Normalize coastline direction to unit vector
len_coast = sqrt(dlat_coast^2 + dlon_coast^2);
coast_dir_lat = dlat_coast / len_coast;
coast_dir_lon = dlon_coast / len_coast;

fprintf('Normalized coastline direction (unit vector):\n');
fprintf('  coast_dir_lat = %.6f (points north-ish)\n', coast_dir_lat);
fprintf('  coast_dir_lon = %.6f (points east-ish)\n', coast_dir_lon);

fprintf('\n========================================\n');
fprintf('STEP 3: CALCULATE OFFSHORE PERPENDICULAR\n');
fprintf('========================================\n');

% The OFFSHORE direction is perpendicular to the coastline
% pointing OUT TO SEA (westward into the Atlantic)
% 
% If coastline points [coast_lon, coast_lat], then perpendicular directions are:
%   Rotate +90° (CCW):  [-coast_lat, coast_lon]  = points left/northwest
%   Rotate -90° (CW):   [+coast_lat, -coast_lon] = points right/southeast
%
% Since coast points from Bojador(NW) to Juby(SE), rotating -90° CW gives SOUTHWEST
% But we want WEST/OFFSHORE. Let's use +90° CCW = NORTHWEST = AWAY FROM LAND

offshore_dir_lon = -coast_dir_lat;   % Perpendicular, pointing seaward
offshore_dir_lat = coast_dir_lon;

fprintf('Offshore direction (perpendicular to coast, +90° rotation):\n');
fprintf('  offshore_dir_lon = %.6f (should be NEGATIVE = westward)\n', offshore_dir_lon);
fprintf('  offshore_dir_lat = %.6f (small value)\n', offshore_dir_lat);

% Double-check the direction
angle_from_east = atan2(offshore_dir_lat, offshore_dir_lon) * 180/pi;
fprintf('  Direction angle: %.1f° from East (0° is E, 90° is N, 180° is W, -90° is S)\n', angle_from_east);
fprintf('  Expected: around 180° (pointing WEST)\n');

% ========================================================================
% STEP 4: READ CCMP GRID
% ========================================================================

fprintf('\n========================================\n');
fprintf('STEP 4: LOAD CCMP DATA GRID\n');
fprintf('========================================\n');

base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/ccmp_data/';
sample_file = fullfile(base_path, 'Y2024', 'M01', ...
    'CCMP_Wind_Analysis_202401_monthly_mean_V03.1_L4.nc');

if isfile(sample_file)
    lon_ccmp = ncread(sample_file, 'longitude');
    lat_ccmp = ncread(sample_file, 'latitude');
    fprintf('CCMP grid: %d lon points, %d lat points\n', ...
        length(lon_ccmp), length(lat_ccmp));
else
    error('Sample CCMP file not found: %s', sample_file);
end

% Subset to bounding box
idx_lon = find(lon_ccmp >= LONMIN & lon_ccmp <= LONMAX);
idx_lat = find(lat_ccmp >= LATMIN & lat_ccmp <= LATMAX);

lon_sub = lon_ccmp(idx_lon);
lat_sub = lat_ccmp(idx_lat);

fprintf('Subset: %d lon points, %d lat points\n', length(lon_sub), length(lat_sub));

% Create meshgrid
[LON_GRID, LAT_GRID] = meshgrid(lon_sub, lat_sub);
LON_GRID_rad = LON_GRID * rad;
LAT_GRID_rad = LAT_GRID * rad;

% ========================================================================
% STEP 5: CALCULATE PERPENDICULAR DISTANCE FOR EACH GRID POINT
% ========================================================================

fprintf('\n========================================\n');
fprintf('STEP 5: CALCULATE PERPENDICULAR DISTANCES\n');
fprintf('========================================\n');

% For each grid point, we calculate the perpendicular distance to the 
% extended coastline that passes through Cape Juby
%
% Using the formula for point-to-line distance:
% We project the vector from Cape Juby to each grid point onto the 
% OFFSHORE normal direction
%
% distance = dot_product(point_vector, offshore_normal)

% Vector from Cape Juby to each grid point
dlon_to_grid = (LON_GRID_rad - cj_lon_rad) .* cos(LAT_GRID_rad);
dlat_to_grid = LAT_GRID_rad - cj_lat_rad;

% Project onto OFFSHORE normal direction (dot product)
% Positive = offshore (seaward), Negative = onshore
distance_offshore_km = (dlon_to_grid .* offshore_dir_lon + dlat_to_grid .* offshore_dir_lat) * 6371;

fprintf('Distance field statistics:\n');
fprintf('  Min distance: %.1f km (most inland/onshore)\n', min(distance_offshore_km(:)));
fprintf('  Max distance: %.1f km (most offshore)\n', max(distance_offshore_km(:)));
fprintf('  Mean distance: %.1f km\n', mean(distance_offshore_km(:)));

% ========================================================================
% STEP 6: CREATE THREE MASKS (50, 150, 300 km)
% ========================================================================

fprintf('\n========================================\n');
fprintf('STEP 6: CREATE THREE OFFSHORE MASKS\n');
fprintf('========================================\n');

% Mask 1: 0-50 km offshore
mask_50km = uint8(distance_offshore_km >= 0 & distance_offshore_km <= 50);
n_50 = sum(mask_50km(:));

% Mask 2: 0-150 km offshore
mask_150km = uint8(distance_offshore_km >= 0 & distance_offshore_km <= 150);
n_150 = sum(mask_150km(:));

% Mask 3: 0-300 km offshore
mask_300km = uint8(distance_offshore_km >= 0 & distance_offshore_km <= 300);
n_300 = sum(mask_300km(:));

fprintf('Mask 50 km:   %d grid points (%.1f%%)\n', n_50, 100*n_50/numel(mask_50km));
fprintf('Mask 150 km:  %d grid points (%.1f%%)\n', n_150, 100*n_150/numel(mask_150km));
fprintf('Mask 300 km:  %d grid points (%.1f%%)\n', n_300, 100*n_300/numel(mask_300km));

% ========================================================================
% STEP 7: SAVE MASKS
% ========================================================================

fprintf('\n========================================\n');
fprintf('STEP 7: SAVE MASKS\n');
fprintf('========================================\n');

mask_file = fullfile(save_mask, 'ccmp_offshore_masks.mat');
save(mask_file, 'mask_50km', 'mask_150km', 'mask_300km', ...
    'lon_sub', 'lat_sub', 'distance_offshore_km', ...
    'cj_lat', 'cj_lon', 'cb_lat', 'cb_lon', ...
    'offshore_dir_lon', 'offshore_dir_lat', 'coast_dir_lon', 'coast_dir_lat');

fprintf('All masks saved to: %s\n', mask_file);

% ========================================================================
% STEP 8: VISUALIZE
% ========================================================================

fprintf('\n========================================\n');
fprintf('STEP 8: CREATE VISUALIZATIONS\n');
fprintf('========================================\n');

figure('Position', [100 100 1600 900], 'Color', 'w');

% Plot 1: Distance field
subplot(2, 3, 1);
contourf(LON_GRID, LAT_GRID, distance_offshore_km, 50, 'LineColor', 'none');
colorbar;
colormap(gca, 'turbo');
caxis([0 300]);
hold on;

% Contour lines
contour(LON_GRID, LAT_GRID, distance_offshore_km, [0 50 150 300], 'k', 'LineWidth', 1.5);

% Mark capes
plot(cj_lon, cj_lat, 'wo', 'MarkerSize', 12, 'LineWidth', 2);
plot(cb_lon, cb_lat, 'ws', 'MarkerSize', 12, 'LineWidth', 2);

% Offshore direction arrow
scale = 5;
quiver(cj_lon, cj_lat, offshore_dir_lon*scale, offshore_dir_lat*scale, 0, ...
    'Color', 'red', 'LineWidth', 3);

xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Perpendicular Distance from Coastline (km)');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');
grid on;

% Plot 2: Mask 50 km
subplot(2, 3, 2);
imagesc(lon_sub, lat_sub, double(mask_50km));
colormap(gca, 'bone');
hold on;
contour(LON_GRID, LAT_GRID, distance_offshore_km, [0 50], 'r', 'LineWidth', 2);
plot(cj_lon, cj_lat, 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(cb_lon, cb_lat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Mask 50 km (white = included)');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');

% Plot 3: Mask 150 km
subplot(2, 3, 3);
imagesc(lon_sub, lat_sub, double(mask_150km));
colormap(gca, 'bone');
hold on;
contour(LON_GRID, LAT_GRID, distance_offshore_km, [0 150], 'r', 'LineWidth', 2);
plot(cj_lon, cj_lat, 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(cb_lon, cb_lat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Mask 150 km (white = included)');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');

% Plot 4: Mask 300 km
subplot(2, 3, 4);
imagesc(lon_sub, lat_sub, double(mask_300km));
colormap(gca, 'bone');
hold on;
contour(LON_GRID, LAT_GRID, distance_offshore_km, [0 300], 'r', 'LineWidth', 2);
plot(cj_lon, cj_lat, 'go', 'MarkerSize', 10, 'LineWidth', 2);
plot(cb_lon, cb_lat, 'gs', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('Mask 300 km (white = included)');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');

% Plot 5: All three contours together
subplot(2, 3, 5);
contourf(LON_GRID, LAT_GRID, distance_offshore_km, 50, 'LineColor', 'none');
colormap(gca, 'turbo');
caxis([0 300]);
hold on;
contour(LON_GRID, LAT_GRID, distance_offshore_km, [50 50], 'b', 'LineWidth', 2.5, ...
    'DisplayName', '50 km');
contour(LON_GRID, LAT_GRID, distance_offshore_km, [150 150], 'g', 'LineWidth', 2.5, ...
    'DisplayName', '150 km');
contour(LON_GRID, LAT_GRID, distance_offshore_km, [300 300], 'r', 'LineWidth', 2.5, ...
    'DisplayName', '300 km');
contour(LON_GRID, LAT_GRID, distance_offshore_km, [0 0], 'k', 'LineWidth', 2, ...
    'DisplayName', 'Coastline');
plot(cj_lon, cj_lat, 'wo', 'MarkerSize', 12, 'LineWidth', 2);
plot(cb_lon, cb_lat, 'ws', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Longitude (°)');
ylabel('Latitude (°)');
title('All Three Contours');
legend('Location', 'best');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');
grid on;

% Plot 6: Distance profile
subplot(2, 3, 6);
lat_profile_idx = find(lat_sub >= 26 & lat_sub <= 28);
lat_profile_idx = lat_profile_idx(ceil(length(lat_profile_idx)/2));
distance_profile = distance_offshore_km(lat_profile_idx, :);
plot(lon_sub, distance_profile, 'b-', 'LineWidth', 2);
hold on;
yline(50, 'b--', 'LineWidth', 1.5, 'DisplayName', '50 km');
yline(150, 'g--', 'LineWidth', 1.5, 'DisplayName', '150 km');
yline(300, 'r--', 'LineWidth', 1.5, 'DisplayName', '300 km');
yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Coastline');
xlabel('Longitude (°)');
ylabel('Distance offshore (km)');
title(sprintf('Distance Profile at Lat = %.2f°N', lat_sub(lat_profile_idx)));
legend('Location', 'best');
xlim([LONMIN LONMAX]);
grid on;

sgtitle('Three Offshore Distance Masks (50, 150, 300 km)', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(save_path_plot, 'three_offshore_masks_visualization.png'));
fprintf('Visualization saved\n');

close(gcf);

fprintf('\n✓ Complete! Three masks created and saved.\n');
fprintf('\nTo use these masks in your CCMP processing:\n');
fprintf('  1. Load: load(''%s'', ''mask_50km'', ''mask_150km'', ''mask_300km'')\n', mask_file);
fprintf('  2. Apply: ccmp_data_masked = ccmp_data .* mask;\n');