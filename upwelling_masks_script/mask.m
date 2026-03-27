clear; clc; close all;

% ========================================================================
% CREATE EXTENDED PARALLEL MASK FOR ENTIRE CCMP DOMAIN (FIXED VERSION)
% ========================================================================

save_path_plot = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/plot/plot_mask';
save_mask= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/upwelling_masks_script'
if ~isfolder(save_path_plot)
    mkdir(save_path_plot);
end

% Full bounding box for all data
LONMIN = -40;
LONMAX = 8;
LATMIN = 8;
LATMAX = 40;

% ========================================================================
% STEP 1: DEFINE TWO CAPE POINTS (Cape Bojador and Cape Juby)
% ========================================================================

% Cape Bojador (northern point)
cb_lat = 26.11;
cb_lon = -14.33;

% Cape Juby (southern point)
cj_lat = 27.96;
cj_lon = -13.20;

fprintf('Cape Bojador: (%.2f°N, %.2f°W)\n', cb_lat, abs(cb_lon));
fprintf('Cape Juby:    (%.2f°N, %.2f°W)\n', cj_lat, abs(cj_lon));

% Convert to radians
rad = pi / 180;

cj_lat_rad = cj_lat * rad;
cj_lon_rad = cj_lon * rad;
cb_lat_rad = cb_lat * rad;
cb_lon_rad = cb_lon * rad;

% ========================================================================
% STEP 2: CALCULATE THE COASTLINE DIRECTION
% ========================================================================

% Vector from Cape Bojador to Cape Juby
dlat = cj_lat_rad - cb_lat_rad;
dlon = (cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad + cb_lat_rad) / 2);

% This vector points ALONG the coast (roughly north-south)
% Normalize it
len_coast = sqrt(dlat^2 + dlon^2);
coast_lat = dlat / len_coast;  % unit vector component (latitude direction)
coast_lon = dlon / len_coast;  % unit vector component (longitude direction)

fprintf('\nCoastline direction (from Bojador to Juby):\n');
fprintf('  Angle: %.2f degrees from East\n', atan2(dlat, dlon) * 180/pi);
fprintf('  Unit vector: [dlon=%.4f, dlat=%.4f]\n', coast_lon, coast_lat);

% OFFSHORE NORMAL: perpendicular to coastline, pointing SEAWARD
% Rotate the coastline vector 90 degrees CLOCKWISE (to point seaward/eastward)
% If coast_vector = [coast_lon, coast_lat]
% Then perpendicular (rotated -90°) = [coast_lat, -coast_lon]
offshore_lon = coast_lat;   % This points roughly eastward (seaward)
offshore_lat = -coast_lon;  % This points roughly southward

fprintf('\nOffshore (seaward) normal:\n');
fprintf('  Unit vector: [dlon=%.4f, dlat=%.4f]\n', offshore_lon, offshore_lat);
fprintf('  This should point roughly EAST (positive lon)\n');

% Convert 50 km to radians
R = 6371;  % Earth radius in km
off_rad = 50 / R;

fprintf('\n50 km converted to radians: %.6f rad\n', off_rad);

% ========================================================================
% STEP 3: READ CCMP DATA TO GET GRID STRUCTURE
% ========================================================================

base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/ccmp_data/';
sample_file = fullfile(base_path, 'Y2024', 'M01', ...
    'CCMP_Wind_Analysis_202401_monthly_mean_V03.1_L4.nc');

if isfile(sample_file)
    lon_ccmp = ncread(sample_file, 'longitude');
    lat_ccmp = ncread(sample_file, 'latitude');
    fprintf('\nCCMP grid loaded: %d lon points, %d lat points\n', ...
        length(lon_ccmp), length(lat_ccmp));
else
    error('Sample CCMP file not found: %s', sample_file);
end

% Find indices within bounding box
idx_lon = find(lon_ccmp >= LONMIN & lon_ccmp <= LONMAX);
idx_lat = find(lat_ccmp >= LATMIN & lat_ccmp <= LATMAX);

lon_ccmp_sub = lon_ccmp(idx_lon);
lat_ccmp_sub = lat_ccmp(idx_lat);

fprintf('Subset grid: %d lon points, %d lat points\n', ...
    length(lon_ccmp_sub), length(lat_ccmp_sub));

% Create mesh grid
[LON_GRID, LAT_GRID] = meshgrid(lon_ccmp_sub, lat_ccmp_sub);

% ========================================================================
% STEP 4: COMPUTE SIGNED DISTANCE FROM COASTLINE
% ========================================================================

% Convert grid to radians
LON_GRID_rad = LON_GRID * rad;
LAT_GRID_rad = LAT_GRID * rad;

% For each grid point, compute vector from Cape Juby to that point
dlon_to_point = (LON_GRID_rad - cj_lon_rad) .* cos(LAT_GRID_rad);
dlat_to_point = LAT_GRID_rad - cj_lat_rad;

% Project onto OFFSHORE normal vector (dot product)
% This gives distance perpendicular to the coastline
% Positive = offshore (seaward), Negative = onshore (toward land)
projection = dlon_to_point .* offshore_lon + dlat_to_point .* offshore_lat;

% Distance in km
distance_km = projection * R;

fprintf('\nDistance statistics (km from coastal line):\n');
fprintf('  Min: %.2f km\n', min(distance_km(:)));
fprintf('  Max: %.2f km\n', max(distance_km(:)));
fprintf('  Mean: %.2f km\n', mean(distance_km(:)));

% ========================================================================
% STEP 5: CREATE BINARY MASK
% ========================================================================

% Mask = 1 if between -10 km (inland buffer) and +50 km (offshore boundary)
% This captures the coastal strip plus some inland region for safety
mask_onshore_buffer = 10;  % km inland
mask_offshore_limit = 50;   % km seaward

mask = zeros(length(lat_ccmp_sub), length(lon_ccmp_sub), 'uint8');

% Mark points in the coastal band
mask(distance_km >= -mask_onshore_buffer & distance_km <= mask_offshore_limit) = 1;

n_coastal = sum(mask(:));
fprintf('\nMask statistics:\n');
fprintf('  Grid points in coastal region: %d (out of %d)\n', ...
    n_coastal, numel(mask));
fprintf('  Percentage: %.1f%%\n', 100*n_coastal/numel(mask));

% ========================================================================
% STEP 6: SAVE MASK
% ========================================================================

mask_file = fullfile(save_mask, 'ccmp_mask_extended_cape_region.mat');
save(mask_file, 'mask', 'lon_ccmp_sub', 'lat_ccmp_sub', 'distance_km', ...
    'cj_lat', 'cj_lon', 'cb_lat', 'cb_lon', ...
    'offshore_lon', 'offshore_lat', 'coast_lon', 'coast_lat', ...
    'R', 'off_rad');

fprintf('\nMask saved to: %s\n', mask_file);

% ========================================================================
% STEP 7: VISUALIZE THE MASK AND PARALLEL LINES
% ========================================================================

figure('Position', [100 100 1400 600], 'Color', 'w');

% Plot 1: Distance field with contours
subplot(1, 2, 1);
contourf(LON_GRID, LAT_GRID, distance_km, 30, 'LineColor', 'none');
colorbar;
colormap(gca, 'jet');
hold on;

% Plot the two parallel lines
% Coastline (distance = 0)
contour(LON_GRID, LAT_GRID, distance_km, [0 0], 'k', 'LineWidth', 2.5, ...
    'DisplayName', 'Coastline');
% 50 km offshore line
contour(LON_GRID, LAT_GRID, distance_km, [50 50], 'r', 'LineWidth', 2.5, ...
    'DisplayName', '50 km offshore');
% Inland buffer
contour(LON_GRID, LAT_GRID, distance_km, [-10 -10], 'b', 'LineStyle', '--', ...
    'LineWidth', 1.5, 'DisplayName', '10 km inland');

% Mark the two capes
plot(cj_lon, cj_lat, 'o', 'Color', [0 0.8 0], 'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', 'Cape Juby');
plot(cb_lon, cb_lat, 's', 'Color', [0.8 0 0], 'MarkerSize', 12, 'LineWidth', 2, ...
    'DisplayName', 'Cape Bojador');

% Add arrow showing offshore direction (at Cape Juby)
arrow_scale = 3;
quiver(cj_lon, cj_lat, offshore_lon*arrow_scale, offshore_lat*arrow_scale, 0, ...
    'Color', 'magenta', 'LineWidth', 2.5, 'MaxHeadSize', 0.5);
text(cj_lon + offshore_lon*4, cj_lat + offshore_lat*4, 'Offshore →', ...
    'Color', 'magenta', 'FontSize', 10, 'FontWeight', 'bold');

xlabel('Longitude', 'FontSize', 12);
ylabel('Latitude', 'FontSize', 12);
title('Distance Field from Coastline (km)');
legend('Location', 'best', 'FontSize', 10);
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');
grid on;

% Plot 2: Binary mask
subplot(1, 2, 2);
mask_display = mask;
imagesc(lon_ccmp_sub, lat_ccmp_sub, mask_display);
colormap(gca, 'gray');
hold on;

% Overlay the boundary lines
contour(LON_GRID, LAT_GRID, distance_km, [0 0], 'r', 'LineWidth', 2, ...
    'DisplayName', 'Coastline');
contour(LON_GRID, LAT_GRID, distance_km, [50 50], 'b', 'LineWidth', 2, ...
    'DisplayName', '50 km offshore');

% Mark capes
plot(cj_lon, cj_lat, 'o', 'Color', [0 0.8 0], 'MarkerSize', 12, 'LineWidth', 2);
plot(cb_lon, cb_lat, 's', 'Color', [0.8 0 0], 'MarkerSize', 12, 'LineWidth', 2);

xlabel('Longitude', 'FontSize', 12);
ylabel('Latitude', 'FontSize', 12);
title('Binary Mask (white = coastal region, black = offshore)');
xlim([LONMIN LONMAX]);
ylim([LATMIN LATMAX]);
set(gca, 'YDir', 'normal');
legend('Location', 'best', 'FontSize', 10);

sgtitle('Extended Coastal Mask - Cape Bojador to Cape Juby', 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, fullfile(save_path_plot, 'mask_extended_visualization_FIXED.png'));
fprintf('Mask visualization saved\n');

close(gcf);

% ========================================================================
% STEP 8: CREATE PROFILE PLOT
% ========================================================================

figure('Position', [100 100 1200 500], 'Color', 'w');

% Profile along a latitude line near the capes
lat_profile = mean([cb_lat cj_lat]);
[~, idx_lat_profile] = min(abs(lat_ccmp_sub - lat_profile));

dist_profile = distance_km(idx_lat_profile, :);
lon_profile = lon_ccmp_sub;

plot(lon_profile, dist_profile, 'b-', 'LineWidth', 2.5);
hold on;

% Reference lines
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Coastline (0 km)');
yline(50, 'r--', 'LineWidth', 1.5, 'DisplayName', '50 km offshore');
yline(-10, 'b--', 'LineWidth', 1.5, 'DisplayName', '10 km inland');

% Shade the mask region
ylim_vals = ylim;
patch([lon_profile(1) lon_profile(end) lon_profile(end) lon_profile(1)], ...
      [-10 -10 50 50], 'green', 'FaceAlpha', 0.1, 'EdgeColor', 'none', ...
      'DisplayName', 'Mask region');

xlabel('Longitude', 'FontSize', 12);
ylabel('Distance from Coastline (km)', 'FontSize', 12);
title(sprintf('Distance Profile at Latitude = %.2f°N', lat_profile));
grid on;
legend('FontSize', 10, 'Location', 'best');
xlim([LONMIN LONMAX]);

saveas(gcf, fullfile(save_path_plot, 'distance_profile_FIXED.png'));
fprintf('Distance profile saved\n');

close(gcf);

fprintf('\n✓ Fixed mask creation complete!\n');