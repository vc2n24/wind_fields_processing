clear; clc; close all;

% ========================================================================
% CREATE THREE OFFSHORE DISTANCE MASKS (50, 150, 300 km)
% ========================================================================

save_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/upwelling_mask_script';

if ~isfolder(save_path)
    mkdir(save_path);
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
cb_lat= 20.91;
cb_lon= -17.06;


% Cape Juby (southernmost)
cj_lat= 27.89;
cj_lon= -13.03;



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
offshore_dir_lon = -coast_dir_lat;
offshore_dir_lat = coast_dir_lon;

fprintf('Offshore direction (perpendicular to coast, +90° rotation):\n');
fprintf('  offshore_dir_lon = %.6f (should be NEGATIVE = westward)\n', offshore_dir_lon);
fprintf('  offshore_dir_lat = %.6f (small value)\n', offshore_dir_lat);

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
lat_sub= lat_ccmp(idx_lat);

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

% Vector from Cape Juby to each grid point
dlon_to_grid = (LON_GRID_rad - cj_lon_rad) .* cos(LAT_GRID_rad);
dlat_to_grid = LAT_GRID_rad - cj_lat_rad;

% Project onto OFFSHORE normal direction (dot product)
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

mask_file = fullfile(save_path, 'ccmp_offshore_masks.mat');
save(mask_file, 'mask_50km', 'mask_150km', 'mask_300km', ...
    'lon_sub', 'lat_sub', 'distance_offshore_km', ...
    'cj_lat', 'cj_lon', 'cb_lat', 'cb_lon', ...
    'offshore_dir_lon', 'offshore_dir_lat', 'coast_dir_lon', 'coast_dir_lat');

fprintf('All masks saved to: %s\n', mask_file);

fprintf('\n✓ Mask creation complete!\n');