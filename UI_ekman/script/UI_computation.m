clear; clc; close all;

% ========================================================================
% CALCULATE EKMAN UPWELLING INDEX (UI_EKMAN) FROM ALL MASKED DATA
% ========================================================================
% Formula: UI_Ekman = tau_along / (rho * f)
% Data dimensions from file: [time=4, latitude=128, longitude=192]
% ========================================================================

% Define paths
base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products/masked';
output_base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman_results_new';

% Create output directory if it doesn't exist
if ~isfolder(output_base_path)
    mkdir(output_base_path);
end

% Constants
rho = 1025;                    % Seawater density [kg/m³]
Omega = 7.2921e-5;            % Earth's angular velocity [rad/s]

% Build up the subset
days = 1:31;
months = 1:12;
years_numbers = 2017:2020;    % All years
mask_distances = {'50km'};    % All masks

% ========================================================================
% PRE-CALCULATE CORIOLIS PARAMETER
% ========================================================================
fprintf('Pre-calculating Coriolis parameter...\n');

% Load sample file to get grid
sample_file = fullfile(base_path, 'masked_50km', 'Y2005', 'M01', ...
    'CCMP_products_20050101_V03.1_L4_masked_50km.nc');

% Read coordinate variables
lon = ncread(sample_file, 'longitude');
lat = ncread(sample_file, 'latitude');
n_lon = length(lon);
n_lat = length(lat);

fprintf('Grid dimensions: lon=%d, lat=%d\n', n_lon, n_lat);

% Create latitude grid for broadcasting [lat, lon]
lat_grid = repmat(lat, 1, n_lon);

% Calculate Coriolis parameter [lat, lon]
f_param = 2 * Omega * sin(deg2rad(lat_grid));
f_param(abs(f_param) < 1e-6) = 1e-6; % Handle near-equatorial values

fprintf('Coriolis f range: [%.2e to %.2e] rad/s\n\n', min(f_param(:)), max(f_param(:)));

% ========================================================================
% LOOP OVER FILES AND MASKS
% ========================================================================

total_files = length(years_numbers) * length(months) * length(days) * length(mask_distances);
file_count = 0;
successful_files = 0;
failed_files = 0;

tic;

for y = 1:length(years_numbers)
    
    year = years_numbers(y);
    Ystr = sprintf('Y%04d', year);
    
    % Determine if leap year
    if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
        days_in_year = 366;
    else
        days_in_year = 365;
    end
    
    for m = 1:length(months)
        
        month = months(m);
        Mstr = sprintf('M%02d', month);
        
        % Get number of days in this month
        if month == 2
            if days_in_year == 366
                n_days = 29;
            else
                n_days = 28;
            end
        elseif ismember(month, [4, 6, 9, 11])
            n_days = 30;
        else
            n_days = 31;
        end
        
        for d = 1:n_days
            
            day = d;
            
            for mask_idx = 1:length(mask_distances)
                
                mask_name = mask_distances{mask_idx};
                file_count = file_count + 1;
                
                % Build filename safely
                year_str = num2str(year);
                month_str = sprintf('%02d', month);
                day_str = sprintf('%02d', day);
                
                input_filename = ['CCMP_products_' year_str month_str day_str ...
                    '_V03.1_L4_masked_' mask_name '.nc'];
                
                mask_dir = ['masked_' mask_name];
                input_filepath = fullfile(base_path, mask_dir, Ystr, Mstr, input_filename);
                
                % Check if file exists
                if ~isfile(input_filepath)
                    if mod(file_count, 1000) == 0
                        fprintf('File %d/%d not found\n', file_count, total_files);
                    end
                    continue;
                end
                
                if mod(file_count, 500) == 0 || file_count == 1
                    elapsed = toc;
                    fprintf('File %d/%d: %s%s%s_%s (%.1f sec, %d successful, %d failed)\n', ...
                        file_count, total_files, year_str, month_str, day_str, mask_name, elapsed, successful_files, failed_files);
                end
                
                % ============================================================
                % READ AND VALIDATE FILE
                % ============================================================
                try
                    % Read tau_along variable
                    tau_along = ncread(input_filepath, 'tau_along');
                    
                    % Get actual dimensions
                    actual_shape = size(tau_along);
                    
                    % Expected values for each dimension
                    expected_dims = [4, 128, 192]; % time, lat, lon
                    
                    % Find which position each dimension occupies
                    time_idx = find(actual_shape == 4, 1);
                    lat_idx = find(actual_shape == 128, 1);
                    lon_idx = find(actual_shape == 192, 1);
                    
                    % Validate we found all dimensions
                    if isempty(time_idx) || isempty(lat_idx) || isempty(lon_idx)
                        error('Dimension mismatch: expected [4, 128, 192] (time, lat, lon), got [%s]', ...
                            sprintf('%d ', actual_shape));
                    end
                    
                    % Permute to standard order [time, lat, lon] = [dim1, dim2, dim3]
                    if ~isequal([time_idx, lat_idx, lon_idx], [1, 2, 3])
                        tau_along = permute(tau_along, [time_idx, lat_idx, lon_idx]);
                    end
                    
                    % Verify final shape
                    [n_time, n_lat_file, n_lon_file] = size(tau_along);
                    
                    if n_lat_file ~= n_lat || n_lon_file ~= n_lon || n_time ~= 4
                        error('After permute: expected [4, %d, %d], got [%d, %d, %d]', ...
                            n_lat, n_lon, n_time, n_lat_file, n_lon_file);
                    end
                    
                catch ME
                    fprintf('  ERROR reading %s: %s\n', input_filename, ME.message);
                    failed_files = failed_files + 1;
                    continue;
                end
                
                % ============================================================
                % CALCULATE UI_EKMAN FOR EACH TIME STEP
                % ============================================================
                
                try
                    % Initialize output array [time, lat, lon]
                    UI_ekman = zeros(n_time, n_lat, n_lon, 'single');
                    
                    % Calculate UI_Ekman for each time step
                    for t = 1:n_time
                        % Extract tau_along for current time step [lat, lon]
                        tau_along_t = squeeze(tau_along(t, :, :));
                        
                        % Calculate UI_Ekman = tau_along / (rho * f)
                        % f_param is [lat, lon], tau_along_t is [lat, lon]
                        UI_ekman(t, :, :) = tau_along_t ./ (rho * f_param);
                    end
                    
                catch ME
                    fprintf('  ERROR calculating UI_ekman for %s: %s\n', input_filename, ME.message);
                    failed_files = failed_files + 1;
                    continue;
                end
                
                % ============================================================
                % SAVE TO OUTPUT NETCDF FILE
                % ============================================================
                
                try
                    % Create output directory
                    output_dir = fullfile(output_base_path, ['UI_ekman_' mask_name], Ystr, Mstr);
                    if ~isfolder(output_dir)
                        mkdir(output_dir);
                    end
                    
                    % Output filename
                    output_filename = ['UI_ekman_' year_str month_str day_str '_' mask_name '.nc'];
                    output_filepath = fullfile(output_dir, output_filename);
                    
                    % Delete file if it already exists
                    if isfile(output_filepath)
                        delete(output_filepath);
                    end
                    
                    % Use netcdf library directly
                    mode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
                    ncid = netcdf.create(output_filepath, mode);
                    
                    % IMPORTANT: Define dimensions in order [time, latitude, longitude]
                    % MATLAB will reverse these due to Fortran column-major order
                    time_dimid = netcdf.defDim(ncid, 'time', n_time);
                    lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
                    lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);
                    
                    % Define coordinate variables
                    time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);
                    lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
                    lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
                    
                    % Define UI_Ekman variable
                    % To get (time, lat, lon) in the file, define as (lon, lat, time)
                    % because MATLAB reverses dimension order when writing
                    ui_varid = netcdf.defVar(ncid, 'UI_ekman', 'NC_FLOAT', ...
                        [lon_dimid, lat_dimid, time_dimid]);
                    
                    % End define mode
                    netcdf.endDef(ncid);
                    
                    % Write coordinate variables
                    netcdf.putVar(ncid, time_varid, int32(1:n_time));
                    netcdf.putVar(ncid, lat_varid, double(lat));
                    netcdf.putVar(ncid, lon_varid, double(lon));
                    
                    % Permute data from [time, lat, lon] to [lon, lat, time]
                    % to compensate for MATLAB's Fortran order
                    UI_ekman_permuted = permute(UI_ekman, [3, 2, 1]);
                    
                    % Write UI_Ekman data
                    netcdf.putVar(ncid, ui_varid, single(UI_ekman_permuted));
                    
                    % Add attributes
                    netcdf.reDef(ncid);
                    
                    netcdf.putAtt(ncid, time_varid, 'long_name', 'Time');
                    netcdf.putAtt(ncid, time_varid, 'units', 'observations per day');
                    
                    netcdf.putAtt(ncid, lat_varid, 'standard_name', 'latitude');
                    netcdf.putAtt(ncid, lat_varid, 'units', 'degrees_north');
                    
                    netcdf.putAtt(ncid, lon_varid, 'standard_name', 'longitude');
                    netcdf.putAtt(ncid, lon_varid, 'units', 'degrees_east');
                    
                    netcdf.putAtt(ncid, ui_varid, 'long_name', 'Ekman Upwelling Index');
                    netcdf.putAtt(ncid, ui_varid, 'units', 'm^2 s^-1');
                    netcdf.putAtt(ncid, ui_varid, 'formula', 'tau_along / (rho * f)');
                    netcdf.putAtt(ncid, ui_varid, 'mask_applied', mask_name);
                    netcdf.putAtt(ncid, ui_varid, 'positive', 'up');
                    
                    % Global attributes
                    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', ...
                        sprintf('UI_Ekman - %s mask', mask_name));
                    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'date', ...
                        sprintf('%04d-%02d-%02d', year, month, day));
                    netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'source', ...
                        'CCMP v03.1 masked data');
                    
                    % Close file
                    netcdf.close(ncid);
                    successful_files = successful_files + 1;
                    
                catch ME
                    fprintf('  ERROR writing %s: %s\n', output_filename, ME.message);
                    failed_files = failed_files + 1;
                    continue;
                end
                
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('✓ UI_EKMAN CALCULATION COMPLETE!\n');
fprintf('Total files attempted: %d\n', file_count);
fprintf('Successful files: %d\n', successful_files);
fprintf('Failed files: %d\n', failed_files);
fprintf('Results saved to: %s\n', output_base_path);
fprintf('Total time: %.1f seconds\n', toc);
fprintf('========================================\n');