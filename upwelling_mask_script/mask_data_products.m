clear; clc; close all;

% ========================================================================
% APPLY OFFSHORE MASKS TO CCMP DATA PRODUCTS AND SAVE TO NEW NETCDF FILES
% ========================================================================

base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products/';
save_mask = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/upwelling_mask_script';
output_base_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products/masked';

% Create output directory if it doesn't exist
if ~isfolder(output_base_path)
    mkdir(output_base_path);
end

% ========================================================================
% LOAD MASKS
% ========================================================================

fprintf('Loading masks...\n');
mask_file = fullfile(save_mask, 'ccmp_offshore_masks.mat');

if ~isfile(mask_file)
    error('Mask file not found at: %s\n', mask_file);
end

load(mask_file, 'mask_50km', 'mask_150km', 'mask_300km', ...
    'lon_sub', 'lat_sub');

fprintf('Masks loaded successfully\n');

% Store masks in a cell array for easier processing
masks = {
    struct('name', '50km', 'mask', mask_50km);
    struct('name', '150km', 'mask', mask_150km);
    struct('name', '300km', 'mask', mask_300km);
};

% ========================================================================
% PROCESSING PARAMETERS
% ========================================================================

days = 1:31;
months = 1:12;
years_numbers = 2024:2024;

% ========================================================================
% MAIN LOOP OVER FILES
% ========================================================================

total_files = length(years_numbers) * length(months) * length(days);
file_count = 0;

for y = 1:length(years_numbers)
    
    year = years_numbers(y);
    Ystr = sprintf('Y%04d', year);
    Ynumber = sprintf('%04d', year);
    
    for m = 1:length(months)
        
        month = months(m);
        Mstr = sprintf('M%02d', month);
        Mnumber = sprintf('%02d', month);
        
        for d = 1:length(days)
            
            day = days(d);
            Dnumber = sprintf('%02d', day);
            
            file_count = file_count + 1;
            
            % ============================================================
            % READ INPUT FILE
            % ============================================================
            
            file_name = sprintf('CCMP_products_%s%s%s_V03.1_L4.nc', Ynumber, Mnumber, Dnumber);
            file_path = fullfile(base_path, Ystr, Mstr, file_name);
            
            if ~isfile(file_path)
                warning('File %d/%d not found: %s', file_count, total_files, file_path);
                continue;
            end
            
            fprintf('\nProcessing file %d/%d: %s\n', file_count, total_files, file_name);
            
            try
                % Read coordinates
                lon = ncread(file_path, 'longitude');
                lat = ncread(file_path, 'latitude');
                
                % READ TIME - handle missing attributes gracefully
                time_original = ncread(file_path, 'time');
                
                % Try to read time units, if it doesn't exist, use a default
                try
                    time_units = ncreadatt(file_path, 'time', 'units');
                catch
                    fprintf('  Warning: time units attribute not found, using default\n');
                    time_units = 'days since 1970-01-01 00:00:00';
                end
                
                % Read data (already rotated in data_products)
                u_along = ncread(file_path, 'u_along');
                u_cross = ncread(file_path, 'u_cross');
                taux = ncread(file_path, 'taux');
                tauy = ncread(file_path, 'tauy');
                tau_along = ncread(file_path, 'tau_along');
                tau_cross = ncread(file_path, 'tau_cross');
                
                fprintf('  Input data size: lon=%d, lat=%d, time=%d\n', ...
                    length(lon), length(lat), size(u_along, 3));
                fprintf('  Time units: %s\n', time_units);
                fprintf('  Time values (first 5): %s\n', sprintf('%g ', time_original(1:min(5,length(time_original)))));
                
            catch ME
                warning('Error reading file: %s', ME.message);
                fprintf('  Skipping file...\n');
                continue;
            end
            
            % ============================================================
            % FIND GRID INDICES MATCHING MASK COORDINATES
            % ============================================================
            
            % Find indices in CCMP grid that match mask grid
            [~, idx_lon_start] = min(abs(lon - min(lon_sub)));
            [~, idx_lon_end] = min(abs(lon - max(lon_sub)));
            [~, idx_lat_start] = min(abs(lat - min(lat_sub)));
            [~, idx_lat_end] = min(abs(lat - max(lat_sub)));
            
            fprintf('  Grid indices: lon=[%d:%d], lat=[%d:%d]\n', ...
                idx_lon_start, idx_lon_end, idx_lat_start, idx_lat_end);
            
            % Extract matching subset from all data
            u_along_sub = u_along(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            u_cross_sub = u_cross(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            taux_sub = taux(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            tauy_sub = tauy(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            tau_along_sub = tau_along(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            tau_cross_sub = tau_cross(idx_lon_start:idx_lon_end, idx_lat_start:idx_lat_end, :);
            
            % Extract corresponding subset of coordinates
            lon_sub_data = lon(idx_lon_start:idx_lon_end);
            lat_sub_data = lat(idx_lat_start:idx_lat_end);
            
            % Extract corresponding subset of time
            time_sub = time_original(1:size(u_along_sub, 3));
            
            fprintf('  Subset data size: lon=%d, lat=%d, time=%d\n', ...
                length(lon_sub_data), length(lat_sub_data), length(time_sub));
            
            % ============================================================
            % PROCESS EACH MASK
            % ============================================================
            
            for mask_idx = 1:length(masks)
                
                current_mask = masks{mask_idx};
                mask_name = current_mask.name;
                mask_data = current_mask.mask;
                
                fprintf('    Applying %s mask...\n', mask_name);
                
                % Ensure mask dimensions match data
                [mask_ny, mask_nx] = size(mask_data);
                [data_nx, data_ny, n_time] = size(u_along_sub);
                
                fprintf('      Mask size: %d x %d, Data size: %d x %d\n', ...
                    mask_nx, mask_ny, data_nx, data_ny);
                
                if mask_ny ~= data_ny || mask_nx ~= data_nx
                    warning('      Mask dimensions (%d x %d) do not match data (%d x %d). Skipping.', ...
                        mask_nx, mask_ny, data_nx, data_ny);
                    continue;
                end

                % Transpose mask to match data structure [lon, lat]
                mask_transpose = mask_data';

                % Expand mask to 3D
                mask_3d = repmat(mask_transpose, [1, 1, n_time]);

                % Create masked versions by setting non-masked regions to NaN
                u_along_masked = u_along_sub;
                u_cross_masked = u_cross_sub;
                taux_masked = taux_sub;
                tauy_masked = tauy_sub;
                tau_along_masked = tau_along_sub;
                tau_cross_masked = tau_cross_sub;

                % Apply mask (set to NaN where mask is 0)
                u_along_masked(mask_3d == 0) = NaN;
                u_cross_masked(mask_3d == 0) = NaN;
                taux_masked(mask_3d == 0) = NaN;
                tauy_masked(mask_3d == 0) = NaN;
                tau_along_masked(mask_3d == 0) = NaN;
                tau_cross_masked(mask_3d == 0) = NaN;
                
                fprintf('      ✓ Mask applied\n');

                % ============================================================
                % CREATE OUTPUT NETCDF FILE
                % ============================================================
                
                % Create output directory for this mask
                output_dir = fullfile(output_base_path, sprintf('masked_%s', mask_name), Ystr, Mstr);
                if ~isfolder(output_dir)
                    mkdir(output_dir);
                end
                
                % Output filename
                output_filename = sprintf('CCMP_products_%s%s%s_V03.1_L4_masked_%s.nc', ...
                    Ynumber, Mnumber, Dnumber, mask_name);
                output_filepath = fullfile(output_dir, output_filename);
                
                % Delete file if it already exists
                if isfile(output_filepath)
                    delete(output_filepath);
                    fprintf('      Deleted existing file\n');
                end
                
                % Create NetCDF file structure with proper dimensions: (lon, lat, time)
                fprintf('      Creating output file: %s\n', output_filename);
                
                % Keep data as [lon, lat, time] - NO PERMUTATION
                u_along_write = u_along_masked;
                u_cross_write = u_cross_masked;
                taux_write = taux_masked;
                tauy_write = tauy_masked;
                tau_along_write = tau_along_masked;
                tau_cross_write = tau_cross_masked;
                
                % Create dimensions in order: lon, lat, time
                nccreate(output_filepath, 'longitude', 'Dimensions', {'longitude', data_nx}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                nccreate(output_filepath, 'latitude', 'Dimensions', {'latitude', data_ny}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                nccreate(output_filepath, 'time', 'Dimensions', {'time', n_time}, ...
                    'Datatype', 'double', 'FillValue', double(-9999));
                
                % Write coordinate variables
                ncwrite(output_filepath, 'longitude', single(lon_sub_data));
                ncwrite(output_filepath, 'latitude', single(lat_sub_data));
                ncwrite(output_filepath, 'time', double(time_sub));
                
                % Write coordinate attributes
                ncwriteatt(output_filepath, 'longitude', 'standard_name', 'longitude');
                ncwriteatt(output_filepath, 'longitude', 'units', 'degrees_east');
                
                ncwriteatt(output_filepath, 'latitude', 'standard_name', 'latitude');
                ncwriteatt(output_filepath, 'latitude', 'units', 'degrees_north');
                
                ncwriteatt(output_filepath, 'time', 'standard_name', 'time');
                ncwriteatt(output_filepath, 'time', 'units', time_units);
                ncwriteatt(output_filepath, 'time', 'calendar', 'gregorian');
                
                % Create and write data variables with dimensions [lon, lat, time]
                nccreate(output_filepath, 'u_along', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'u_along', single(u_along_write));
                ncwriteatt(output_filepath, 'u_along', 'long_name', 'Masked along-shore wind component');
                ncwriteatt(output_filepath, 'u_along', 'units', 'm s-1');
                ncwriteatt(output_filepath, 'u_along', 'comment', 'Positive from Cape Juby to Cape Blanc');
                ncwriteatt(output_filepath, 'u_along', 'mask_applied', mask_name);
                
                nccreate(output_filepath, 'u_cross', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'u_cross', single(u_cross_write));
                ncwriteatt(output_filepath, 'u_cross', 'long_name', 'Masked cross-shore wind component');
                ncwriteatt(output_filepath, 'u_cross', 'units', 'm s-1');
                ncwriteatt(output_filepath, 'u_cross', 'comment', 'Positive from on-shore to offshore');
                ncwriteatt(output_filepath, 'u_cross', 'mask_applied', mask_name);
                
                nccreate(output_filepath, 'taux', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'taux', single(taux_write));
                ncwriteatt(output_filepath, 'taux', 'long_name', 'Masked eastward wind stress');
                ncwriteatt(output_filepath, 'taux', 'units', 'N m-2');
                ncwriteatt(output_filepath, 'taux', 'mask_applied', mask_name);
                
                nccreate(output_filepath, 'tauy', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'tauy', single(tauy_write));
                ncwriteatt(output_filepath, 'tauy', 'long_name', 'Masked northward wind stress');
                ncwriteatt(output_filepath, 'tauy', 'units', 'N m-2');
                ncwriteatt(output_filepath, 'tauy', 'mask_applied', mask_name);
                
                nccreate(output_filepath, 'tau_along', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'tau_along', single(tau_along_write));
                ncwriteatt(output_filepath, 'tau_along', 'long_name', 'Masked along-shore wind stress');
                ncwriteatt(output_filepath, 'tau_along', 'units', 'N m-2');
                ncwriteatt(output_filepath, 'tau_along', 'comment', 'Positive from Cape Juby toward Cape Blanc');
                ncwriteatt(output_filepath, 'tau_along', 'mask_applied', mask_name);
                
                nccreate(output_filepath, 'tau_cross', ...
                    'Dimensions', {'longitude', 'latitude', 'time'}, ...
                    'Datatype', 'single', 'FillValue', single(-9999));
                ncwrite(output_filepath, 'tau_cross', single(tau_cross_write));
                ncwriteatt(output_filepath, 'tau_cross', 'long_name', 'Masked cross-shore wind stress');
                ncwriteatt(output_filepath, 'tau_cross', 'units', 'N m-2');
                ncwriteatt(output_filepath, 'tau_cross', 'comment', 'Positive from on-shore to offshore');
                ncwriteatt(output_filepath, 'tau_cross', 'mask_applied', mask_name);
                
                % Add global attributes
                ncwriteatt(output_filepath, '/', 'mask_type', mask_name);
                ncwriteatt(output_filepath, '/', 'mask_description', sprintf('Data filtered to %s offshore band', mask_name));
                ncwriteatt(output_filepath, '/', 'creation_date', string(datetime('now')));
                ncwriteatt(output_filepath, '/', 'source', 'CCMP v03.1');
                
                fprintf('      ✓ Saved successfully\n');
                
            end
            
        end
    end
end

fprintf('\n========================================\n');
fprintf('✓ Processing complete!\n');
fprintf('Masked data saved to: %s\n', output_base_path);
fprintf('========================================\n');