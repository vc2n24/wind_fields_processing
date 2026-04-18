% ========================================================================
% CALCULATE UI_EKMAN CLIMATOLOGY, ANOMALIES & STATISTICS (2005-2024)
% ========================================================================
clear all; close all; clc;

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman_results';
output_path = fullfile(ui_results_path, 'analysis');
plot_path = fullfile(output_path, 'plots_update');

% Create output directories
if ~isfolder(output_path)
    mkdir(output_path);
end
if ~isfolder(plot_path)
    mkdir(plot_path);
end

% PARAMETERS

mask_name = '50km';
years = 2005:2024;  % All years
n_years = length(years);


% LOAD ALL UI DATA


fprintf('Loading UI_ekman data for 2005-2024...\n');

all_data = {};
dates_all = {};
file_count = 0;

for y = 1:n_years
    current_year = years(y);
    Ystr = sprintf('Y%04d', current_year);
    
    % Determine if leap year
    if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
        days_in_year = 366;
    else
        days_in_year = 365;
    end
    
    for m = 1:12
        current_month = m;
        Mstr = sprintf('M%02d', current_month);
        
        % Get number of days in this month and remove leap year days in
        % February by not picking them

        if current_month == 2
                n_days = 28;
        elseif ismember(current_month, [4, 6, 9, 11])
            n_days = 30;
        else
            n_days = 31;
        end
        
        for d = 1:n_days
            current_day = d;
            
            % Build filename
            year_str = sprintf('%04d', current_year);
            month_str = sprintf('%02d', current_month);
            day_str = sprintf('%02d', current_day);
            date_str = [year_str month_str day_str];
            
            input_filename = sprintf('UI_ekman_%s%s%s_%s.nc', year_str, month_str, day_str, mask_name);
            
            ui_dir = sprintf('UI_ekman_%s', mask_name);
            input_filepath = fullfile(ui_results_path, ui_dir, Ystr, Mstr, input_filename);
            
            if ~isfile(input_filepath)
                continue;
            end
            
            % Read UI data
            try
                UI_data = ncread(input_filepath, 'UI_ekman');  % [time, lat, lon]
                lon = ncread(input_filepath, 'longitude');
                lat = ncread(input_filepath, 'latitude');
                
                file_count = file_count + 1;
                all_data{file_count} = UI_data;
                dates_all{file_count} = date_str;
                
                if mod(file_count, 500) == 0
                    fprintf('  Loaded file %d: %s\n', file_count, date_str);
                end
                
            catch ME
                fprintf('  Error reading %s: %s\n', input_filename, ME.message);
                continue;
            end
        end
    end
end

fprintf('✓ Loaded %d files\n', file_count);
fprintf('  Grid size: lon=%d, lat=%d\n', length(lon), length(lat));

[n_lon, n_lat, n_time] = size(UI_data);
fprintf('  Time steps per day: %d\n\n', n_time);

% ========================================================================
% CALCULATE MONTHLY CLIMATOLOGY
% ========================================================================

fprintf('Calculating monthly climatology (calendar month average)...\n');

% Get month for each file
month_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    month_all(f) = month_val;
end

% Initialize: store ALL observations for each month (not averages)
monthly_data_pool = cell(12, 1);
for m = 1:12
    monthly_data_pool{m} = [];
end

% Collect all individual observations from all files
fprintf('  Pooling all observations by month...\n');
for f = 1:file_count
    m = month_all(f);
    ui_data = all_data{f};  % [lon, lat, time] - DO NOT average yet
    
    % Reshape to [lon*lat, time] and append to pool
    ui_reshaped = reshape(ui_data, [], size(ui_data, 3));
    monthly_data_pool{m} = [monthly_data_pool{m}, ui_reshaped];
    
    if mod(f, 500) == 0
        fprintf('    Processed file %d/%d\n', f, file_count);
    end
end

% Initialize monthly climatology: [lon, lat, 12]
monthly_climatology = zeros(n_lon, n_lat, 12, 'single');

% Compute true mean for each month from pooled observations
fprintf('  Computing monthly means from pooled observations...\n');
for m = 1:12
    all_obs = monthly_data_pool{m};  % [lon*lat, all_time_steps_in_month]
    
    if ~isempty(all_obs)
        monthly_mean = mean(all_obs, 2, 'omitnan');  % Mean over all time steps
        monthly_climatology(:, :, m) = reshape(monthly_mean, n_lon, n_lat);
        n_obs = size(all_obs, 2);
        fprintf('    Month %2d: %d observations pooled\n', m, n_obs);
    else
        monthly_climatology(:, :, m) = NaN;
        fprintf('    Month %2d: No data\n', m);
    end
end

fprintf('✓ Monthly climatology calculated for 12 months\n');
fprintf('  Monthly climatology size: %s\n\n', mat2str(size(monthly_climatology)));

% Clear memory
clear monthly_data_pool;
% ========================================================================
% CREATE SMOOTH DAILY CLIMATOLOGY (390 DAYS WITH PADDING)
% ========================================================================

fprintf('Creating smooth daily climatology with point-wise interpolation...\n');

% Days per month 
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
half_month = 15;  % Days to pad before Jan and after Dec

% Calculate midpoint of each month (in the 365-day year)
% Day 1 = Jan 1, so midpoint of Jan = day 15.5, midpoint of Feb = day 46.5, etc.
month_midpoints_365 = zeros(12, 1);
cumulative_days = 0;
for m = 1:12
    month_midpoints_365(m) = cumulative_days + (days_per_month(m) / 2);
    cumulative_days = cumulative_days + days_per_month(m);
end

fprintf('  Month midpoints in 365-day year:\n');
for m = 1:12
    fprintf('    Month %2d: day %.1f\n', m, month_midpoints_365(m));
end

% Create extended 390-day calendar with midpoints + padding
% In extended calendar: Dec midpoint at position 15, Jan midpoint at 15+15=30, etc.
extended_midpoints = zeros(14, 1);  % 12 months + Dec repeat + Jan repeat
extended_values_lon_lat = zeros(n_lon, n_lat, 14, 'single');

% Dec of previous year (padding at start) - position 15
extended_midpoints(1) = 15;
extended_values_lon_lat(:, :, 1) = monthly_climatology(:, :, 12);

% All 12 months of the year - positions 15, 30, 45, 60, ...
for m = 1:12
    extended_midpoints(m + 1) = 15 + month_midpoints_365(m);
    extended_values_lon_lat(:, :, m + 1) = monthly_climatology(:, :, m);
end

% Jan of next year (padding at end) - position 390
extended_midpoints(14) = 390;
extended_values_lon_lat(:, :, 14) = monthly_climatology(:, :, 1);

fprintf('\n  Extended calendar midpoints (390-day period):\n');
for i = 1:14
    fprintf('    Point %2d: day %6.1f\n', i, extended_midpoints(i));
end

% ========================================================================
% SMOOTH INTERPOLATION ACROSS DAYS - INTERPOLATE ONLY BETWEEN MIDPOINTS
% ========================================================================

fprintf('\nSmoothing climatology with cubic spline interpolation between month midpoints...\n');

% For each spatial point, interpolate smoothly across days
climatology_smooth = zeros(n_lon, n_lat, 390, 'single');

for i = 1:n_lon
    if mod(i, 50) == 0
        fprintf('  Processing longitude point %d/%d\n', i, n_lon);
    end
    
    for j = 1:n_lat
        % Extract values at the 14 midpoints for this spatial point
        midpoint_values = squeeze(extended_values_lon_lat(i, j, :));  % [14]
        
        % Check if this point has enough valid data
        valid_idx = ~isnan(midpoint_values);
        n_valid = sum(valid_idx);
        
        if n_valid >= 2
            % Interpolate from the 14 midpoints to all 390 days
            days = 1:390;
            valid_midpoints = extended_midpoints(valid_idx);
            valid_values = midpoint_values(valid_idx);
            
            % Cubic spline interpolation between midpoints
            climatology_smooth(i, j, :) = interp1(valid_midpoints, valid_values, days, 'spline');
        else
            % Not enough valid points
            climatology_smooth(i, j, :) = NaN;
        end
    end
end

fprintf('✓ Smoothing complete\n\n');

% ========================================================================
% EXTRACT FINAL CLIMATOLOGY (365 DAYS, REMOVING PADDING)
% ========================================================================

fprintf('Extracting final daily climatology (365 days)...\n');

% Remove the padding: keep days 16-380 (which corresponds to full year)
% Days 1-15 are Dec padding, days 381-390 are Jan padding
climatology = climatology_smooth(:, :, (half_month+1):(half_month+365));

fprintf('✓ Final climatology extracted: %s\n', mat2str(size(climatology)));
fprintf('  Dimensions: [lon=%d, lat=%d, doy=365]\n');
fprintf('  Days 1-15 removed (Dec padding)\n');
fprintf('  Days 16-380 extracted (full year)\n');
fprintf('  Days 381-390 removed (Jan padding)\n\n');

% ========================================================================
% VERIFY CLIMATOLOGY
% ========================================================================

fprintf('Climatology verification:\n');
n_valid_points = 0;
n_nan_points = 0;

for i = 1:n_lon
    for j = 1:n_lat
        clim_ts = climatology(i, j, :);
        if all(isnan(clim_ts))
            n_nan_points = n_nan_points + 1;
        else
            n_valid_points = n_valid_points + 1;
        end
    end
end

fprintf('  Valid spatial points: %d\n', n_valid_points);
fprintf('  All-NaN spatial points: %d\n', n_nan_points);
fprintf('  Total spatial points: %d\n\n', n_lon * n_lat);

% Show range of values
all_clim_data = climatology(~isnan(climatology));
if ~isempty(all_clim_data)
    fprintf('  Overall climatology range:\n');
    fprintf('    Min: %.6f m²/s\n', min(all_clim_data));
    fprintf('    Max: %.6f m²/s\n', max(all_clim_data));
    fprintf('    Mean: %.6f m²/s\n', mean(all_clim_data));
    fprintf('    Std: %.6f m²/s\n\n', std(all_clim_data));
end

fprintf('✓ Verification complete\n\n');

% ========================================================================
% CALCULATE DAY OF YEAR FOR EACH FILE
% ========================================================================

fprintf('Calculating day of year for each file...\n');

day_of_year_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    day_val = str2double(date_str(7:8));
    
    % Calculate day of year (always 365 days, no Feb 29)
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    day_of_year = sum(days_in_months(1:month_val-1)) + day_val;
    day_of_year_all(f) = day_of_year;
end

fprintf('✓ Day of year calculated\n\n');

% ========================================================================
% CALCULATE ANOMALIES FROM DAILY CLIMATOLOGY
% ========================================================================

fprintf('Calculating anomalies from daily climatology...\n');

% Initialize anomalies array: [lon, lat, file_count]
anomalies_data = zeros(n_lon, n_lat, file_count, 'single');

for f = 1:file_count
    doy = day_of_year_all(f);
    
    % all_data{f} is [lon, lat, time]
    % Average over time (dimension 3)
    ui_mean_time = squeeze(mean(all_data{f}, 3, 'omitnan'));  % [lon, lat]
    
    % Get climatology for this day of year
    clim_field = climatology(:, :, doy);  % [lon, lat]
    
    % Calculate anomaly: daily - climatology
    anomalies_data(:, :, f) = ui_mean_time - clim_field;
    
    if mod(f, 200) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

fprintf('✓ Anomalies calculated: %s\n\n', mat2str(size(anomalies_data)));


% ========================================================================
% CALCULATE SPATIAL AND TEMPORAL STATISTICS
% ========================================================================

fprintf('Calculating spatial and temporal statistics...\n');

% Initialize arrays for statistics
daily_ui_mean = zeros(file_count, 1);
daily_clim_mean = zeros(file_count, 1);
daily_anom_mean = zeros(file_count, 1);
dates_datetime = zeros(file_count, 1);

for f = 1:file_count
    doy = day_of_year_all(f);
    
    % Average over space (ignoring NaN): daily UI
    ui_spatial = all_data{f};  % [lon, lat, time]
    ui_spatial_mean = squeeze(mean(ui_spatial, 3, 'omitnan'));  % [lon, lat]
    daily_ui_mean(f) = mean(ui_spatial_mean(~isnan(ui_spatial_mean)));
    
    % Average over space (ignoring NaN): climatology
    clim_spatial = climatology(:, :, doy);  % [lon, lat]
    daily_clim_mean(f) = mean(clim_spatial(~isnan(clim_spatial)));
    
    % Anomaly
    daily_anom_mean(f) = daily_ui_mean(f) - daily_clim_mean(f);
    
    % Date
    date_str = dates_all{f};
    dates_datetime(f) = datenum(date_str, 'yyyymmdd');
end

% Convert to datetime
dates_datetime = datetime(dates_datetime, 'ConvertFrom', 'datenum');

fprintf('✓ Spatial statistics calculated\n\n');

% ========================================================================
% CALCULATE MONTHLY STATISTICS
% ========================================================================

fprintf('Calculating monthly statistics...\n');

monthly_stats = struct();
monthly_stats.month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
                              'July', 'August', 'September', 'October', 'November', 'December'};
monthly_stats.month_numbers = 1:12;
monthly_stats.ui_mean = zeros(12, 1);
monthly_stats.ui_std = zeros(12, 1);
monthly_stats.ui_min = zeros(12, 1);
monthly_stats.ui_max = zeros(12, 1);
monthly_stats.ui_median = zeros(12, 1);
monthly_stats.anom_mean = zeros(12, 1);
monthly_stats.anom_std = zeros(12, 1);
monthly_stats.anom_min = zeros(12, 1);
monthly_stats.anom_max = zeros(12, 1);
monthly_stats.count = zeros(12, 1);

for m = 1:12
    % Find all files in this month
    month_mask = month_all == m;
    month_ui = daily_ui_mean(month_mask);
    month_anom = daily_anom_mean(month_mask);
    
    % Remove NaN values
    month_ui = month_ui(~isnan(month_ui));
    month_anom = month_anom(~isnan(month_anom));
    
    if ~isempty(month_ui)
        monthly_stats.ui_mean(m) = mean(month_ui);
        monthly_stats.ui_std(m) = std(month_ui);
        monthly_stats.ui_min(m) = min(month_ui);
        monthly_stats.ui_max(m) = max(month_ui);
        monthly_stats.ui_median(m) = median(month_ui);
        monthly_stats.anom_mean(m) = mean(month_anom);
        monthly_stats.anom_std(m) = std(month_anom);
        monthly_stats.anom_min(m) = min(month_anom);
        monthly_stats.anom_max(m) = max(month_anom);
        monthly_stats.count(m) = length(month_ui);
    else
        monthly_stats.ui_mean(m) = NaN;
        monthly_stats.ui_std(m) = NaN;
        monthly_stats.ui_min(m) = NaN;
        monthly_stats.ui_max(m) = NaN;
        monthly_stats.ui_median(m) = NaN;
        monthly_stats.anom_mean(m) = NaN;
        monthly_stats.anom_std(m) = NaN;
        monthly_stats.anom_min(m) = NaN;
        monthly_stats.anom_max(m) = NaN;
        monthly_stats.count(m) = 0;
    end
end

fprintf('✓ Monthly statistics calculated\n\n');

% Display monthly statistics
fprintf('Monthly Statistics Summary:\n');
fprintf('%-12s %10s %10s %10s %10s %10s | %10s %10s %10s %10s | %6s\n', ...
    'Month', 'UI Mean', 'UI Std', 'UI Min', 'UI Max', 'UI Median', ...
    'Anom Mean', 'Anom Std', 'Anom Min', 'Anom Max', 'Count');
fprintf(repmat('-', 1, 130));
fprintf('\n');

for m = 1:12
    fprintf('%-12s %10.4f %10.4f %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f %10.4f | %6d\n', ...
        monthly_stats.month_names{m}, ...
        monthly_stats.ui_mean(m), monthly_stats.ui_std(m), ...
        monthly_stats.ui_min(m), monthly_stats.ui_max(m), monthly_stats.ui_median(m), ...
        monthly_stats.anom_mean(m), monthly_stats.anom_std(m), ...
        monthly_stats.anom_min(m), monthly_stats.anom_max(m), monthly_stats.count(m));
end

fprintf('\n✓ Monthly statistics display complete\n\n');

% ========================================================================
% CALCULATE YEARLY STATISTICS
% ========================================================================

fprintf('Calculating yearly statistics...\n');

yearly_stats = struct();
yearly_stats.years = years';
yearly_stats.ui_mean = zeros(n_years, 1);
yearly_stats.ui_std = zeros(n_years, 1);
yearly_stats.ui_min = zeros(n_years, 1);
yearly_stats.ui_max = zeros(n_years, 1);
yearly_stats.anom_mean = zeros(n_years, 1);
yearly_stats.anom_std = zeros(n_years, 1);
yearly_stats.count = zeros(n_years, 1);

for y = 1:n_years
    current_year = years(y);
    year_mask = year(dates_datetime) == current_year;
    year_ui = daily_ui_mean(year_mask);
    year_anom = daily_anom_mean(year_mask);
    
    year_ui = year_ui(~isnan(year_ui));
    year_anom = year_anom(~isnan(year_anom));
    
    if ~isempty(year_ui)
        yearly_stats.ui_mean(y) = mean(year_ui);
        yearly_stats.ui_std(y) = std(year_ui);
        yearly_stats.ui_min(y) = min(year_ui);
        yearly_stats.ui_max(y) = max(year_ui);
        yearly_stats.anom_mean(y) = mean(year_anom);
        yearly_stats.anom_std(y) = std(year_anom);
        yearly_stats.count(y) = length(year_ui);
    else
        yearly_stats.ui_mean(y) = NaN;
        yearly_stats.ui_std(y) = NaN;
        yearly_stats.ui_min(y) = NaN;
        yearly_stats.ui_max(y) = NaN;
        yearly_stats.anom_mean(y) = NaN;
        yearly_stats.anom_std(y) = NaN;
        yearly_stats.count(y) = 0;
    end
end

fprintf('✓ Yearly statistics calculated\n\n');

% ========================================================================
% CALCULATE OVERALL STATISTICS
% ========================================================================

fprintf('Calculating overall statistics...\n');

overall_stats = struct();
overall_ui_valid = daily_ui_mean(~isnan(daily_ui_mean));
overall_anom_valid = daily_anom_mean(~isnan(daily_anom_mean));

overall_stats.ui_mean = mean(overall_ui_valid);
overall_stats.ui_std = std(overall_ui_valid);
overall_stats.ui_min = min(overall_ui_valid);
overall_stats.ui_max = max(overall_ui_valid);
overall_stats.ui_median = median(overall_ui_valid);
overall_stats.anom_mean = mean(overall_anom_valid);
overall_stats.anom_std = std(overall_anom_valid);
overall_stats.anom_min = min(overall_anom_valid);
overall_stats.anom_max = max(overall_anom_valid);
overall_stats.n_observations = length(overall_ui_valid);
overall_stats.period = sprintf('%d-%d', years(1), years(end));

fprintf('✓ Overall statistics calculated\n');
fprintf('  Period: %s\n', overall_stats.period);
fprintf('  Total observations: %d\n', overall_stats.n_observations);
fprintf('  UI mean: %.4f m²/s\n', overall_stats.ui_mean);
fprintf('  UI std: %.4f m²/s\n\n', overall_stats.ui_std);

% ========================================================================
% CREATE STATISTICS DATASET
% ========================================================================

fprintf('Creating statistics dataset...\n');

statistics_dataset = struct();
statistics_dataset.name = 'UI_Ekman_Statistics';
statistics_dataset.description = 'Comprehensive statistics of Ekman Upwelling Index';
statistics_dataset.units = 'm^2 s^-1';

% Overall statistics with descriptions
statistics_dataset.overall = struct();
statistics_dataset.overall.ui_mean = struct(...
    'value', overall_stats.ui_mean, ...
    'description', 'Overall mean of UI_Ekman across all days and years', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.ui_std = struct(...
    'value', overall_stats.ui_std, ...
    'description', 'Overall standard deviation of UI_Ekman', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.ui_min = struct(...
    'value', overall_stats.ui_min, ...
    'description', 'Minimum UI_Ekman value (most downwelling-favorable)', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.ui_max = struct(...
    'value', overall_stats.ui_max, ...
    'description', 'Maximum UI_Ekman value (most upwelling-favorable)', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.ui_median = struct(...
    'value', overall_stats.ui_median, ...
    'description', 'Median UI_Ekman value', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.anom_mean = struct(...
    'value', overall_stats.anom_mean, ...
    'description', 'Mean anomaly (should be ~0 by definition)', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.anom_std = struct(...
    'value', overall_stats.anom_std, ...
    'description', 'Standard deviation of anomalies', ...
    'units', 'm^2 s^-1');
statistics_dataset.overall.period = overall_stats.period;
statistics_dataset.overall.n_observations = overall_stats.n_observations;

% Monthly statistics with descriptions
statistics_dataset.monthly = struct();
statistics_dataset.monthly.months = monthly_stats.month_names';
statistics_dataset.monthly.ui_mean = struct(...
    'data', monthly_stats.ui_mean, ...
    'description', 'Mean UI_Ekman for each calendar month', ...
    'units', 'm^2 s^-1', ...
    'interpretation', 'Positive = upwelling-favorable; Negative = downwelling-favorable');
statistics_dataset.monthly.ui_std = struct(...
    'data', monthly_stats.ui_std, ...
    'description', 'Standard deviation of UI_Ekman within each month', ...
    'units', 'm^2 s^-1', ...
    'interpretation', 'Variability of upwelling index during each month');
statistics_dataset.monthly.ui_min = struct(...
    'data', monthly_stats.ui_min, ...
    'description', 'Minimum UI_Ekman value in each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.ui_max = struct(...
    'data', monthly_stats.ui_max, ...
    'description', 'Maximum UI_Ekman value in each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.ui_median = struct(...
    'data', monthly_stats.ui_median, ...
    'description', 'Median UI_Ekman value in each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.anom_mean = struct(...
    'data', monthly_stats.anom_mean, ...
    'description', 'Mean anomaly for each calendar month', ...
    'units', 'm^2 s^-1', ...
    'interpretation', 'Deviations from the long-term daily climatology');
statistics_dataset.monthly.anom_std = struct(...
    'data', monthly_stats.anom_std, ...
    'description', 'Standard deviation of anomalies within each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.anom_min = struct(...
    'data', monthly_stats.anom_min, ...
    'description', 'Minimum anomaly in each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.anom_max = struct(...
    'data', monthly_stats.anom_max, ...
    'description', 'Maximum anomaly in each month', ...
    'units', 'm^2 s^-1');
statistics_dataset.monthly.count = struct(...
    'data', monthly_stats.count, ...
    'description', 'Number of observations in each month', ...
    'units', 'count');

% Yearly statistics with descriptions
statistics_dataset.yearly = struct();
statistics_dataset.yearly.years = yearly_stats.years;
statistics_dataset.yearly.ui_mean = struct(...
    'data', yearly_stats.ui_mean, ...
    'description', 'Mean UI_Ekman for each year', ...
    'units', 'm^2 s^-1', ...
    'interpretation', 'Interannual variability in upwelling intensity');
statistics_dataset.yearly.ui_std = struct(...
    'data', yearly_stats.ui_std, ...
    'description', 'Standard deviation of UI_Ekman within each year', ...
    'units', 'm^2 s^-1');
statistics_dataset.yearly.ui_min = struct(...
    'data', yearly_stats.ui_min, ...
    'description', 'Minimum UI_Ekman in each year', ...
    'units', 'm^2 s^-1');
statistics_dataset.yearly.ui_max = struct(...
    'data', yearly_stats.ui_max, ...
    'description', 'Maximum UI_Ekman in each year', ...
    'units', 'm^2 s^-1');
statistics_dataset.yearly.anom_mean = struct(...
    'data', yearly_stats.anom_mean, ...
    'description', 'Mean anomaly for each year', ...
    'units', 'm^2 s^-1');
statistics_dataset.yearly.anom_std = struct(...
    'data', yearly_stats.anom_std, ...
    'description', 'Standard deviation of anomalies within each year', ...
    'units', 'm^2 s^-1');
statistics_dataset.yearly.count = struct(...
    'data', yearly_stats.count, ...
    'description', 'Number of observations in each year', ...
    'units', 'count');

% Add metadata
statistics_dataset.calculation_method = struct(...
    'ui_mean', 'Spatially averaged (lon, lat) mean for each day, then aggregated by month/year', ...
    'ui_std', 'Standard deviation of spatially-averaged daily values', ...
    'ui_min', 'Minimum of spatially-averaged daily values', ...
    'ui_max', 'Maximum of spatially-averaged daily values', ...
    'anom_mean', 'Mean of daily anomalies (daily - climatology)', ...
    'anom_std', 'Standard deviation of daily anomalies');

statistics_dataset.mask_name = mask_name;
statistics_dataset.period = sprintf('%d-%d', years(1), years(end));
statistics_dataset.n_years = n_years;
statistics_dataset.processing_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

fprintf('✓ Statistics dataset created\n\n');

% ========================================================================
% SAVE STATISTICS DATASET TO .MAT
% ========================================================================

fprintf('Saving statistics dataset...\n');

stats_filepath = fullfile(output_path, sprintf('UI_ekman_statistics_%s.mat', mask_name));
save(stats_filepath, 'statistics_dataset', '-v7.3');

fprintf('✓ Saved statistics: %s\n\n', stats_filepath);

% ========================================================================
% SAVE ALL DATASETS TO NETCDF WITH FULL DOCUMENTATION
% ========================================================================

fprintf('Saving all datasets to NetCDF files...\n');

% ========================================================================
% SAVE CLIMATOLOGY TO NETCDF - [lon, lat, time=365]
% ========================================================================

climatology_filepath = fullfile(output_path, sprintf('UI_ekman_daily_climatology_%s.nc', mask_name));
if isfile(climatology_filepath)
    delete(climatology_filepath);
end

mode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(climatology_filepath, mode);

% Define dimensions: [lon, lat, time]
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
time_dimid = netcdf.defDim(ncid, 'time', 365);

% Define coordinate variables
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);

% Define climatology data variable: [lon, lat, time]
climatology_varid = netcdf.defVar(ncid, 'UI_ekman_climatology', 'NC_FLOAT', ...
    [lon_dimid, lat_dimid, time_dimid]);

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, lat_varid, double(lat));
netcdf.putVar(ncid, time_varid, int32(1:365));
netcdf.putVar(ncid, climatology_varid, single(climatology));

% Add attributes
netcdf.reDef(ncid);

netcdf.putAtt(ncid, lon_varid, 'standard_name', 'longitude');
netcdf.putAtt(ncid, lon_varid, 'long_name', 'Longitude');
netcdf.putAtt(ncid, lon_varid, 'units', 'degrees_east');

netcdf.putAtt(ncid, lat_varid, 'standard_name', 'latitude');
netcdf.putAtt(ncid, lat_varid, 'long_name', 'Latitude');
netcdf.putAtt(ncid, lat_varid, 'units', 'degrees_north');

netcdf.putAtt(ncid, time_varid, 'long_name', 'Day of Year');
netcdf.putAtt(ncid, time_varid, 'units', 'days (1-365)');
netcdf.putAtt(ncid, time_varid, 'calendar', 'standard (365 days, no leap day)');

netcdf.putAtt(ncid, climatology_varid, 'long_name', 'Ekman Upwelling Index - Daily Climatology');
netcdf.putAtt(ncid, climatology_varid, 'units', 'm^2 s^-1');
netcdf.putAtt(ncid, climatology_varid, 'description', ...
    'Smoothly interpolated daily climatology of UI_Ekman');

% ========================================================================
% GLOBAL ATTRIBUTES FOR CLIMATOLOGY
% ========================================================================

netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', ...
    'Ekman Upwelling Index - Daily Climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'summary', ...
    'Daily climatology of Ekman Upwelling Index for the Canary Current region');

% CALCULATION METHOD
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'calculation_method', ...
    'MONTHLY_CLIMATOLOGY -> EXTENDED_390_DAYS -> SPLINE_INTERPOLATION -> EXTRACT_365_DAYS');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'calculation_steps', ...
    ['STEP 1: Calculate mean UI_Ekman for each calendar month (average all Januaries, all Februaries, etc). ' ...
     'STEP 2: Extend to 390 days by adding 15-day padding (December values before Jan, January values after Dec). ' ...
     'STEP 3: Apply cubic spline interpolation across the 390-day period for smooth transitions between months. ' ...
     'STEP 4: Extract middle 365 days (removing the 15-day padding on each side) for final climatology.']);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'interpolation_method', 'cubic_spline');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'padding_days', int32(15));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'extended_period_days', int32(390));

% FORMULA
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'formula', ...
    'UI_Ekman = tau_along / (rho * f)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'formula_components', ...
    'tau_along = along-shore wind stress component; rho = seawater density (1025 kg/m³); f = Coriolis parameter');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'rho_value', 1025);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'rho_units', 'kg/m^3');

% SOURCE AND MASK
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'source_data', 'CCMP v03.1 L4 wind stress');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask_applied', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask_description', ...
    'Offshore distance mask - 0 to 50 km from Canary Current coast');

% TEMPORAL COVERAGE
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'time_period', ...
    sprintf('%d-%d', years(1), years(end)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'years_included', sprintf('%d-%d', years(1), years(end)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'n_years', int32(n_years));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days_handling', 'removed - always 365 days');

% SPATIAL COVERAGE
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'spatial_coverage', ...
    'Canary Current region (Western Africa)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lon_min', double(min(lon)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lon_max', double(max(lon)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lat_min', double(min(lat)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lat_max', double(max(lat)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'n_grid_points', int32(n_lon * n_lat));

% DATA QUALITY
all_clim_valid = climatology(~isnan(climatology));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_min', single(min(all_clim_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_max', single(max(all_clim_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_mean', single(mean(all_clim_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_std', single(std(all_clim_valid)));

% PROCESSING
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'processing_date', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'processed_by', 'MATLAB climatology calculation script');

% REFERENCES
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'references', ...
    'Atlas, R. et al. (2011) The Cross-Calibrated Multi-Platform (CCMP) Ocean Surface Wind Velocity Product');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'notes', ...
    'Positive UI values indicate upwelling-favorable wind conditions. Negative = downwelling-favorable. NaN = masked regions.');

netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.7');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'institution', 'SOTON/NOC');

netcdf.close(ncid);

fprintf('✓ Saved climatology: %s\n', climatology_filepath);

% ========================================================================
% SAVE ANOMALIES TO NETCDF - [lon, lat, time]
% ========================================================================

anomalies_filepath = fullfile(output_path, sprintf('UI_ekman_anomalies_%s.nc', mask_name));
if isfile(anomalies_filepath)
    delete(anomalies_filepath);
end

ncid = netcdf.create(anomalies_filepath, mode);

% Define dimensions: [lon, lat, time]
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
time_dimid = netcdf.defDim(ncid, 'time', file_count);

% Define coordinate variables
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
date_varid = netcdf.defVar(ncid, 'date', 'NC_CHAR', [netcdf.defDim(ncid, 'date_strlen', 8), time_dimid]);

% Define anomalies data variable: [lon, lat, time]
anomalies_varid = netcdf.defVar(ncid, 'UI_ekman_anomalies', 'NC_FLOAT', ...
    [lon_dimid, lat_dimid, time_dimid]);

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, lat_varid, double(lat));

% Create time array (days since 2000-01-01)
time_values = datenum(dates_all, 'yyyymmdd') - datenum('2000-01-01');
netcdf.putVar(ncid, time_varid, double(time_values));

% Write dates
dates_char = char(dates_all);
netcdf.putVar(ncid, date_varid, dates_char');

% Write anomalies
netcdf.putVar(ncid, anomalies_varid, single(anomalies_data));

% Add attributes
netcdf.reDef(ncid);

netcdf.putAtt(ncid, lon_varid, 'standard_name', 'longitude');
netcdf.putAtt(ncid, lon_varid, 'long_name', 'Longitude');
netcdf.putAtt(ncid, lon_varid, 'units', 'degrees_east');

netcdf.putAtt(ncid, lat_varid, 'standard_name', 'latitude');
netcdf.putAtt(ncid, lat_varid, 'long_name', 'Latitude');
netcdf.putAtt(ncid, lat_varid, 'units', 'degrees_north');

netcdf.putAtt(ncid, time_varid, 'standard_name', 'time');
netcdf.putAtt(ncid, time_varid, 'long_name', 'Time');
netcdf.putAtt(ncid, time_varid, 'units', 'days since 2000-01-01 00:00:00');
netcdf.putAtt(ncid, time_varid, 'calendar', 'standard');

netcdf.putAtt(ncid, date_varid, 'long_name', 'Date in YYYYMMDD format');

netcdf.putAtt(ncid, anomalies_varid, 'long_name', 'Ekman Upwelling Index - Daily Anomalies');
netcdf.putAtt(ncid, anomalies_varid, 'units', 'm^2 s^-1');

% ========================================================================
% GLOBAL ATTRIBUTES FOR ANOMALIES
% ========================================================================

netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', ...
    'Ekman Upwelling Index - Daily Anomalies');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'summary', ...
    'Daily anomalies of Ekman Upwelling Index for the Canary Current region');

% CALCULATION METHOD
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'calculation_method', ...
    'DAILY_MEAN - DAILY_CLIMATOLOGY');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'calculation_description', ...
    ['STEP 1: Calculate daily mean UI_Ekman by averaging the 4 daily time steps (observations per day). ' ...
     'STEP 2: For each day, find the corresponding day-of-year in the 365-day climatology. ' ...
     'STEP 3: Subtract the climatological value from the daily mean: anomaly = daily_mean - climatology(DOY). ' ...
     'STEP 4: Result = deviation from long-term normal for that day of year.']);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'interpretation', ...
    'Positive anomalies = stronger than normal upwelling. Negative anomalies = weaker than normal upwelling.');

% REFERENCE CLIMATOLOGY
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'climatology_dataset', ...
    sprintf('UI_ekman_daily_climatology_%s.nc', mask_name));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'climatology_period', ...
    sprintf('%d-%d', years(1), years(end)));

% FORMULA
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'ui_formula', ...
    'UI_Ekman = tau_along / (rho * f)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'rho_value', 1025);

% SOURCE AND MASK
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'source_data', 'CCMP v03.1 L4 wind stress');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask_applied', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask_description', ...
    'Offshore distance mask - 0 to 50 km from Canary Current coast');

% TEMPORAL COVERAGE
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'time_period', ...
    sprintf('%s to %s', dates_all{1}, dates_all{end}));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'n_days', int32(file_count));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'years_covered', ...
    sprintf('%d-%d', years(1), years(end)));

% SPATIAL COVERAGE
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'spatial_coverage', ...
    'Canary Current region (Western Africa)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lon_min', double(min(lon)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lon_max', double(max(lon)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lat_min', double(min(lat)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'lat_max', double(max(lat)));

% DATA QUALITY
all_anom_valid = anomalies_data(~isnan(anomalies_data));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_min', single(min(all_anom_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_max', single(max(all_anom_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_mean', single(mean(all_anom_valid)));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'data_std', single(std(all_anom_valid)));

% PROCESSING
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'processing_date', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'processed_by', 'MATLAB climatology calculation script');

% REFERENCES
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'references', ...
    'Atlas, R. et al. (2011) The Cross-Calibrated Multi-Platform (CCMP) Ocean Surface Wind Velocity Product');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'notes', ...
    'NaN values indicate masked regions outside the offshore domain.');

netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'Conventions', 'CF-1.7');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'institution', 'SOTON/NOC');

netcdf.close(ncid);

fprintf('✓ Saved anomalies: %s\n\n', anomalies_filepath);

% ========================================================================
% SUMMARY
% ========================================================================

fprintf('\n');
fprintf('========================================\n');
fprintf('✓ ALL DATASETS SAVED SUCCESSFULLY!\n');
fprintf('========================================\n');
fprintf('\nOutput Files:\n\n');

fprintf('1. STATISTICS DATASET (.MAT):\n');
fprintf('   Path: %s\n', stats_filepath);
fprintf('   Contents:\n');
fprintf('     - Overall statistics (mean, std, min, max, median)\n');
fprintf('     - Monthly statistics (12 months)\n');
fprintf('     - Yearly statistics (%d years)\n', n_years);
fprintf('     - Full descriptions and interpretations\n\n');

fprintf('2. CLIMATOLOGY (.NC):\n');
fprintf('   Path: %s\n', climatology_filepath);
fprintf('   Dimensions: [lon=%d, lat=%d, time=365]\n', n_lon, n_lat);
fprintf('   Calculation: Monthly mean -> 390-day extended -> Cubic spline interpolation -> 365-day extraction\n');
fprintf('   Period: %d-%d years\n\n', years(1), years(end));

fprintf('3. ANOMALIES (.NC):\n');
fprintf('   Path: %s\n', anomalies_filepath);
fprintf('   Dimensions: [lon=%d, lat=%d, time=%d days]\n', n_lon, n_lat, file_count);
fprintf('   Calculation: Daily mean - Daily climatology (DOY)\n');
fprintf('   Period: %s to %s\n\n', dates_all{1}, dates_all{end});

fprintf('Metadata Documentation:\n');
fprintf('   ✓ Calculation methods fully documented\n');
fprintf('   ✓ Formula: UI_Ekman = tau_along / (rho * f)\n');
fprintf('   ✓ Mask: %s (0-50 km offshore)\n', mask_name);
fprintf('   ✓ Interpolation: Cubic spline with 15-day monthly padding\n');
fprintf('   ✓ Data quality statistics included\n');
fprintf('   ✓ Temporal and spatial coverage documented\n');
fprintf('   ✓ Processing timestamps recorded\n');
fprintf('========================================\n\n');



%% Subset Latitudes

% ========================================================================
% DEFINE LATITUDE SUBSETS
% ========================================================================

fprintf('Defining latitude subsets...\n');

% Find latitude indices for subsets
lat_min_all = min(lat);
lat_max_all = max(lat);
lat_25_idx = find(lat >= lat_min_all & lat < 25);
lat_25_30_idx = find(lat >= 25 & lat < 30);
lat_30_max_idx = find(lat >= 30 & lat <= lat_max_all);

fprintf('  Subset 1: %.2f°N to 25°N - %d lat points\n', lat_min_all, length(lat_25_idx));
fprintf('  Subset 2: 25°N to 30°N - %d lat points\n', length(lat_25_30_idx));
fprintf('  Subset 3: 30°N to %.2f°N - %d lat points\n\n', lat_max_all, length(lat_30_max_idx));

% Define subsets
lat_subsets = struct();
lat_subsets(1).name = sprintf('%.2f-25N', lat_min_all);
lat_subsets(1).lat_range = [lat_min_all, 25];
lat_subsets(1).idx = lat_25_idx;
lat_subsets(2).name = '25-30N';
lat_subsets(2).lat_range = [25, 30];
lat_subsets(2).idx = lat_25_30_idx;
lat_subsets(3).name = sprintf('30-%.2fN', lat_max_all);
lat_subsets(3).lat_range = [30, lat_max_all];
lat_subsets(3).idx = lat_30_max_idx;

n_subsets = length(lat_subsets);

% ========================================================================
% CALCULATE STATISTICS FOR LATITUDE SUBSETS
% ========================================================================

fprintf('Calculating statistics for latitude subsets...\n\n');

% Initialize subset statistics
subset_stats = struct();

for subset_idx = 1:n_subsets
    subset_name = lat_subsets(subset_idx).name;
    lat_idx = lat_subsets(subset_idx).idx;
    
    fprintf('Processing subset: %s\n', subset_name);
    
    % Calculate daily means for this latitude subset
    daily_ui_subset = zeros(file_count, 1);
    daily_anom_subset = zeros(file_count, 1);
    
    for f = 1:file_count
        doy = day_of_year_all(f);
        
        % Extract data for this latitude subset [lon, n_lat_subset, time]
        ui_subset = all_data{f}(:, lat_idx, :);  % [lon, n_lat_subset, time]
        
        % Average over time, then over lon and lat for this subset
        ui_subset_mean = squeeze(mean(ui_subset, 3, 'omitnan'));  % [lon, n_lat_subset]
        daily_ui_subset(f) = mean(ui_subset_mean(~isnan(ui_subset_mean)));
        
        % Anomaly
        clim_subset = climatology(:, lat_idx, doy);  % [lon, n_lat_subset]
        daily_anom_subset(f) = daily_ui_subset(f) - mean(clim_subset(~isnan(clim_subset)));
    end
    
    % Store overall statistics
    subset_stats(subset_idx).name = subset_name;
    subset_stats(subset_idx).lat_range = lat_subsets(subset_idx).lat_range;
    
    % Overall stats
    ui_valid = daily_ui_subset(~isnan(daily_ui_subset));
    anom_valid = daily_anom_subset(~isnan(daily_anom_subset));
    
    subset_stats(subset_idx).overall.ui_mean = mean(ui_valid);
    subset_stats(subset_idx).overall.ui_std = std(ui_valid);
    subset_stats(subset_idx).overall.ui_min = min(ui_valid);
    subset_stats(subset_idx).overall.ui_max = max(ui_valid);
    subset_stats(subset_idx).overall.ui_median = median(ui_valid);
    subset_stats(subset_idx).overall.anom_mean = mean(anom_valid);
    subset_stats(subset_idx).overall.anom_std = std(anom_valid);
    subset_stats(subset_idx).overall.n_observations = length(ui_valid);
    
    % Monthly statistics
    subset_stats(subset_idx).monthly.months = monthly_stats.month_names';
    subset_stats(subset_idx).monthly.ui_mean = zeros(12, 1);
    subset_stats(subset_idx).monthly.ui_std = zeros(12, 1);
    subset_stats(subset_idx).monthly.ui_min = zeros(12, 1);
    subset_stats(subset_idx).monthly.ui_max = zeros(12, 1);
    subset_stats(subset_idx).monthly.ui_median = zeros(12, 1);
    subset_stats(subset_idx).monthly.anom_mean = zeros(12, 1);
    subset_stats(subset_idx).monthly.anom_std = zeros(12, 1);
    subset_stats(subset_idx).monthly.anom_min = zeros(12, 1);
    subset_stats(subset_idx).monthly.anom_max = zeros(12, 1);
    subset_stats(subset_idx).monthly.count = zeros(12, 1);
    
    for m = 1:12
        month_mask = month_all == m;
        month_ui_subset = daily_ui_subset(month_mask);
        month_anom_subset = daily_anom_subset(month_mask);
        
        month_ui_subset = month_ui_subset(~isnan(month_ui_subset));
        month_anom_subset = month_anom_subset(~isnan(month_anom_subset));
        
        if ~isempty(month_ui_subset)
            subset_stats(subset_idx).monthly.ui_mean(m) = mean(month_ui_subset);
            subset_stats(subset_idx).monthly.ui_std(m) = std(month_ui_subset);
            subset_stats(subset_idx).monthly.ui_min(m) = min(month_ui_subset);
            subset_stats(subset_idx).monthly.ui_max(m) = max(month_ui_subset);
            subset_stats(subset_idx).monthly.ui_median(m) = median(month_ui_subset);
            subset_stats(subset_idx).monthly.anom_mean(m) = mean(month_anom_subset);
            subset_stats(subset_idx).monthly.anom_std(m) = std(month_anom_subset);
            subset_stats(subset_idx).monthly.anom_min(m) = min(month_anom_subset);
            subset_stats(subset_idx).monthly.anom_max(m) = max(month_anom_subset);
            subset_stats(subset_idx).monthly.count(m) = length(month_ui_subset);
        else
            subset_stats(subset_idx).monthly.ui_mean(m) = NaN;
            subset_stats(subset_idx).monthly.ui_std(m) = NaN;
            subset_stats(subset_idx).monthly.ui_min(m) = NaN;
            subset_stats(subset_idx).monthly.ui_max(m) = NaN;
            subset_stats(subset_idx).monthly.ui_median(m) = NaN;
            subset_stats(subset_idx).monthly.anom_mean(m) = NaN;
            subset_stats(subset_idx).monthly.anom_std(m) = NaN;
            subset_stats(subset_idx).monthly.anom_min(m) = NaN;
            subset_stats(subset_idx).monthly.anom_max(m) = NaN;
            subset_stats(subset_idx).monthly.count(m) = 0;
        end
    end
    
    % Yearly statistics
    subset_stats(subset_idx).yearly.years = yearly_stats.years;
    subset_stats(subset_idx).yearly.ui_mean = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.ui_std = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.ui_min = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.ui_max = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.anom_mean = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.anom_std = zeros(n_years, 1);
    subset_stats(subset_idx).yearly.count = zeros(n_years, 1);
    
    for y = 1:n_years
        year_mask = year(dates_datetime) == years(y);
        year_ui_subset = daily_ui_subset(year_mask);
        year_anom_subset = daily_anom_subset(year_mask);
        
        year_ui_subset = year_ui_subset(~isnan(year_ui_subset));
        year_anom_subset = year_anom_subset(~isnan(year_anom_subset));
        
        if ~isempty(year_ui_subset)
            subset_stats(subset_idx).yearly.ui_mean(y) = mean(year_ui_subset);
            subset_stats(subset_idx).yearly.ui_std(y) = std(year_ui_subset);
            subset_stats(subset_idx).yearly.ui_min(y) = min(year_ui_subset);
            subset_stats(subset_idx).yearly.ui_max(y) = max(year_ui_subset);
            subset_stats(subset_idx).yearly.anom_mean(y) = mean(year_anom_subset);
            subset_stats(subset_idx).yearly.anom_std(y) = std(year_anom_subset);
            subset_stats(subset_idx).yearly.count(y) = length(year_ui_subset);
        else
            subset_stats(subset_idx).yearly.ui_mean(y) = NaN;
            subset_stats(subset_idx).yearly.ui_std(y) = NaN;
            subset_stats(subset_idx).yearly.ui_min(y) = NaN;
            subset_stats(subset_idx).yearly.ui_max(y) = NaN;
            subset_stats(subset_idx).yearly.anom_mean(y) = NaN;
            subset_stats(subset_idx).yearly.anom_std(y) = NaN;
            subset_stats(subset_idx).yearly.count(y) = 0;
        end
    end
    
    fprintf('  ✓ Statistics calculated for %s\n\n', subset_name);
end

% ========================================================================
% DISPLAY SUBSET STATISTICS
% ========================================================================

fprintf('LATITUDE SUBSET STATISTICS SUMMARY\n');
fprintf('========================================\n\n');

for subset_idx = 1:n_subsets
    subset_name = subset_stats(subset_idx).name;
    lat_range = subset_stats(subset_idx).lat_range;
    
    fprintf('SUBSET %d: %s (%.2f°N to %.2f°N)\n', subset_idx, subset_name, lat_range(1), lat_range(2));
    fprintf('========================================\n');
    
    fprintf('OVERALL STATISTICS:\n');
    fprintf('  UI Mean:        %.4f m²/s\n', subset_stats(subset_idx).overall.ui_mean);
    fprintf('  UI Std:         %.4f m²/s\n', subset_stats(subset_idx).overall.ui_std);
    fprintf('  UI Min:         %.4f m²/s\n', subset_stats(subset_idx).overall.ui_min);
    fprintf('  UI Max:         %.4f m²/s\n', subset_stats(subset_idx).overall.ui_max);
    fprintf('  UI Median:      %.4f m²/s\n', subset_stats(subset_idx).overall.ui_median);
    fprintf('  Anom Mean:      %.4f m²/s\n', subset_stats(subset_idx).overall.anom_mean);
    fprintf('  Anom Std:       %.4f m²/s\n', subset_stats(subset_idx).overall.anom_std);
    fprintf('  Observations:   %d\n\n', subset_stats(subset_idx).overall.n_observations);
    
    fprintf('MONTHLY STATISTICS:\n');
    fprintf('%-12s %10s %10s %10s %10s %10s | %10s %10s %10s %10s | %6s\n', ...
        'Month', 'UI Mean', 'UI Std', 'UI Min', 'UI Max', 'UI Median', ...
        'Anom Mean', 'Anom Std', 'Anom Min', 'Anom Max', 'Count');
    fprintf(repmat('-', 1, 130));
    fprintf('\n');
    
    for m = 1:12
        fprintf('%-12s %10.4f %10.4f %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f %10.4f | %6d\n', ...
            subset_stats(subset_idx).monthly.months{m}, ...
            subset_stats(subset_idx).monthly.ui_mean(m), subset_stats(subset_idx).monthly.ui_std(m), ...
            subset_stats(subset_idx).monthly.ui_min(m), subset_stats(subset_idx).monthly.ui_max(m), ...
            subset_stats(subset_idx).monthly.ui_median(m), ...
            subset_stats(subset_idx).monthly.anom_mean(m), subset_stats(subset_idx).monthly.anom_std(m), ...
            subset_stats(subset_idx).monthly.anom_min(m), subset_stats(subset_idx).monthly.anom_max(m), ...
            subset_stats(subset_idx).monthly.count(m));
    end
    fprintf('\n');
end

fprintf('========================================\n\n');

% ========================================================================
% CREATE EXTENDED STATISTICS DATASET WITH SUBSETS
% ========================================================================

fprintf('Creating extended statistics dataset with latitude subsets...\n');

statistics_dataset_extended = statistics_dataset;  % Start with original stats

% Add subset statistics
statistics_dataset_extended.latitude_subsets = struct();
statistics_dataset_extended.latitude_subsets.n_subsets = n_subsets;
statistics_dataset_extended.latitude_subsets.subset_definitions = lat_subsets;

for subset_idx = 1:n_subsets
    subset_name = ['subset_' num2str(subset_idx)];
    
    statistics_dataset_extended.latitude_subsets.(subset_name) = struct();
    statistics_dataset_extended.latitude_subsets.(subset_name).name = subset_stats(subset_idx).name;
    statistics_dataset_extended.latitude_subsets.(subset_name).lat_range = subset_stats(subset_idx).lat_range;
    statistics_dataset_extended.latitude_subsets.(subset_name).lat_range_description = ...
        sprintf('%.2f°N to %.2f°N', subset_stats(subset_idx).lat_range(1), subset_stats(subset_idx).lat_range(2));
    
    % Overall stats with descriptions
    statistics_dataset_extended.latitude_subsets.(subset_name).overall = struct();
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.ui_mean = struct(...
        'value', subset_stats(subset_idx).overall.ui_mean, ...
        'description', 'Mean UI_Ekman for this latitude subset', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.ui_std = struct(...
        'value', subset_stats(subset_idx).overall.ui_std, ...
        'description', 'Standard deviation of UI_Ekman', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.ui_min = struct(...
        'value', subset_stats(subset_idx).overall.ui_min, ...
        'description', 'Minimum UI_Ekman (most downwelling-favorable)', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.ui_max = struct(...
        'value', subset_stats(subset_idx).overall.ui_max, ...
        'description', 'Maximum UI_Ekman (most upwelling-favorable)', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.ui_median = struct(...
        'value', subset_stats(subset_idx).overall.ui_median, ...
        'description', 'Median UI_Ekman value', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.anom_mean = struct(...
        'value', subset_stats(subset_idx).overall.anom_mean, ...
        'description', 'Mean anomaly', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.anom_std = struct(...
        'value', subset_stats(subset_idx).overall.anom_std, ...
        'description', 'Standard deviation of anomalies', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).overall.n_observations = ...
        subset_stats(subset_idx).overall.n_observations;
    
    % Monthly stats with descriptions
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly = struct();
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.months = subset_stats(subset_idx).monthly.months;
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.ui_mean = struct(...
        'data', subset_stats(subset_idx).monthly.ui_mean, ...
        'description', 'Mean UI_Ekman for each calendar month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.ui_std = struct(...
        'data', subset_stats(subset_idx).monthly.ui_std, ...
        'description', 'Standard deviation within each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.ui_min = struct(...
        'data', subset_stats(subset_idx).monthly.ui_min, ...
        'description', 'Minimum UI_Ekman in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.ui_max = struct(...
        'data', subset_stats(subset_idx).monthly.ui_max, ...
        'description', 'Maximum UI_Ekman in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.ui_median = struct(...
        'data', subset_stats(subset_idx).monthly.ui_median, ...
        'description', 'Median UI_Ekman in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.anom_mean = struct(...
        'data', subset_stats(subset_idx).monthly.anom_mean, ...
        'description', 'Mean anomaly for each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.anom_std = struct(...
        'data', subset_stats(subset_idx).monthly.anom_std, ...
        'description', 'Standard deviation of anomalies in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.anom_min = struct(...
        'data', subset_stats(subset_idx).monthly.anom_min, ...
        'description', 'Minimum anomaly in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.anom_max = struct(...
        'data', subset_stats(subset_idx).monthly.anom_max, ...
        'description', 'Maximum anomaly in each month', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).monthly.count = struct(...
        'data', subset_stats(subset_idx).monthly.count, ...
        'description', 'Number of observations in each month', ...
        'units', 'count');
    
    % Yearly stats with descriptions
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly = struct();
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.years = subset_stats(subset_idx).yearly.years;
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.ui_mean = struct(...
        'data', subset_stats(subset_idx).yearly.ui_mean, ...
        'description', 'Mean UI_Ekman for each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.ui_std = struct(...
        'data', subset_stats(subset_idx).yearly.ui_std, ...
        'description', 'Standard deviation within each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.ui_min = struct(...
        'data', subset_stats(subset_idx).yearly.ui_min, ...
        'description', 'Minimum UI_Ekman in each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.ui_max = struct(...
        'data', subset_stats(subset_idx).yearly.ui_max, ...
        'description', 'Maximum UI_Ekman in each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.anom_mean = struct(...
        'data', subset_stats(subset_idx).yearly.anom_mean, ...
        'description', 'Mean anomaly for each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.anom_std = struct(...
        'data', subset_stats(subset_idx).yearly.anom_std, ...
        'description', 'Standard deviation of anomalies within each year', ...
        'units', 'm^2 s^-1');
    statistics_dataset_extended.latitude_subsets.(subset_name).yearly.count = struct(...
        'data', subset_stats(subset_idx).yearly.count, ...
        'description', 'Number of observations in each year', ...
        'units', 'count');
end

fprintf('Extended statistics dataset created with latitude subsets\n\n');

% ========================================================================
% SAVE EXTENDED STATISTICS DATASET
% ========================================================================

fprintf('Saving extended statistics dataset...\n');

stats_extended_filepath = fullfile(output_path, sprintf('UI_ekman_statistics_with_subsets_%s.mat', mask_name));
save(stats_extended_filepath, 'statistics_dataset_extended', '-v7.3');

fprintf('Saved extended statistics: %s\n\n', stats_extended_filepath);

fprintf('========================================\n');
fprintf('LATITUDE SUBSET ANALYSIS COMPLETE!\n');


% ========================================================================
% CREATE PLOT 1: TIME SERIES WITH GREY/PINK SHADING AND BLACK LINE
% ========================================================================

fprintf('Creating plot 1: Time series with grey/pink shading...\n\n');

set(0, 'DefaultFigureVisible', 'off');

% Define subset names for plot titles
subset_names_plot = {
    sprintf('Northern Area (30°N - %.2f°N)', lat_max_all);
    'Middle Area (25°N - 30°N)';
    sprintf('Southern Area (%.2f°N - 25°N)', lat_min_all)
};

% Reorder the subsets: Northern (3), Middle (2), Southern (1)
subset_order = [3, 2, 1];

fig1 = figure('Position', [100 100 1600 1000], 'Color', 'w', 'Visible', 'off');

for plot_idx = 1:n_subsets
    subplot(3, 1, plot_idx);
    
    % Get the actual subset index based on reordered plot position
    subset_idx = subset_order(plot_idx);
    
    % Get latitude indices for this subset
    lat_idx = lat_subsets(subset_idx).idx;
    
    % Calculate daily UI for this subset
    daily_ui_subset = zeros(file_count, 1);
    
    for f = 1:file_count
        % Extract data for this latitude subset
        ui_subset = all_data{f}(:, lat_idx, :);  % [lon, n_lat_subset, time]
        ui_subset_mean = squeeze(mean(ui_subset, 3, 'omitnan'));  % [lon, n_lat_subset]
        daily_ui_subset(f) = mean(ui_subset_mean(~isnan(ui_subset_mean)));
    end
    
    hold on
    
    % ====================================================================
    % Fill shaded regions: light grey for negative, light red for positive
    % ====================================================================
    for f = 1:file_count-1
        y1 = 0;
        y2_ui = daily_ui_subset(f);
        y2_next_ui = daily_ui_subset(f+1);
        
        if y2_ui >= 0
            % Positive - light red
            fill([dates_datetime(f), dates_datetime(f+1), dates_datetime(f+1), dates_datetime(f)], ...
                 [y1, y1, y2_next_ui, y2_ui], ...
                 [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        else
            % Negative - light grey
            fill([dates_datetime(f), dates_datetime(f+1), dates_datetime(f+1), dates_datetime(f)], ...
                 [y1, y1, y2_next_ui, y2_ui], ...
                 [0.85, 0.85, 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        end
    end
    
    % Plot zero line
    plot(dates_datetime([1, end]), [0, 0], 'k-', 'LineWidth', 1.5);
    
    % ====================================================================
    % Plot time series with thin black line
    % ====================================================================
    plot(dates_datetime, daily_ui_subset, 'k-', 'LineWidth', 0.8);
    
    % ====================================================================
    % Formatting
    % ====================================================================
    
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('UI_{Ekman} (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
    title(subset_names_plot{plot_idx}, 'FontSize', 12, 'FontWeight', 'bold');
    grid on
    
    % Format x-axis: MMM/YYYY (e.g., Jan/2006) horizontal
    datetick(gca, 'x', 'mmm/yyyy', 'keeplimits', 'keepticks');
    set(gca, 'XTickLabelRotation', 0);
    
    hold off
end

xlabel('Time', 'FontSize', 11, 'FontWeight', 'bold');
sgtitle('UI_{Ekman}: Daily Time Series', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(fig1, fullfile(plot_path, 'UI_subsets_timeseries.png'), 'png');
saveas(fig1, fullfile(plot_path, 'UI_subsets_timeseries.fig'), 'fig');
print(fig1, fullfile(plot_path, 'UI_subsets_timeseries.jpg'), '-djpeg', '-r150');
close(fig1);

fprintf('✓ Saved plot 1: UI_subsets_timeseries\n');
fprintf('  - Black thin line: Daily mean UI_Ekman\n');
fprintf('  - Light red shading: Positive values (upwelling-favorable)\n');
fprintf('  - Light grey shading: Negative values (downwelling-favorable)\n');
fprintf('  - X-axis: MMM/YYYY format, horizontal orientation\n\n');

% ========================================================================
% CREATE PLOT 1: HOVMÖLLER DIAGRAM (TIME vs LATITUDE)
% ========================================================================

fprintf('Creating plot 1: Hovmöller diagram (time vs latitude)...\n\n');

set(0, 'DefaultFigureVisible', 'off');

% Calculate anomalies for all grid points
anomalies_full = zeros(n_lon, n_lat, file_count, 'single');

for f = 1:file_count
    doy = day_of_year_all(f);
    
    for i = 1:n_lon
        for j = 1:n_lat
            ui_point = all_data{f}(i, j, :);
            ui_mean = mean(ui_point(~isnan(ui_point)), 'omitnan');
            clim_point = climatology(i, j, doy);
            anomalies_full(i, j, f) = ui_mean - clim_point;
        end
    end
end

% Average over longitude for Hovmöller
hovmoller_data = squeeze(mean(anomalies_full, 1, 'omitnan'));  % [lat, time]

% Create figure
fig1 = figure('Position', [100 100 1400 800], 'Color', 'w', 'Visible', 'off');

% Plot Hovmöller
imagesc(1:file_count, lat, hovmoller_data);
colorbar_h = colorbar;
ylabel(colorbar_h, 'Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');

% Set colormap
colormap(gca, 'jet');
caxis([-4 4]);

% Formatting
set(gca, 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Latitude (°N)', 'FontSize', 12, 'FontWeight', 'bold');
title('Hovmöller Diagram: UI_{Ekman} Anomalies (Longitude-averaged)', ...
    'FontSize', 13, 'FontWeight', 'bold');
grid on

% ====================================================================
% Add time markers: January and June of every year
% ====================================================================

% Find indices for January and June of each year
time_tick_indices = [];
time_tick_labels = {};

for y = 1:n_years
    current_year = years(y);
    
    % Find January of this year
    jan_idx = find(year(dates_datetime) == current_year & month(dates_datetime) == 1, 1, 'first');
    if ~isempty(jan_idx)
        time_tick_indices = [time_tick_indices, jan_idx];
        time_tick_labels{end+1} = sprintf('Jan/%d', current_year);
    end
    
    % Find June of this year
    jun_idx = find(year(dates_datetime) == current_year & month(dates_datetime) == 6, 1, 'first');
    if ~isempty(jun_idx)
        time_tick_indices = [time_tick_indices, jun_idx];
        time_tick_labels{end+1} = sprintf('Jun/%d', current_year);
    end
end

% Set time ticks and labels
set(gca, 'XTick', time_tick_indices);
set(gca, 'XTickLabel', time_tick_labels);
set(gca, 'XTickLabelRotation', 0);

% Flip Y-axis to show lower latitudes at bottom and higher at top
set(gca, 'YDir', 'normal');

saveas(fig1, fullfile(plot_path, 'UI_hovmoller.png'), 'png');
saveas(fig1, fullfile(plot_path, 'UI_hovmoller.fig'), 'fig');
print(fig1, fullfile(plot_path, 'UI_hovmoller.jpg'), '-djpeg', '-r150');
close(fig1);

fprintf('✓ Saved plot 1: Hovmöller diagram\n');
fprintf('  - Latitude orientation: Low (bottom) → High (top)\n');
fprintf('  - Time markers: January and June of every year\n');
fprintf('  - Y-axis direction: Normal (lower to higher)\n\n');


% ========================================================================
% CREATE PLOT 4: TIME SERIES BY LATITUDE BAND
% ========================================================================

fprintf('Creating plot 4: Time series by latitude band...\n\n');

fig4 = figure('Position', [100 100 1600 1000], 'Color', 'w', 'Visible', 'off');

% Create 3 subplots for the 3 latitude subsets
subset_order = [3, 2, 1];

for plot_idx = 1:n_subsets
    subplot(3, 1, plot_idx);
    
    subset_idx = subset_order(plot_idx);
    lat_idx = lat_subsets(subset_idx).idx;
    
    % Extract time series for this latitude band from hovmoller
    lat_subset_ts = mean(hovmoller_data(lat_idx, :), 1);  % Average over latitude band
    
    hold on
    
    % Fill shading: negative grey, positive pink
    for f = 1:file_count-1
        y1 = 0;
        y2 = lat_subset_ts(f);
        y2_next = lat_subset_ts(f+1);
        
        if y2 >= 0
            % Positive - light pink
            fill([dates_datetime(f), dates_datetime(f+1), dates_datetime(f+1), dates_datetime(f)], ...
                 [y1, y1, y2_next, y2], ...
                 [1, 0.75, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        else
            % Negative - light grey
            fill([dates_datetime(f), dates_datetime(f+1), dates_datetime(f+1), dates_datetime(f)], ...
                 [y1, y1, y2_next, y2], ...
                 [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end
    end
    
    % Plot zero line
    plot(dates_datetime([1, end]), [0, 0], 'k-', 'LineWidth', 1.5);
    
    % Plot time series with thin black line
    plot(dates_datetime, lat_subset_ts, 'k-', 'LineWidth', 0.8);
    
    % Formatting
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s - Anomaly Time Series', subset_names_plot{plot_idx}), ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on
    
    % Format x-axis
    datetick(gca, 'x', 'mmm/yyyy', 'keeplimits', 'keepticks');
    set(gca, 'XTickLabelRotation', 0);
    
    hold off
end

xlabel('Time', 'FontSize', 11, 'FontWeight', 'bold');
sgtitle('Anomaly Time Series by Latitude Subset', 'FontSize', 13, 'FontWeight', 'bold');

saveas(fig4, fullfile(plot_path, 'UI_anomaly_timeseries_by_latitude.png'), 'png');
saveas(fig4, fullfile(plot_path, 'UI_anomaly_timeseries_by_latitude.fig'), 'fig');
print(fig4, fullfile(plot_path, 'UI_anomaly_timeseries_by_latitude.jpg'), '-djpeg', '-r150');
close(fig4);

fprintf('Saved plot 4: Time series by latitude band\n\n');

% ========================================================================
% CREATE PLOT 5: MONTHLY BOX PLOT (DETAILED VERSION)
% ========================================================================

fprintf('Creating plot 5: Monthly box plot with variability analysis...\n\n');

fig5 = figure('Position', [100 100 1600 1000], 'Color', 'w', 'Visible', 'off');

subset_order = [3, 2, 1];

for plot_idx = 1:n_subsets
    subplot(3, 1, plot_idx);
    
    subset_idx = subset_order(plot_idx);
    lat_idx = lat_subsets(subset_idx).idx;
    
    % Prepare data for each month
    bp_data_matrix = [];
    bp_groups = [];
    monthly_variability = zeros(12, 1);  % IQR (variability)
    
    for m = 1:12
        month_mask = month_all == m;
        month_files = find(month_mask);
        
        month_anom = [];
        for f = month_files'
            doy = day_of_year_all(f);
            
            ui_subset = all_data{f}(:, lat_idx, :);
            ui_subset_mean = squeeze(mean(ui_subset, 3, 'omitnan'));
            daily_ui = mean(ui_subset_mean(~isnan(ui_subset_mean)));
            
            clim_subset = climatology(:, lat_idx, doy);
            daily_clim = mean(clim_subset(~isnan(clim_subset)));
            
            month_anom = [month_anom; daily_ui - daily_clim];
        end
        
        % Append data
        bp_data_matrix = [bp_data_matrix; month_anom];
        bp_groups = [bp_groups; repmat(m, length(month_anom), 1)];
        
        % Calculate IQR (Interquartile Range) - measure of variability
        q75 = quantile(month_anom, 0.75);
        q25 = quantile(month_anom, 0.25);
        monthly_variability(m) = q75 - q25;
    end
    
    % Create box plot with correct data format
    hold on
    
    month_labels_box = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
                        'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
    
    bp = boxplot(bp_data_matrix, bp_groups, 'Labels', month_labels_box, ...
                 'Positions', 1:12, 'Width', 0.6, 'OutlierSize', 5);
    
    % Color the boxes with gradient based on variability
    for i = 1:12
        % Normalize variability to 0-1
        norm_var = (monthly_variability(i) - min(monthly_variability)) / ...
                   (max(monthly_variability) - min(monthly_variability) + eps);
        
        % Color: more variable = more red, less variable = more blue
        box_color = [norm_var, 0.7*(1-norm_var), 1-norm_var];
        
        patch(get(bp(5,i), 'XData'), get(bp(5,i), 'YData'), ...
              box_color, 'FaceAlpha', 0.7);
    end
    
    % Plot zero line
    plot([0.5, 12.5], [0, 0], 'r--', 'LineWidth', 2);
    
    % Formatting
    set(gca, 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s - Monthly Anomaly Distribution (Colored by Variability)', subset_names_plot{plot_idx}), ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on
    xlim([0.5, 12.5]);
    set(gca, 'XTickLabelRotation', 0);
    
    % Add text annotation
    text(0.02, 0.98, sprintf('More Variable (Red) → Less Variable (Blue)'), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'white', 'FontSize', 9);
    
    hold off
end

xlabel('Month', 'FontSize', 11, 'FontWeight', 'bold');
sgtitle('Monthly Distribution of Anomalies by Latitude Subset (Color = Variability)', ...
    'FontSize', 13, 'FontWeight', 'bold');

saveas(fig5, fullfile(plot_path, 'UI_monthly_boxplot.png'), 'png');
saveas(fig5, fullfile(plot_path, 'UI_monthly_boxplot.fig'), 'fig');
print(fig5, fullfile(plot_path, 'UI_monthly_boxplot.jpg'), '-djpeg', '-r150');
close(fig5);

fprintf('✓ Saved plot 5: Monthly box plot\n');
fprintf('  - Box size = intra-monthly variability\n');
fprintf('  - Red boxes = più variabili\n');
fprintf('  - Blue boxes = meno variabili\n');
fprintf('  - Outliers (●) = valori estremi\n\n');