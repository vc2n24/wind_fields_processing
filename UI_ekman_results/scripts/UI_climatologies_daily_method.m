% ========================================================================
% SCRIPT: TRUE DAILY CLIMATOLOGY WITH GRIDPOINT-LEVEL DATA
% Latitude band: 25-29°N
% With spatial anomaly maps
% ========================================================================
clear all; close all; clc;

ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman';
output_path = fullfile(ui_results_path, 'analysis_true_daily');
plot_path = fullfile(output_path, 'plots');

if ~isfolder(output_path), mkdir(output_path); end
if ~isfolder(plot_path), mkdir(plot_path); end

mask_name = '50km';
years = 2005:2024;
n_years = length(years);
lat_min = 25;
lat_max = 29;

fprintf('========================================================================\n');
fprintf('TRUE DAILY CLIMATOLOGY ANALYSIS WITH GRIDPOINT DATA\n');
fprintf('========================================================================\n');
fprintf('Analysis Parameters:\n');
fprintf('  Latitude band: %.1f - %.1f°N\n', lat_min, lat_max);
fprintf('  Years: %d - %d\n', years(1), years(end));
fprintf('  Mask: %s\n\n', mask_name);

% ========================================================================
% LOAD ALL UI DATA
% ========================================================================

fprintf('Loading UI_ekman data for 2005-2024...\n');

all_data = {};
dates_all = {};
file_count = 0;

for y = 1:n_years
    current_year = years(y);
    Ystr = sprintf('Y%04d', current_year);
    
    for m = 1:12
        current_month = m;
        Mstr = sprintf('M%02d', current_month);
        
        if current_month == 2
            n_days = 28;
        elseif ismember(current_month, [4, 6, 9, 11])
            n_days = 30;
        else
            n_days = 31;
        end
        
        for d = 1:n_days
            year_str = sprintf('%04d', current_year);
            month_str = sprintf('%02d', current_month);
            day_str = sprintf('%02d', d);
            date_str = [year_str month_str day_str];
            
            input_filename = sprintf('UI_ekman_%s%s%s_%s.nc', year_str, month_str, day_str, mask_name);
            ui_dir = sprintf('UI_ekman_%s', mask_name);
            input_filepath = fullfile(ui_results_path, ui_dir, Ystr, Mstr, input_filename);
            
            if ~isfile(input_filepath)
                continue;
            end
            
            try
                UI_data = ncread(input_filepath, 'UI_ekman');
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
% FIND LATITUDE INDICES FOR BAND 25-29°N
% ========================================================================

fprintf('Finding latitude band %.1f - %.1f°N...\n', lat_min, lat_max);

lat_idx = find(lat >= lat_min & lat <= lat_max);
if isempty(lat_idx)
    error('No data found in latitude band %.1f - %.1f°N', lat_min, lat_max);
end

n_lat_band = length(lat_idx);
fprintf('✓ Found %d latitude points in band\n', n_lat_band);
fprintf('  Latitude indices: %d to %d\n\n', lat_idx(1), lat_idx(end));

% ========================================================================
% CALCULATE DAY OF YEAR FOR EACH FILE
% ========================================================================

fprintf('Calculating day of year for each file...\n');

day_of_year_all = zeros(file_count, 1);
month_all = zeros(file_count, 1);
year_all = zeros(file_count, 1);

for f = 1:file_count
    date_str = dates_all{f};
    year_all(f) = str2double(date_str(1:4));
    month_val = str2double(date_str(5:6));
    day_val = str2double(date_str(7:8));
    
    month_all(f) = month_val;
    
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    day_of_year = sum(days_in_months(1:month_val-1)) + day_val;
    day_of_year_all(f) = day_of_year;
end

fprintf('✓ Day of year calculated\n\n');

% ========================================================================
% CALCULATE TRUE DAILY CLIMATOLOGY (365 DAYS)
% Including both LAT-BAND AVERAGED and GRIDPOINT-LEVEL
% ========================================================================

fprintf('Calculating TRUE daily climatology (day-by-day from all years)...\n');

% Initialize: store observations by day
daily_data_pool = cell(365, 1);
daily_years_pool = cell(365, 1);

% Initialize gridpoint-level storage
daily_gridpoint_pool = cell(365, n_lon, n_lat_band);

for doy = 1:365
    daily_data_pool{doy} = [];
    daily_years_pool{doy} = [];
    for ilon = 1:n_lon
        for ilat = 1:n_lat_band
            daily_gridpoint_pool{doy, ilon, ilat} = [];
        end
    end
end

% Collect all observations by day-of-year from latitude band
fprintf('  Pooling all observations by day-of-year...\n');
for f = 1:file_count
    doy = day_of_year_all(f);
    ui_data = all_data{f};  % [lon, lat, time]
    
    % Extract latitude band: [lon, lat_band, time]
    ui_lat_band = ui_data(:, lat_idx, :);  % [n_lon, n_lat_band, n_time]
    
    % Store lat-band average (for overall climatology)
    ui_reshaped = reshape(ui_lat_band, [], size(ui_lat_band, 3));
    daily_data_pool{doy} = [daily_data_pool{doy}, ui_reshaped];
    daily_years_pool{doy} = [daily_years_pool{doy}; repmat(year_all(f), size(ui_reshaped, 2), 1)];
    
    % Store gridpoint-level data
    for ilon = 1:n_lon
        for ilat = 1:n_lat_band
            % Get all time steps for this gridpoint on this day
            gridpoint_values = squeeze(ui_lat_band(ilon, ilat, :));
            daily_gridpoint_pool{doy, ilon, ilat} = [daily_gridpoint_pool{doy, ilon, ilat}; gridpoint_values'];
        end
    end
    
    if mod(f, 200) == 0
        fprintf('    Processed file %d/%d\n', f, file_count);
    end
end

% Initialize LAT-BAND averaged climatology
daily_climatology_raw = zeros(365, 1);
daily_climatology_raw_std = zeros(365, 1);
daily_obs_count = zeros(365, 1);
daily_years_count = zeros(365, 1);

% Initialize GRIDPOINT-LEVEL climatology
daily_gridpoint_clim = zeros(365, n_lon, n_lat_band);
daily_gridpoint_std = zeros(365, n_lon, n_lat_band);
daily_gridpoint_count = zeros(365, n_lon, n_lat_band);

% Compute climatologies
fprintf('  Computing daily means, std, and observation counts...\n');
for doy = 1:365
    % LAT-BAND AVERAGED
    all_obs = daily_data_pool{doy};
    years_obs = daily_years_pool{doy};
    
    if ~isempty(all_obs)
        all_obs_flat = all_obs(:);
        all_obs_flat = all_obs_flat(~isnan(all_obs_flat));
        
        if ~isempty(all_obs_flat)
            daily_climatology_raw(doy) = mean(all_obs_flat, 'omitnan');
            daily_climatology_raw_std(doy) = std(all_obs_flat, 'omitnan');
            unique_years = unique(years_obs);
            daily_years_count(doy) = length(unique_years);
            daily_obs_count(doy) = length(all_obs_flat);
        else
            daily_climatology_raw(doy) = NaN;
            daily_climatology_raw_std(doy) = NaN;
            daily_obs_count(doy) = 0;
            daily_years_count(doy) = 0;
        end
    else
        daily_climatology_raw(doy) = NaN;
        daily_climatology_raw_std(doy) = NaN;
        daily_obs_count(doy) = 0;
        daily_years_count(doy) = 0;
    end
    
    % GRIDPOINT-LEVEL
    for ilon = 1:n_lon
        for ilat = 1:n_lat_band
            gridpoint_obs = daily_gridpoint_pool{doy, ilon, ilat};
            gridpoint_obs = gridpoint_obs(~isnan(gridpoint_obs));
            
            if ~isempty(gridpoint_obs)
                daily_gridpoint_clim(doy, ilon, ilat) = mean(gridpoint_obs, 'omitnan');
                daily_gridpoint_std(doy, ilon, ilat) = std(gridpoint_obs, 'omitnan');
                daily_gridpoint_count(doy, ilon, ilat) = length(gridpoint_obs);
            else
                daily_gridpoint_clim(doy, ilon, ilat) = NaN;
                daily_gridpoint_std(doy, ilon, ilat) = NaN;
                daily_gridpoint_count(doy, ilon, ilat) = 0;
            end
        end
    end
    
    if mod(doy, 50) == 0
        fprintf('    DOY %3d: %.4f ± %.4f m²/s (%d data points, %d years)\n', ...
            doy, daily_climatology_raw(doy), daily_climatology_raw_std(doy), ...
            daily_obs_count(doy), daily_years_count(doy));
    end
end

fprintf('✓ True daily climatology calculated for 365 days\n');
fprintf('  Average years per day: %.1f\n', mean(daily_years_count(daily_years_count > 0)));
fprintf('  Min years per day: %d\n', min(daily_years_count(daily_years_count > 0)));
fprintf('  Max years per day: %d\n', max(daily_years_count));
fprintf('  Average data points per day: %.0f\n', mean(daily_obs_count(daily_obs_count > 0)));
fprintf('  Gridpoint climatology shape: [%d days, %d lon, %d lat]\n\n', ...
    size(daily_gridpoint_clim, 1), size(daily_gridpoint_clim, 2), size(daily_gridpoint_clim, 3));

clear daily_data_pool daily_years_pool daily_gridpoint_pool;

% ========================================================================
% CALCULATE SMOOTHED DAILY CLIMATOLOGY
% ========================================================================

fprintf('Creating smoothed daily climatology...\n');

% Gaussian smoothing directly on 365 days (LAT-BAND)
climatology_smooth_mean = smoothdata(daily_climatology_raw, 'gaussian', 15);
climatology_smooth_std = smoothdata(daily_climatology_raw_std, 'gaussian', 15);

fprintf('  Gaussian smoothing applied (15-day window) to lat-band average\n');

% Use smoothed climatology directly
daily_climatology = climatology_smooth_mean;
daily_climatology_std = climatology_smooth_std;

fprintf('✓ Smoothed daily climatology created: %d days\n\n', length(daily_climatology));

% Gaussian smoothing for GRIDPOINTS
fprintf('  Applying Gaussian smoothing to gridpoint climatology...\n');
daily_gridpoint_clim_smooth = zeros(365, n_lon, n_lat_band);
daily_gridpoint_std_smooth = zeros(365, n_lon, n_lat_band);

for ilon = 1:n_lon
    for ilat = 1:n_lat_band
        % Smooth mean
        clim_profile = squeeze(daily_gridpoint_clim(:, ilon, ilat));
        daily_gridpoint_clim_smooth(:, ilon, ilat) = smoothdata(clim_profile, 'gaussian', 15);
        
        % Smooth std
        std_profile = squeeze(daily_gridpoint_std(:, ilon, ilat));
        daily_gridpoint_std_smooth(:, ilon, ilat) = smoothdata(std_profile, 'gaussian', 15);
    end
    
    if mod(ilon, 50) == 0
        fprintf('    Processed %d/%d longitude points\n', ilon, n_lon);
    end
end

fprintf('✓ Gridpoint climatology smoothed\n\n');

% ========================================================================
% CALCULATE MONTHLY CLIMATOLOGY (FOR REFERENCE)
% ========================================================================

fprintf('Calculating monthly climatology for reference...\n');

% Initialize: store ALL observations for each month
monthly_data_pool = cell(12, 1);
for m = 1:12
    monthly_data_pool{m} = [];
end

% Collect observations from latitude band by month
fprintf('  Pooling observations by month...\n');
for f = 1:file_count
    m = month_all(f);
    ui_data = all_data{f};
    ui_lat_band = ui_data(:, lat_idx, :);
    ui_reshaped = reshape(ui_lat_band, [], size(ui_lat_band, 3));
    monthly_data_pool{m} = [monthly_data_pool{m}, ui_reshaped];
    
    if mod(f, 500) == 0
        fprintf('    Processed file %d/%d\n', f, file_count);
    end
end

% Compute monthly climatology
monthly_climatology = zeros(12, 1);
monthly_climatology_std = zeros(12, 1);

fprintf('  Computing monthly means and std...\n');
for m = 1:12
    all_obs = monthly_data_pool{m};
    
    if ~isempty(all_obs)
        all_obs_flat = all_obs(:);
        all_obs_flat = all_obs_flat(~isnan(all_obs_flat));
        
        monthly_climatology(m) = mean(all_obs_flat, 'omitnan');
        monthly_climatology_std(m) = std(all_obs_flat, 'omitnan');
        n_obs = length(all_obs_flat);
        fprintf('    Month %2d: %.4f ± %.4f m²/s (%d obs)\n', m, monthly_climatology(m), monthly_climatology_std(m), n_obs);
    else
        monthly_climatology(m) = NaN;
        monthly_climatology_std(m) = NaN;
        fprintf('    Month %2d: No data\n', m);
    end
end

fprintf('✓ Monthly climatology calculated\n\n');

clear monthly_data_pool;

% ========================================================================
% CALCULATE DAILY ANOMALIES AND STATISTICS (LAT-BAND AVERAGED)
% ========================================================================

fprintf('Calculating daily statistics and anomalies (lat-band averaged)...\n');

daily_observed_mean = zeros(file_count, 1);
daily_anom_from_daily_clim = zeros(file_count, 1);
dates_datetime = zeros(file_count, 1);
month_all_final = zeros(file_count, 1);

for f = 1:file_count
    doy = day_of_year_all(f);
    
    ui_data = all_data{f};
    ui_lat_band = ui_data(:, lat_idx, :);
    
    % Average over time and space for this day
    daily_observed_mean(f) = mean(ui_lat_band, 'all', 'omitnan');
    
    % Get daily climatology for this day
    daily_clim_value = daily_climatology(doy);
    
    % Daily anomaly = observed - daily climatology
    daily_anom_from_daily_clim(f) = daily_observed_mean(f) - daily_clim_value;
    
    % Date
    date_str = dates_all{f};
    dates_datetime(f) = datenum(date_str, 'yyyymmdd');
    month_all_final(f) = str2double(date_str(5:6));
    
    if mod(f, 500) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

dates_datetime = datetime(dates_datetime, 'ConvertFrom', 'datenum');

fprintf('✓ Daily statistics and anomalies calculated (lat-band averaged)\n\n');

% ========================================================================
% CALCULATE GRIDPOINT-LEVEL OBSERVED DATA & ANOMALIES
% ========================================================================

fprintf('Calculating gridpoint-level observed data and anomalies...\n');

% Initialize gridpoint observed data [file_count, n_lon, n_lat_band]
gridpoint_observed = zeros(file_count, n_lon, n_lat_band);

% Initialize gridpoint anomalies [file_count, n_lon, n_lat_band]
gridpoint_anom = zeros(file_count, n_lon, n_lat_band);

for f = 1:file_count
    doy = day_of_year_all(f);
    
    ui_data = all_data{f};
    ui_lat_band = ui_data(:, lat_idx, :);  % [n_lon, n_lat_band, n_time]
    
    % Average over time for each gridpoint
    ui_gridpoint_mean = mean(ui_lat_band, 3, 'omitnan');  % [n_lon, n_lat_band]
    
    % Store the observed data
    gridpoint_observed(f, :, :) = ui_gridpoint_mean;
    
    % Get climatology for this day
    clim_gridpoint = daily_gridpoint_clim_smooth(doy, :, :);  % [1, n_lon, n_lat_band]
    
    % Anomaly for each gridpoint
    gridpoint_anom(f, :, :) = ui_gridpoint_mean - squeeze(clim_gridpoint);
    
    if mod(f, 500) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

fprintf('✓ Gridpoint observed data and anomalies calculated\n');
fprintf('  Observed data array shape: [%d observations, %d lon, %d lat]\n', ...
    size(gridpoint_observed, 1), size(gridpoint_observed, 2), size(gridpoint_observed, 3));
fprintf('  Anomaly array shape: [%d observations, %d lon, %d lat]\n\n', ...
    size(gridpoint_anom, 1), size(gridpoint_anom, 2), size(gridpoint_anom, 3));

% ========================================================================
% CALCULATE MONTHLY STATISTICS AND ANOMALIES
% ========================================================================

fprintf('Calculating monthly statistics and anomalies...\n');

monthly_stats = struct();
monthly_stats.month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
                              'July', 'August', 'September', 'October', 'November', 'December'};
monthly_stats.observed_mean = zeros(12, 1);
monthly_stats.observed_std = zeros(12, 1);
monthly_stats.anom_mean = zeros(12, 1);
monthly_stats.anom_std = zeros(12, 1);
monthly_stats.count = zeros(12, 1);

for m = 1:12
    month_mask = month_all_final == m;
    month_observed = daily_observed_mean(month_mask);
    month_anom = daily_anom_from_daily_clim(month_mask);
    
    month_observed = month_observed(~isnan(month_observed));
    month_anom = month_anom(~isnan(month_anom));
    
    if ~isempty(month_observed)
        monthly_stats.observed_mean(m) = mean(month_observed);
        monthly_stats.observed_std(m) = std(month_observed);
        monthly_stats.anom_mean(m) = mean(month_anom);
        monthly_stats.anom_std(m) = std(month_anom);
        monthly_stats.count(m) = length(month_observed);
    else
        monthly_stats.observed_mean(m) = NaN;
        monthly_stats.observed_std(m) = NaN;
        monthly_stats.anom_mean(m) = NaN;
        monthly_stats.anom_std(m) = NaN;
        monthly_stats.count(m) = 0;
    end
end

fprintf('✓ Monthly statistics and anomalies calculated\n\n');

% ========================================================================
% CALCULATE YEARLY STATISTICS
% ========================================================================

fprintf('Calculating yearly statistics...\n');

yearly_stats = struct();
yearly_stats.years = years';
yearly_stats.observed_mean = zeros(n_years, 1);
yearly_stats.observed_std = zeros(n_years, 1);
yearly_stats.anom_mean = zeros(n_years, 1);
yearly_stats.anom_std = zeros(n_years, 1);

for y = 1:n_years
    current_year = years(y);
    year_mask = year(dates_datetime) == current_year;
    year_observed = daily_observed_mean(year_mask);
    year_anom = daily_anom_from_daily_clim(year_mask);
    
    year_observed = year_observed(~isnan(year_observed));
    year_anom = year_anom(~isnan(year_anom));
    
    if ~isempty(year_observed)
        yearly_stats.observed_mean(y) = mean(year_observed);
        yearly_stats.observed_std(y) = std(year_observed);
        yearly_stats.anom_mean(y) = mean(year_anom);
        yearly_stats.anom_std(y) = std(year_anom);
    else
        yearly_stats.observed_mean(y) = NaN;
        yearly_stats.observed_std(y) = NaN;
        yearly_stats.anom_mean(y) = NaN;
        yearly_stats.anom_std(y) = NaN;
    end
end

fprintf('✓ Yearly statistics calculated\n\n');

% ========================================================================
% CALCULATE OVERALL STATISTICS
% ========================================================================

fprintf('Calculating overall statistics...\n');

overall_stats = struct();
daily_observed_valid = daily_observed_mean(~isnan(daily_observed_mean));
daily_anom_valid = daily_anom_from_daily_clim(~isnan(daily_anom_from_daily_clim));

overall_stats.daily_observed_mean = mean(daily_observed_valid);
overall_stats.daily_observed_std = std(daily_observed_valid);
overall_stats.daily_observed_min = min(daily_observed_valid);
overall_stats.daily_observed_max = max(daily_observed_valid);
overall_stats.daily_observed_median = median(daily_observed_valid);

overall_stats.daily_anom_mean = mean(daily_anom_valid);
overall_stats.daily_anom_std = std(daily_anom_valid);
overall_stats.daily_anom_min = min(daily_anom_valid);
overall_stats.daily_anom_max = max(daily_anom_valid);

overall_stats.n_daily_observations = length(daily_observed_valid);

fprintf('✓ Overall statistics calculated\n\n');

% ========================================================================
% PRINT STATISTICS SUMMARY
% ========================================================================

fprintf('========================================================================\n');
fprintf('STATISTICS SUMMARY\n');
fprintf('========================================================================\n\n');

fprintf('DAILY CLIMATOLOGY (EMPIRICAL - from real daily observations):\n');
fprintf('  Mean: %.4f m²/s\n', mean(daily_climatology_raw(~isnan(daily_climatology_raw))));
fprintf('  Std Dev: %.4f m²/s\n', std(daily_climatology_raw(~isnan(daily_climatology_raw))));
fprintf('  Min: %.4f m²/s\n', min(daily_climatology_raw(~isnan(daily_climatology_raw))));
fprintf('  Max: %.4f m²/s\n\n', max(daily_climatology_raw(~isnan(daily_climatology_raw))));

fprintf('DAILY CLIMATOLOGY (SMOOTHED - for plotting):\n');
fprintf('  Mean: %.4f m²/s\n', mean(daily_climatology(~isnan(daily_climatology))));
fprintf('  Std Dev: %.4f m²/s\n', std(daily_climatology(~isnan(daily_climatology))));
fprintf('  Min: %.4f m²/s\n', min(daily_climatology(~isnan(daily_climatology))));
fprintf('  Max: %.4f m²/s\n\n', max(daily_climatology(~isnan(daily_climatology))));

fprintf('DAILY OBSERVED DATA:\n');
fprintf('  Mean: %.4f m²/s\n', overall_stats.daily_observed_mean);
fprintf('  Std Dev: %.4f m²/s\n', overall_stats.daily_observed_std);
fprintf('  Min: %.4f m²/s\n', overall_stats.daily_observed_min);
fprintf('  Max: %.4f m²/s\n\n', overall_stats.daily_observed_max);

fprintf('DAILY ANOMALIES:\n');
fprintf('  Mean: %.4f m²/s\n', overall_stats.daily_anom_mean);
fprintf('  Std Dev: %.4f m²/s\n', overall_stats.daily_anom_std);
fprintf('  Min: %.4f m²/s\n', overall_stats.daily_anom_min);
fprintf('  Max: %.4f m²/s\n\n', overall_stats.daily_anom_max);

fprintf('MONTHLY CLIMATOLOGY:\n');
fprintf('  Mean: %.4f m²/s\n', mean(monthly_climatology));
fprintf('  Std Dev: %.4f m²/s\n', std(monthly_climatology));
fprintf('  Min: %.4f m²/s\n', min(monthly_climatology));
fprintf('  Max: %.4f m²/s\n\n', max(monthly_climatology));

fprintf('MONTHLY ANOMALIES:\n');
fprintf('  Mean: %.4f m²/s\n', mean(monthly_stats.anom_mean));
fprintf('  Std Dev: %.4f m²/s\n', std(monthly_stats.anom_mean));
fprintf('  Min: %.4f m²/s\n', min(monthly_stats.anom_mean));
fprintf('  Max: %.4f m²/s\n\n', max(monthly_stats.anom_mean));

fprintf('Monthly Statistics Table:\n');
fprintf('%-12s %10s %10s | %10s %10s | %6s\n', ...
    'Month', 'Obs Mean', 'Obs Std', 'Anom Mean', 'Anom Std', 'Count');
fprintf(repmat('-', 1, 80));
fprintf('\n');

for m = 1:12
    fprintf('%-12s %10.4f %10.4f | %10.4f %10.4f | %6d\n', ...
        monthly_stats.month_names{m}, ...
        monthly_stats.observed_mean(m), monthly_stats.observed_std(m), ...
        monthly_stats.anom_mean(m), monthly_stats.anom_std(m), monthly_stats.count(m));
end

fprintf('\n========================================================================\n\n');

% ========================================================================
% PLOTS (1-6) - NO GRID
% ========================================================================

fprintf('Creating Plot 1: Daily Climatology with Uncertainty Envelope...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

daily_clim_upper = daily_climatology + daily_climatology_std;
daily_clim_lower = daily_climatology - daily_climatology_std;

hold off;

fill([1:365, 365:-1:1], [daily_clim_upper', flipud(daily_clim_lower')], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'DisplayName', '±1 Std Dev');
hold on;

plot(1:365, daily_climatology, '-', 'LineWidth', 3, 'Color', [1, 0, 0], 'DisplayName', 'Daily climatology (smoothed)');

month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
plot(month_starts, monthly_climatology, 'o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0, 0.5, 1], 'DisplayName', 'Monthly climatology points');

hold off;

xlabel('Day of Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Daily Climatology with Variability Range - Latitude Band %.1f - %.1f°N (%d-%d)', ...
    lat_min, lat_max, years(1), years(end)), ...
    'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'south', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'black');

set(gca, 'FontSize', 11, 'LineWidth', 1.5);
set(gca, 'XLim', [1 365], 'YLim', [min(daily_clim_lower)-0.1 max(daily_clim_upper)+0.1]);

savefig(fullfile(plot_path, '01_daily_climatology_envelope.fig'));
print(fullfile(plot_path, '01_daily_climatology_envelope.png'), '-dpng', '-r300');
fprintf('✓ Plot 1 saved (300 dpi)\n');
close all;

fprintf('Creating Plot 2: Daily Anomalies Time Series...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

plot(dates_datetime, daily_anom_from_daily_clim, '-', 'LineWidth', 0.3, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Daily anomalies');
hold on;
plot(dates_datetime, movmean(daily_anom_from_daily_clim, 10, 'omitnan'), '-', 'LineWidth', 2, 'Color', [1, 0, 0], 'DisplayName', '10-day moving mean');
plot(dates_datetime([1, end]), [0, 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'Zero anomaly', 'HandleVisibility', 'off');
hold off;

xlabel('Date', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Daily Anomalies Time Series - Latitude Band %.1f - %.1f°N', lat_min, lat_max), ...
    'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'south', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'black');

set(gca, 'FontSize', 11, 'LineWidth', 1.5);

savefig(fullfile(plot_path, '02_daily_anomalies_timeseries.fig'));
print(fullfile(plot_path, '02_daily_anomalies_timeseries.png'), '-dpng', '-r300');
fprintf('✓ Plot 2 saved (300 dpi)\n');
close all;

fprintf('Creating Plot 3: Monthly Climatology with Error Bars...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

errorbar(1:12, monthly_climatology, monthly_climatology_std, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, 'Color', [0, 0.5, 1]);

xlabel('Month', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Monthly Climatology ± Std Dev - Latitude Band %.1f - %.1f°N (%d-%d)', ...
    lat_min, lat_max, years(1), years(end)), ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'XTick', 1:12);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

savefig(fullfile(plot_path, '03_monthly_climatology.fig'));
print(fullfile(plot_path, '03_monthly_climatology.png'), '-dpng', '-r300');
fprintf('✓ Plot 3 saved (300 dpi)\n');
close all;

fprintf('Creating Plot 4: Monthly Anomalies with Error Bars...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

errorbar(1:12, monthly_stats.anom_mean, monthly_stats.anom_std, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, 'Color', [0, 0.7, 0]);
hold on;
plot([0 13], [0 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;

xlabel('Month', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Monthly Mean Anomaly ± Std Dev - Latitude Band %.1f - %.1f°N (%d-%d)', ...
    lat_min, lat_max, years(1), years(end)), ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'XTick', 1:12);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

savefig(fullfile(plot_path, '04_monthly_anomalies.fig'));
print(fullfile(plot_path, '04_monthly_anomalies.png'), '-dpng', '-r300');
fprintf('✓ Plot 4 saved (300 dpi)\n');
close all;

fprintf('Creating Plot 5: Daily vs Monthly Climatology Comparison...\n\n');

figure('Position', [100, 100, 1400, 650], 'Color', 'W');

subplot(2, 1, 1);
hold off;
fill([1:365, 365:-1:1], [daily_clim_upper', flipud(daily_clim_lower')], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'DisplayName', '±1 Std Dev (daily)');
hold on;
plot(1:365, daily_climatology, '-', 'LineWidth', 2.5, 'Color', [1, 0, 0], 'DisplayName', 'Daily climatology (smoothed)');
plot(month_starts, monthly_climatology, 'o-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', [0, 0.5, 1], 'DisplayName', 'Monthly climatology (reference)');
hold off;
xlabel('Day of Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Daily vs Monthly Climatology', 'FontSize', 13, 'FontWeight', 'bold');

legend('Location', 'south', 'FontSize', 11);
set(gca, 'XLim', [1 365]);

subplot(2, 1, 2);
daily_anom_movmean = movmean(daily_anom_from_daily_clim, 10, 'omitnan');
plot(dates_datetime, daily_anom_movmean, '-', 'LineWidth', 1.5, 'Color', [0.7, 0.7, 0.7], 'DisplayName', '10-day moving mean');
hold on;
plot(dates_datetime([1, end]), [0, 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
hold off;
xlabel('Date', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Daily Anomalies (10-day smoothed)', 'FontSize', 13, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 11);

sgtitle(sprintf('Climatology Comparison - Latitude Band %.1f - %.1f°N', lat_min, lat_max), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '05_climatology_comparison.fig'));
print(fullfile(plot_path, '05_climatology_comparison.png'), '-dpng', '-r300');
fprintf('✓ Plot 5 saved (300 dpi)\n');
close all;

fprintf('Creating Plot 6: Yearly Trends and Statistics...\n\n');

figure('Position', [100, 100, 1400, 750], 'Color', 'W');

% Subplot 1: Yearly mean observed
subplot(2, 2, 1);
errorbar(yearly_stats.years, yearly_stats.observed_mean, yearly_stats.observed_std, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [1, 0, 0]);
xlabel('Year', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Yearly Mean ± Std Dev (Observed)', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 2: Yearly anomalies
subplot(2, 2, 2);
errorbar(yearly_stats.years, yearly_stats.anom_mean, yearly_stats.anom_std, 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'Color', [0, 0.7, 0]);
hold on;
plot([years(1)-1 years(end)+1], [0 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1);
hold off;
xlabel('Year', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Yearly Mean Anomaly ± Std Dev', 'FontSize', 12, 'FontWeight', 'bold');

% Subplot 3: Monthly comparison
subplot(2, 2, 3);
bar(1:12, monthly_stats.observed_mean, 'FaceColor', [0.5, 0.5, 1]);
xlabel('Month', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Monthly Observed Mean', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:12);

% Subplot 4: Monthly anomalies comparison
subplot(2, 2, 4);
bar(1:12, monthly_stats.anom_mean, 'FaceColor', [0.5, 1, 0.5]);
hold on;
plot([0 13], [0 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1);
hold off;
xlabel('Month', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Monthly Mean Anomaly', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:12);

sgtitle(sprintf('Yearly & Monthly Trends - Latitude Band %.1f - %.1f°N', lat_min, lat_max), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '06_yearly_monthly_trends.fig'));
print(fullfile(plot_path, '06_yearly_monthly_trends.png'), '-dpng', '-r300');
fprintf('✓ Plot 6 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 7: SPATIAL ANOMALY MAPS (Sample days)
% ========================================================================

fprintf('Creating Plot 7: Spatial Anomaly Maps...\n\n');

% Select 4 representative days throughout the year
sample_days = [1, 100, 200, 300];  % Jan 1, Apr 10, Jul 19, Oct 27

figure('Position', [100, 100, 1400, 900], 'Color', 'W');

for i = 1:length(sample_days)
    doy = sample_days(i);
    
    % Calculate mean anomaly for this day across all years
    anom_idx = (day_of_year_all == doy);
    if sum(anom_idx) > 0
        anom_map = squeeze(mean(gridpoint_anom(anom_idx, :, :), 1, 'omitnan'));  % [n_lon, n_lat_band]
        
        subplot(2, 2, i);
        
        % Get longitude and latitude for the latitude band
        lon_grid = lon;
        lat_grid = lat(lat_idx);
        
        % Create meshgrid
        [LON, LAT] = meshgrid(lon_grid, lat_grid);
        LON = LON';
        LAT = LAT';
        
        % Plot contourf with symmetric colormap
        anom_max = max(abs(anom_map(:)));
        contourf(LON, LAT, anom_map, 20, 'LineStyle', 'none');
        caxis([-anom_max anom_max]);
        colormap(gca, flipud(jet));
        colorbar('FontSize', 10);
        
        hold on;
        contour(LON, LAT, anom_map, [0 0], 'k-', 'LineWidth', 1.5);
        hold off;
        
        xlabel('Longitude (°E)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Latitude (°N)', 'FontSize', 11, 'FontWeight', 'bold');
        
        % Calculate month and day from DOY
        days_in_months = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
        month_idx = find(doy >= days_in_months, 1, 'last');
        day_in_month = doy - days_in_months(month_idx);
        month_names = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
        
        title(sprintf('DOY %d (%s %d) - Mean Anomaly', doy, month_names{month_idx}, day_in_month), ...
            'FontSize', 12, 'FontWeight', 'bold');
    end
end

sgtitle(sprintf('Spatial Anomaly Maps - Latitude Band %.1f - %.1f°N', lat_min, lat_max), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '07_spatial_anomaly_maps.fig'));
print(fullfile(plot_path, '07_spatial_anomaly_maps.png'), '-dpng', '-r300');
fprintf('✓ Plot 7 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 8: SEASONAL ANOMALY MAPS (Winter, Spring, Summer, Fall)
% ========================================================================

fprintf('Creating Plot 8: Seasonal Anomaly Maps...\n\n');

% Define seasons
seasons = {
    'Winter', [355:365, 1:59];           % Dec 21 - Mar 20
    'Spring', 60:151;                    % Mar 21 - Jun 20
    'Summer', 152:243;                   % Jun 21 - Sep 22
    'Fall', 244:334;                     % Sep 23 - Dec 20
};

figure('Position', [100, 100, 1400, 900], 'Color', 'W');

for i = 1:4
    season_name = seasons{i, 1};
    season_doys = seasons{i, 2};
    
    % Get longitude and latitude for the latitude band
    lon_grid = lon;
    lat_grid = lat(lat_idx);
    
    % Create meshgrid
    [LON, LAT] = meshgrid(lon_grid, lat_grid);
    LON = LON';
    LAT = LAT';
    
    % Calculate mean anomaly for all observations in this season
    season_anom_map = zeros(n_lon, n_lat_band);
    count_map = zeros(n_lon, n_lat_band);
    
    for doy = season_doys
        anom_idx = (day_of_year_all == doy);
        if sum(anom_idx) > 0
            anom_data = gridpoint_anom(anom_idx, :, :);
            season_anom_map = season_anom_map + squeeze(sum(anom_data, 1, 'omitnan'));
            count_map = count_map + squeeze(sum(~isnan(anom_data), 1, 'omitnan'));
        end
    end
    
    % Average
    season_anom_map = season_anom_map ./ (count_map + (count_map == 0));
    
    subplot(2, 2, i);
    
    % Plot contourf with symmetric colormap
    anom_max = max(abs(season_anom_map(:)));
    contourf(LON, LAT, season_anom_map, 20, 'LineStyle', 'none');
    caxis([-anom_max anom_max]);
    colormap(gca, flipud(jet));
    colorbar('FontSize', 10);
    
    hold on;
    contour(LON, LAT, season_anom_map, [0 0], 'k-', 'LineWidth', 1.5);
    hold off;
    
    xlabel('Longitude (°E)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Latitude (°N)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s - Mean Anomaly', season_name), 'FontSize', 12, 'FontWeight', 'bold');
end

sgtitle(sprintf('Seasonal Anomaly Maps - Latitude Band %.1f - %.1f°N', lat_min, lat_max), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '08_seasonal_anomaly_maps.fig'));
print(fullfile(plot_path, '08_seasonal_anomaly_maps.png'), '-dpng', '-r300');
fprintf('✓ Plot 8 saved (300 dpi)\n');
close all;

% ========================================================================
% SAVE RESULTS
% ========================================================================

fprintf('\nSaving comprehensive results...\n');

save(fullfile(output_path, 'true_daily_analysis_lat25to29.mat'), ...
    'daily_climatology_raw', 'daily_climatology_raw_std', ...
    'daily_climatology', 'daily_climatology_std', ...
    'daily_obs_count', 'daily_years_count', ...
    'daily_gridpoint_clim', 'daily_gridpoint_clim_smooth', ...
    'daily_gridpoint_std', 'daily_gridpoint_std_smooth', ...
    'daily_gridpoint_count', ...
    'gridpoint_observed', ...
    'gridpoint_anom', ...
    'monthly_climatology', 'monthly_climatology_std', ...
    'daily_observed_mean', 'daily_anom_from_daily_clim', ...
    'monthly_stats', 'yearly_stats', 'overall_stats', ...
    'dates_datetime', 'dates_all', 'day_of_year_all', 'month_all_final', ...
    'lat_min', 'lat_max', 'lon', 'lat', 'years', '-v7.3');

fprintf('✓ Results saved to: %s\n\n', fullfile(output_path, 'true_daily_analysis_lat25to29.mat'));

fprintf('✓ TRUE DAILY CLIMATOLOGY ANALYSIS COMPLETE\n');
