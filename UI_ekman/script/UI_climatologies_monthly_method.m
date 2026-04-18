% ========================================================================
% SCRIPT: COMPREHENSIVE CLIMATOLOGY & ANOMALIES ANALYSIS WITH GRIDPOINT DATA
% Latitude band: 25-29°N, Longitude: -20 to -10°W
% Seasonal grouping with shared colorbars
% ========================================================================
clear all; close all; clc;

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman';
output_path = fullfile(ui_results_path, 'analysis_monthly_climatologies');
plot_path = fullfile(output_path, 'plots');

% Create output directories
if ~isfolder(output_path)
    mkdir(output_path);
end
if ~isfolder(plot_path)
    mkdir(plot_path);
end

% PARAMETERS
mask_name = '50km';
years = 2005:2024;
n_years = length(years);
lat_min = 25;
lat_max = 29;
lon_min = -20;
lon_max = -10;

fprintf('Analysis Parameters:\n');
fprintf('  Latitude band: %.1f - %.1f°N\n', lat_min, lat_max);
fprintf('  Longitude range: %.1f - %.1f°W\n', abs(lon_max), abs(lon_min));
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
    
    if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
        days_in_year = 366;
    else
        days_in_year = 365;
    end
    
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
            current_day = d;
            
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

fprintf(' Loaded %d files\n', file_count);
fprintf('  Grid size: lon=%d, lat=%d\n', length(lon), length(lat));

[n_lon, n_lat, n_time] = size(UI_data);
fprintf('  Time steps per day: %d\n\n', n_time);

% ========================================================================
% FIND LATITUDE AND LONGITUDE INDICES FOR REGION
% ========================================================================

fprintf('Finding latitude band %.1f - %.1f°N and longitude %.1f to %.1f°W...\n', lat_min, lat_max, abs(lon_max), abs(lon_min));

lat_idx = find(lat >= lat_min & lat <= lat_max);
lon_idx = find(lon >= lon_min & lon <= lon_max);

if isempty(lat_idx)
    error('No data found in latitude band %.1f - %.1f°N', lat_min, lat_max);
end
if isempty(lon_idx)
    error('No data found in longitude range %.1f - %.1f', lon_min, lon_max);
end

n_lat_band = length(lat_idx);
n_lon_region = length(lon_idx);

fprintf('✓ Found %d latitude points in band\n', n_lat_band);
fprintf('✓ Found %d longitude points in region\n', n_lon_region);
fprintf('  Latitude indices: %d to %d\n', lat_idx(1), lat_idx(end));
fprintf('  Longitude indices: %d to %d\n\n', lon_idx(1), lon_idx(end));

% ========================================================================
% CALCULATE CORRECTED MONTHLY CLIMATOLOGY (FOR LAT BAND AND LON REGION)
% ========================================================================

fprintf('Calculating monthly climatology for region...\n');

% Get month for each file
month_all = zeros(file_count, 1);
year_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    month_all(f) = month_val;
    year_all(f) = str2double(date_str(1:4));
end

% Initialize: store ALL observations for each month
monthly_data_pool = cell(12, 1);
monthly_gridpoint_pool = cell(12, n_lon_region, n_lat_band);

for m = 1:12
    monthly_data_pool{m} = [];
    for ilon = 1:n_lon_region
        for ilat = 1:n_lat_band
            monthly_gridpoint_pool{m, ilon, ilat} = [];
        end
    end
end

% Collect observations from region
fprintf('  Pooling observations by month...\n');
for f = 1:file_count
    m = month_all(f);
    ui_data = all_data{f};
    
    % Extract region: [lon_region, lat_band, time]
    ui_region = ui_data(lon_idx, lat_idx, :);
    
    % Reshape to [lon_region*lat_band, time] and append (for region average)
    ui_reshaped = reshape(ui_region, [], size(ui_region, 3));
    monthly_data_pool{m} = [monthly_data_pool{m}, ui_reshaped];
    
    % Store gridpoint-level data
    for ilon = 1:n_lon_region
        for ilat = 1:n_lat_band
            gridpoint_values = squeeze(ui_region(ilon, ilat, :));
            monthly_gridpoint_pool{m, ilon, ilat} = [monthly_gridpoint_pool{m, ilon, ilat}; gridpoint_values'];
        end
    end
    
    if mod(f, 500) == 0
        fprintf('    Processed file %d/%d\n', f, file_count);
    end
end

% Initialize: [12 months]
monthly_climatology = zeros(12, 1);
monthly_climatology_std = zeros(12, 1);

% Initialize gridpoint-level monthly climatology [12 months, n_lon_region, n_lat_band]
monthly_gridpoint_clim = zeros(12, n_lon_region, n_lat_band);
monthly_gridpoint_std = zeros(12, n_lon_region, n_lat_band);

% Compute monthly climatology
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
    
    % Gridpoint-level monthly climatology
    for ilon = 1:n_lon_region
        for ilat = 1:n_lat_band
            gridpoint_obs = monthly_gridpoint_pool{m, ilon, ilat};
            gridpoint_obs = gridpoint_obs(~isnan(gridpoint_obs));
            
            if ~isempty(gridpoint_obs)
                monthly_gridpoint_clim(m, ilon, ilat) = mean(gridpoint_obs, 'omitnan');
                monthly_gridpoint_std(m, ilon, ilat) = std(gridpoint_obs, 'omitnan');
            else
                monthly_gridpoint_clim(m, ilon, ilat) = NaN;
                monthly_gridpoint_std(m, ilon, ilat) = NaN;
            end
        end
    end
end

fprintf(' Monthly climatology calculated\n');
fprintf('  Gridpoint climatology shape: [12 months, %d lon, %d lat]\n\n', n_lon_region, n_lat_band);

clear monthly_data_pool monthly_gridpoint_pool;

% ========================================================================
% CREATE SMOOTH DAILY CLIMATOLOGY WITH PADDING
% Linear Interpolation + Gaussian Smoothing
% ========================================================================

% Days per month
days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
half_month = 15;

% Calculate midpoint of each month (in the 365-day year)
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
extended_midpoints = zeros(14, 1);
extended_values_mean = zeros(14, 1);
extended_values_std = zeros(14, 1);

extended_midpoints(1) = 15;
extended_values_mean(1) = monthly_climatology(12);
extended_values_std(1) = monthly_climatology_std(12);

for m = 1:12
    extended_midpoints(m + 1) = 15 + month_midpoints_365(m);
    extended_values_mean(m + 1) = monthly_climatology(m);
    extended_values_std(m + 1) = monthly_climatology_std(m);
end

extended_midpoints(14) = 390;
extended_values_mean(14) = monthly_climatology(1);
extended_values_std(14) = monthly_climatology_std(1);

fprintf('\n  Extended calendar midpoints (390-day period):\n');
for i = 1:14
    fprintf('    Point %2d: day %6.1f, mean: %.4f, std: %.4f m²/s\n', i, extended_midpoints(i), extended_values_mean(i), extended_values_std(i));
end

% ========================================================================
% SMOOTH INTERPOLATION ACROSS DAYS
% ========================================================================

fprintf('\nSmoothing climatology with linear interpolation + Gaussian smoothing...\n');

climatology_smooth_mean = zeros(390, 1);
climatology_smooth_std = zeros(390, 1);

for doy_extended = 1:390
    valid_idx = ~isnan(extended_values_mean);
    valid_midpoints = extended_midpoints(valid_idx);
    valid_values_mean = extended_values_mean(valid_idx);
    valid_values_std = extended_values_std(valid_idx);
    
    if length(valid_idx) >= 2
        interp_val_mean = interp1(valid_midpoints, valid_values_mean, doy_extended, 'linear');
        interp_val_std = interp1(valid_midpoints, valid_values_std, doy_extended, 'linear');
        
        climatology_smooth_mean(doy_extended) = interp_val_mean;
        climatology_smooth_std(doy_extended) = interp_val_std;
    else
        climatology_smooth_mean(doy_extended) = NaN;
        climatology_smooth_std(doy_extended) = NaN;
    end
end

climatology_smooth_mean = smoothdata(climatology_smooth_mean, 'gaussian', 15);
climatology_smooth_std = smoothdata(climatology_smooth_std, 'gaussian', 15);

fprintf(' Interpolation and smoothing complete\n\n');

% ========================================================================
% EXTRACT FINAL CLIMATOLOGY (365 DAYS, REMOVING PADDING)
% ========================================================================

fprintf('Extracting final daily climatology (365 days)...\n');

daily_climatology = climatology_smooth_mean((half_month+1):(half_month+365));
daily_climatology_std = climatology_smooth_std((half_month+1):(half_month+365));

fprintf(' Final climatology extracted: %d days\n\n', length(daily_climatology));

% ========================================================================
% CALCULATE DAY OF YEAR FOR EACH FILE
% ========================================================================

fprintf('Calculating day of year for each file...\n');

day_of_year_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    day_val = str2double(date_str(7:8));
    
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    day_of_year = sum(days_in_months(1:month_val-1)) + day_val;
    day_of_year_all(f) = day_of_year;
end

fprintf(' Day of year calculated\n\n');

% ========================================================================
% CALCULATE DAILY ANOMALIES AND STATISTICS
% ========================================================================

fprintf('Calculating daily statistics and anomalies...\n');

daily_observed_mean = zeros(file_count, 1);
daily_anom_from_daily_clim = zeros(file_count, 1);
dates_datetime = zeros(file_count, 1);
month_all_final = zeros(file_count, 1);

for f = 1:file_count
    doy = day_of_year_all(f);
    
    ui_data = all_data{f};
    ui_region = ui_data(lon_idx, lat_idx, :);
    
    daily_observed_mean(f) = mean(ui_region, 'all', 'omitnan');
    daily_clim_value = daily_climatology(doy);
    daily_anom_from_daily_clim(f) = daily_observed_mean(f) - daily_clim_value;
    
    date_str = dates_all{f};
    dates_datetime(f) = datenum(date_str, 'yyyymmdd');
    month_all_final(f) = str2double(date_str(5:6));
    
    if mod(f, 500) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

dates_datetime = datetime(dates_datetime, 'ConvertFrom', 'datenum');

fprintf(' Daily statistics and anomalies calculated\n\n');

% ========================================================================
% CALCULATE GRIDPOINT-LEVEL OBSERVED DATA & ANOMALIES
% ========================================================================

fprintf('Calculating gridpoint-level observed data and anomalies...\n');

gridpoint_observed = zeros(file_count, n_lon_region, n_lat_band);
gridpoint_climatology = zeros(file_count, n_lon_region, n_lat_band);
gridpoint_anom = zeros(file_count, n_lon_region, n_lat_band);

for f = 1:file_count
    m = month_all_final(f);
    
    ui_data = all_data{f};
    ui_region = ui_data(lon_idx, lat_idx, :);
    
    ui_gridpoint_mean = mean(ui_region, 3, 'omitnan');
    gridpoint_observed(f, :, :) = ui_gridpoint_mean;
    
    gridpoint_clim_month = squeeze(monthly_gridpoint_clim(m, :, :));
    gridpoint_climatology(f, :, :) = gridpoint_clim_month;
    gridpoint_anom(f, :, :) = ui_gridpoint_mean - gridpoint_clim_month;
    
    if mod(f, 500) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

fprintf('✓ Gridpoint observed data and anomalies calculated\n\n');

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

fprintf(' Monthly statistics and anomalies calculated\n\n');

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

fprintf(' Yearly statistics calculated\n\n');

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

fprintf(' Overall statistics calculated\n\n');

% Display statistics

fprintf('STATISTICS SUMMARY\n');

fprintf('DAILY CLIMATOLOGY:\n');
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
% PLOT 1: DAILY CLIMATOLOGY WITH UNCERTAINTY ENVELOPE
% ========================================================================

fprintf('Creating Plot 1: Daily Climatology with Uncertainty Envelope...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

daily_clim_upper = daily_climatology + daily_climatology_std;
daily_clim_lower = daily_climatology - daily_climatology_std;

hold off;

fill([1:365, 365:-1:1], [daily_clim_upper', flipud(daily_clim_lower')], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'DisplayName', '±1 Std Dev');
hold on;

plot(1:365, daily_climatology, '-', 'LineWidth', 3, 'Color', [1, 0, 0], 'DisplayName', 'Daily climatology');

month_starts = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
plot(month_starts, monthly_climatology, 'o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [0, 0.5, 1], 'DisplayName', 'Monthly climatology points');

hold off;

xlabel('Day of Year', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Daily Climatology with Variability Range - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'south', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'black');

set(gca, 'FontSize', 11, 'LineWidth', 1.5);
set(gca, 'XLim', [1 365], 'YLim', [min(daily_clim_lower)-0.1 max(daily_clim_upper)+0.1]);

savefig(fullfile(plot_path, '01_daily_climatology_envelope.fig'));
print(fullfile(plot_path, '01_daily_climatology_envelope.png'), '-dpng', '-r300');
fprintf(' Plot 1 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 2: DAILY ANOMALIES TIME SERIES
% ========================================================================

fprintf('Creating Plot 2: Daily Anomalies Time Series...\n\n');

figure('Position', [100, 100, 1000, 600],'Color', 'W');

plot(dates_datetime, daily_anom_from_daily_clim, '-', 'LineWidth', 0.3, 'Color', [0.7, 0.7, 0.7], 'DisplayName', 'Daily anomalies');
hold on;
plot(dates_datetime, movmean(daily_anom_from_daily_clim, 10, 'omitnan'), '-', 'LineWidth', 2, 'Color', [1, 0, 0], 'DisplayName', '10-day moving mean');
plot(dates_datetime([1, end]), [0, 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'Zero anomaly', 'HandleVisibility', 'off');
hold off;

xlabel('Date', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Daily Anomalies Time Series - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

legend('Location', 'south', 'FontSize', 12, 'Box', 'on', 'EdgeColor', 'black');

set(gca, 'FontSize', 11, 'LineWidth', 1.5);

savefig(fullfile(plot_path, '02_daily_anomalies_timeseries.fig'));
print(fullfile(plot_path, '02_daily_anomalies_timeseries.png'), '-dpng', '-r300');
fprintf(' Plot 2 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 3: MONTHLY CLIMATOLOGY WITH ERROR BARS
% ========================================================================

fprintf('Creating Plot 3: Monthly Climatology with Error Bars...\n\n');

figure('Position', [100, 100, 1000, 600],'Color', 'W');

errorbar(1:12, monthly_climatology, monthly_climatology_std, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, 'Color', [0, 0.5, 1]);

xlabel('Month', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Monthly Climatology ± Std Dev - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'XTick', 1:12);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

savefig(fullfile(plot_path, '03_monthly_climatology.fig'));
print(fullfile(plot_path, '03_monthly_climatology.png'), '-dpng', '-r300');
fprintf(' Plot 3 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 4: MONTHLY ANOMALIES WITH ERROR BARS
% ========================================================================

fprintf('Creating Plot 4: Monthly Anomalies with Error Bars...\n\n');

figure('Position', [100, 100, 1000, 600], 'Color', 'W');

errorbar(1:12, monthly_stats.anom_mean, monthly_stats.anom_std, 'o-', 'LineWidth', 2.5, 'MarkerSize', 10, 'Color', [0, 0.7, 0]);
hold on;
plot([0 13], [0 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1.5, 'DisplayName', 'Zero anomaly', 'HandleVisibility', 'off');
hold off;

xlabel('Month', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Monthly Mean Anomaly ± Std Dev - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

set(gca, 'XTick', 1:12);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

legend('Location', 'best', 'FontSize', 12);

savefig(fullfile(plot_path, '04_monthly_anomalies.fig'));
print(fullfile(plot_path, '04_monthly_anomalies.png'), '-dpng', '-r300');
fprintf(' Plot 4 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 5: DAILY VS MONTHLY CLIMATOLOGY COMPARISON
% ========================================================================

fprintf('Creating Plot 5: Daily vs Monthly Climatology Comparison...\n\n');

figure('Position', [100, 100, 1400, 650], 'Color', 'W');

subplot(2, 1, 1);
hold off;
fill([1:365, 365:-1:1], [daily_clim_upper', flipud(daily_clim_lower')], ...
    [0.7, 0.85, 1], 'EdgeColor', 'none', 'DisplayName', '±1 Std Dev (daily)');
hold on;
plot(1:365, daily_climatology, '-', 'LineWidth', 2.5, 'Color', [1, 0, 0], 'DisplayName', 'Daily climatology (interpolated & smoothed)');
plot(month_starts, monthly_climatology, 'o-', 'LineWidth', 2, 'MarkerSize', 10, 'Color', [0, 0.5, 1], 'DisplayName', 'Monthly climatology (reference points)');
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
plot([1, length(dates_datetime)], [mean(monthly_stats.anom_mean), mean(monthly_stats.anom_mean)], '--', 'Color', [1, 0, 0], 'LineWidth', 2, 'DisplayName', 'Mean monthly anomaly', 'HandleVisibility', 'off');
hold off;
xlabel('Date', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Daily Anomalies (10-day smoothed)', 'FontSize', 13, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 11);

sgtitle(sprintf('Climatology Comparison - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '05_climatology_comparison.fig'));
print(fullfile(plot_path, '05_climatology_comparison.png'), '-dpng', '-r300');
fprintf(' Plot 5 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 6: YEARLY TRENDS AND STATISTICS
% ========================================================================

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
bar(1:12, monthly_stats.observed_mean, 'FaceColor', [0.5, 0.5, 1], 'DisplayName', 'Monthly observed mean');
xlabel('Month', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('UI Ekman (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Monthly Observed Mean', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:12);

% Subplot 4: Monthly anomalies comparison
subplot(2, 2, 4);
bar(1:12, monthly_stats.anom_mean, 'FaceColor', [0.5, 1, 0.5], 'DisplayName', 'Monthly anomaly mean');
hold on;
plot([0 13], [0 0], '--', 'Color', [0, 0, 0], 'LineWidth', 1);
hold off;
xlabel('Month', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Anomaly (m²/s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Monthly Mean Anomaly', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTick', 1:12);

sgtitle(sprintf('Yearly & Monthly Trends - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '06_yearly_monthly_trends.fig'));
print(fullfile(plot_path, '06_yearly_monthly_trends.png'), '-dpng', '-r300');
fprintf(' Plot 6 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 7: SEASONAL CLIMATOLOGY MAPS (4 Seasons - Shared Colorbar)
% ========================================================================

fprintf('Creating Plot 7: Seasonal Climatology Maps (Shared Colorbar)...\n\n');

% Define seasons and months
seasons = struct();
seasons(1).name = 'Winter';
seasons(1).months = [12, 1, 2];
seasons(2).name = 'Spring';
seasons(2).months = [3, 4, 5];
seasons(3).name = 'Summer';
seasons(3).months = [6, 7, 8];
seasons(4).name = 'Fall';
seasons(4).months = [9, 10, 11];

% Calculate seasonal climatology maps
seasonal_clim = zeros(4, n_lon_region, n_lat_band);

for s = 1:4
    clim_map = zeros(n_lon_region, n_lat_band);
    count_map = zeros(n_lon_region, n_lat_band);
    
    for m = seasons(s).months
        month_clim = squeeze(monthly_gridpoint_clim(m, :, :));
        clim_map = clim_map + month_clim;
        count_map = count_map + ~isnan(month_clim);
    end
    
    seasonal_clim(s, :, :) = clim_map ./ (count_map + (count_map == 0));
end

% Find global min/max for shared colorbar
global_clim_min = min(seasonal_clim(~isnan(seasonal_clim)));
global_clim_max = max(seasonal_clim(~isnan(seasonal_clim)));

% Plot
figure('Position', [100, 100, 1600, 900], 'Color', 'W');

lon_region = lon(lon_idx);
lat_region = lat(lat_idx);

for s = 1:4
    subplot(2, 2, s);
    
    clim_map = squeeze(seasonal_clim(s, :, :));
    
    [LON, LAT] = meshgrid(lon_region, lat_region);
    LON = LON';
    LAT = LAT';
    
    contourf(LON, LAT, clim_map, 25, 'LineStyle', 'none');
    colormap(gca, 'jet');
    caxis([global_clim_min global_clim_max]);
    cb = colorbar('FontSize', 10);
    cb.Label.String = 'UI Climatology (m²/s)';
    
    xlabel('Longitude (°W)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Latitude (°N)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s', seasons(s).name), 'FontSize', 12, 'FontWeight', 'bold');
    
    axis tight;
    set(gca, 'FontSize', 10);
    
    xticks_val = get(gca, 'XTick');
    set(gca, 'XTickLabel', num2str(abs(xticks_val')));
end

sgtitle(sprintf('Seasonal Climatology Maps (Shared Colorbar) - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '07_seasonal_climatology_maps.fig'));
print(fullfile(plot_path, '07_seasonal_climatology_maps.png'), '-dpng', '-r300');
fprintf(' Plot 7 saved (300 dpi)\n');
close all;

% ========================================================================
% PLOT 8: SEASONAL ANOMALY MAPS (4 Seasons - Shared Colorbar)
% ========================================================================

fprintf('Creating Plot 8: Seasonal Anomaly Maps (Shared Colorbar)...\n\n');

% Calculate seasonal anomaly maps
seasonal_anom = zeros(4, n_lon_region, n_lat_band);

for s = 1:4
    anom_map = zeros(n_lon_region, n_lat_band);
    count_map = zeros(n_lon_region, n_lat_band);
    
    for m = seasons(s).months
        month_mask = (month_all_final == m);
        if sum(month_mask) > 0
            month_anom = squeeze(mean(gridpoint_anom(month_mask, :, :), 1, 'omitnan'));
            anom_map = anom_map + month_anom;
            count_map = count_map + ~isnan(month_anom);
        end
    end
    
    seasonal_anom(s, :, :) = anom_map ./ (count_map + (count_map == 0));
end

% Find global min/max for shared colorbar (symmetric)
global_anom_max = max(abs(seasonal_anom(~isnan(seasonal_anom))));

% Plot
figure('Position', [100, 100, 1600, 900], 'Color', 'W');

for s = 1:4
    subplot(2, 2, s);
    
    anom_map = squeeze(seasonal_anom(s, :, :));
    
    [LON, LAT] = meshgrid(lon_region, lat_region);
    LON = LON';
    LAT = LAT';
    
    contourf(LON, LAT, anom_map, 25, 'LineStyle', 'none');
    colormap(gca, 'jet');
    caxis([-global_anom_max global_anom_max]);
    cb = colorbar('FontSize', 10);
    cb.Label.String = 'Anomaly (m²/s)';
    
    hold on;
    contour(LON, LAT, anom_map, [0 0], 'k-', 'LineWidth', 2);
    hold off;
    
    xlabel('Longitude (°W)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('Latitude (°N)', 'FontSize', 11, 'FontWeight', 'bold');
    title(sprintf('%s', seasons(s).name), 'FontSize', 12, 'FontWeight', 'bold');
    
    axis tight;
    set(gca, 'FontSize', 10);
    
    xticks_val = get(gca, 'XTick');
    set(gca, 'XTickLabel', num2str(abs(xticks_val')));
end

sgtitle(sprintf('Seasonal Anomaly Maps (Shared Colorbar) - Region (Lat %.1f-%.1f°N, Lon %.1f-%.1f°W)', ...
    lat_min, lat_max, abs(lon_max), abs(lon_min)), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, '08_seasonal_anomaly_maps.fig'));
print(fullfile(plot_path, '08_seasonal_anomaly_maps.png'), '-dpng', '-r300');
fprintf(' Plot 8 saved (300 dpi)\n');
close all;

% ========================================================================
% SAVE RESULTS
% ========================================================================

fprintf('\nSaving comprehensive results...\n');

save(fullfile(output_path, 'comprehensive_analysis_lat25to29.mat'), ...
    'daily_climatology', 'daily_climatology_std', ...
    'monthly_climatology', 'monthly_climatology_std', ...
    'monthly_gridpoint_clim', 'monthly_gridpoint_std', ...
    'gridpoint_observed', 'gridpoint_climatology', 'gridpoint_anom', ...
    'daily_observed_mean', 'daily_anom_from_daily_clim', ...
    'monthly_stats', 'yearly_stats', 'overall_stats', ...
    'dates_datetime', 'dates_all', 'day_of_year_all', 'month_all_final', ...
    'lat_min', 'lat_max', 'lon_min', 'lon_max', 'lon', 'lat', 'years', '-v7.3');

fprintf('✓ Results saved to: %s\n\n', fullfile(output_path, 'comprehensive_analysis_lat25to29.mat'));

fprintf('========================================================================\n');
fprintf('✓ COMPREHENSIVE CLIMATOLOGY ANALYSIS COMPLETE\n');
