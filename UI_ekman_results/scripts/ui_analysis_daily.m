% ========================================================================
% SCRIPT 2: DAILY CLIMATOLOGY ANALYSIS
% Calculate daily climatology (365 days) with smoothing & statistics
% ========================================================================
clear all; close all; clc;

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman_results';
output_path = fullfile(ui_results_path, 'analysis_daily');
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

fprintf('✓ Loaded %d files\n', file_count);
fprintf('  Grid size: lon=%d, lat=%d\n', length(lon), length(lat));

[n_lon, n_lat, n_time] = size(UI_data);
fprintf('  Time steps per day: %d\n\n', n_time);

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

fprintf('✓ Day of year calculated\n\n');

% ========================================================================
% CALCULATE TRUE DAILY CLIMATOLOGY (365 DAYS)
% ========================================================================

fprintf('Calculating true daily climatology (day-of-year average)...\n');

% Initialize: store ALL observations for each day of year
daily_data_pool = cell(365, 1);
for doy = 1:365
    daily_data_pool{doy} = [];
end

% Collect all observations by day-of-year
fprintf('  Pooling all observations by day-of-year...\n');
for f = 1:file_count
    doy = day_of_year_all(f);
    ui_data = all_data{f};  % [lon, lat, time]
    
    % Reshape to [lon*lat, time] and append to pool for this DOY
    ui_reshaped = reshape(ui_data, [], size(ui_data, 3));
    daily_data_pool{doy} = [daily_data_pool{doy}, ui_reshaped];
    
    if mod(f, 200) == 0
        fprintf('    Processed file %d/%d\n', f, file_count);
    end
end

% Initialize daily climatology: [lon, lat, 365]
climatology = zeros(n_lon, n_lat, 365, 'single');
daily_count = zeros(365, 1);

% Compute true mean for each day of year from pooled observations
fprintf('  Computing daily means from pooled observations...\n');
for doy = 1:365
    all_obs = daily_data_pool{doy};  % [lon*lat, all_observations_on_this_doy]
    
    if ~isempty(all_obs)
        daily_mean = mean(all_obs, 2, 'omitnan');  % Mean over all 20 years
        climatology(:, :, doy) = reshape(daily_mean, n_lon, n_lat);
        daily_count(doy) = size(all_obs, 2);
        if mod(doy, 50) == 0
            fprintf('    DOY %3d: %d observations pooled\n', doy, daily_count(doy));
        end
    else
        climatology(:, :, doy) = NaN;
        daily_count(doy) = 0;
        fprintf('    DOY %3d: No data\n', doy);
    end
end

fprintf('✓ True daily climatology calculated for 365 days\n');
fprintf('  Daily climatology size: %s\n\n', mat2str(size(climatology)));

clear daily_data_pool;

% ========================================================================
% APPLY GAUSSIAN SMOOTHING TO REDUCE NOISE
% ========================================================================

fprintf('Applying Gaussian smoothing to reduce noise...\n');

climatology_smooth = zeros(n_lon, n_lat, 365, 'single');
window_length = 7;  % 7-day smoothing window

for i = 1:n_lon
    if mod(i, 50) == 0
        fprintf('  Processing longitude point %d/%d\n', i, n_lon);
    end
    
    for j = 1:n_lat
        ts = squeeze(climatology(i, j, :));
        if ~all(isnan(ts))
            % Apply Gaussian filter for smoothing
            climatology_smooth(i, j, :) = smoothdata(ts, 'gaussian', window_length);
        else
            climatology_smooth(i, j, :) = NaN;
        end
    end
end

climatology = climatology_smooth;

fprintf('✓ Climatology smoothing complete (window: %d days)\n\n', window_length);

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

all_clim_data = climatology(~isnan(climatology));
if ~isempty(all_clim_data)
    fprintf('  Overall climatology range:\n');
    fprintf('    Min: %.6f m²/s\n', min(all_clim_data));
    fprintf('    Max: %.6f m²/s\n', max(all_clim_data));
    fprintf('    Mean: %.6f m²/s\n', mean(all_clim_data));
    fprintf('    Std: %.6f m²/s\n\n', std(all_clim_data));
end

% ========================================================================
% CALCULATE SPATIAL AND TEMPORAL STATISTICS
% ========================================================================

fprintf('Calculating spatial and temporal statistics...\n');

daily_ui_mean = zeros(file_count, 1);
daily_clim_mean = zeros(file_count, 1);
daily_anom_mean = zeros(file_count, 1);
dates_datetime = zeros(file_count, 1);

for f = 1:file_count
    doy = day_of_year_all(f);
    
    ui_spatial = all_data{f};
    ui_spatial_mean = squeeze(mean(ui_spatial, 3, 'omitnan'));
    daily_ui_mean(f) = mean(ui_spatial_mean(~isnan(ui_spatial_mean)));
    
    clim_spatial = climatology(:, :, doy);
    daily_clim_mean(f) = mean(clim_spatial(~isnan(clim_spatial)));
    
    daily_anom_mean(f) = daily_ui_mean(f) - daily_clim_mean(f);
    
    date_str = dates_all{f};
    dates_datetime(f) = datenum(date_str, 'yyyymmdd');
    
    if mod(f, 200) == 0
        fprintf('  Processed file %d/%d\n', f, file_count);
    end
end

dates_datetime = datetime(dates_datetime, 'ConvertFrom', 'datenum');

fprintf('✓ Spatial statistics calculated\n\n');

% ========================================================================
% GET MONTH FOR EACH FILE
% ========================================================================

month_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    month_all(f) = month_val;
end

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
    month_mask = month_all == m;
    month_ui = daily_ui_mean(month_mask);
    month_anom = daily_anom_mean(month_mask);
    
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
% PLOTTING
% ========================================================================

fprintf('Creating plots...\n\n');

% Plot 1: Climatology time series at center point
i_plot = round(n_lon / 2);
j_plot = round(n_lat / 2);

figure('Position', [100, 100, 1400, 900]);

% Subplot 1: Daily climatology (full year)
subplot(2, 3, 1);
daily_clim_ts = squeeze(climatology(i_plot, j_plot, :));
plot(1:365, daily_clim_ts, '-', 'LineWidth', 2, 'Color', [0, 0.7, 0.3]);
xlabel('Day of Year', 'FontSize', 12);
ylabel('UI Ekman (m²/s)', 'FontSize', 12);
title('Daily Climatology (365 days, smoothed)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Subplot 2: Number of observations per day
subplot(2, 3, 2);
bar(1:365, daily_count, 'FaceColor', [0.7, 0.7, 1], 'EdgeColor', 'none');
xlabel('Day of Year', 'FontSize', 12);
ylabel('Number of Observations', 'FontSize', 12);
title('Data Coverage per Day', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Subplot 3: Monthly UI statistics
subplot(2, 3, 3);
errorbar(1:12, monthly_stats.ui_mean, monthly_stats.ui_std, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', [1, 0, 0]);
xlabel('Month', 'FontSize', 12);
ylabel('UI Ekman (m²/s)', 'FontSize', 12);
title('Monthly Mean ± Std Dev', 'FontSize', 12, 'FontWeight', 'bold');
grid on;
set(gca, 'XTick', 1:12);

% Subplot 4: Yearly time series with std bands (no Alpha property)
subplot(2, 3, 4);
upper_bound = yearly_stats.ui_mean + yearly_stats.ui_std;
lower_bound = yearly_stats.ui_mean - yearly_stats.ui_std;

% Plot std bands as thin dashed lines
plot(yearly_stats.years, upper_bound, '--', 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
hold on;
plot(yearly_stats.years, lower_bound, '--', 'LineWidth', 1, 'Color', [0.7, 0.7, 0.7]);
% Plot mean as main line
plot(yearly_stats.years, yearly_stats.ui_mean, 'o-', 'LineWidth', 2.5, 'MarkerSize', 6, 'Color', [0.5, 0, 0.5]);
hold off;

xlabel('Year', 'FontSize', 12);
ylabel('UI Ekman (m²/s)', 'FontSize', 12);
title('Yearly Mean UI ± Std Dev', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Subplot 5: Daily anomalies time series
subplot(2, 3, 5);
plot(dates_datetime, daily_anom_mean, '-', 'LineWidth', 0.5, 'Color', [0.7, 0.7, 0.7]);
hold on;
plot(dates_datetime, movmean(daily_anom_mean, 30), '-', 'LineWidth', 2, 'Color', [1, 0, 0]);
hold off;
xlabel('Date', 'FontSize', 12);
ylabel('Anomaly (m²/s)', 'FontSize', 12);
title('Daily Anomalies (red = 30-day moving mean)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

% Subplot 6: Statistics summary text
subplot(2, 3, 6);
axis off;
stats_text = sprintf(['Overall Statistics (2005-2024)\n\n' ...
    'UI Ekman:\n' ...
    '  Mean: %.4f m²/s\n' ...
    '  Std Dev: %.4f m²/s\n' ...
    '  Min: %.4f m²/s\n' ...
    '  Max: %.4f m²/s\n' ...
    '  Median: %.4f m²/s\n\n' ...
    'Anomalies:\n' ...
    '  Mean: %.4f m²/s\n' ...
    '  Std Dev: %.4f m²/s\n' ...
    '  Min: %.4f m²/s\n' ...
    '  Max: %.4f m²/s\n\n' ...
    'Total Observations: %d'], ...
    overall_stats.ui_mean, overall_stats.ui_std, ...
    overall_stats.ui_min, overall_stats.ui_max, overall_stats.ui_median, ...
    overall_stats.anom_mean, overall_stats.anom_std, ...
    overall_stats.anom_min, overall_stats.anom_max, ...
    overall_stats.n_observations);
text(0.1, 0.5, stats_text, 'FontSize', 11, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

sgtitle(sprintf('Daily Climatology Analysis - Grid Point (lon=%d, lat=%d)', i_plot, j_plot), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig(fullfile(plot_path, 'daily_climatology_analysis.fig'));
print(fullfile(plot_path, 'daily_climatology_analysis.png'), '-dpng', '-r150');
fprintf('✓ Main analysis plot saved\n');

% ========================================================================
% SAVE RESULTS
% ========================================================================

fprintf('\nSaving results...\n');

save(fullfile(output_path, 'daily_climatology_data.mat'), ...
    'climatology', 'daily_count', 'monthly_stats', 'yearly_stats', 'overall_stats', ...
    'daily_ui_mean', 'daily_clim_mean', 'daily_anom_mean', 'dates_datetime', ...
    'lon', 'lat', 'years', '-v7.3');

fprintf('✓ Results saved to daily_climatology_data.mat\n\n');

fprintf('========================================================================\n');
fprintf('✓ DAILY CLIMATOLOGY ANALYSIS COMPLETE\n');
fprintf('========================================================================\n');