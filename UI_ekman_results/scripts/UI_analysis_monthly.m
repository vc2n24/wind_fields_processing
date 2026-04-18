clear; clc; close all;

% ========================================================================
% CALCULATE UI_EKMAN MONTHLY CLIMATOLOGY, ANOMALIES & STATISTICS (2005-2024)
% ========================================================================

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman_results';
output_path = fullfile(ui_results_path, 'analysis');
plot_path = fullfile(output_path, 'plots', 'monthly_climatologies');

% Create output directories
if ~isfolder(output_path)
    mkdir(output_path);
end
if ~isfolder(plot_path)
    mkdir(plot_path);
end

% ========================================================================
% PARAMETERS
% ========================================================================

mask_name = '50km';
years = 2005:2024;  % All years
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
    
    % Determine if leap year
    if mod(current_year, 4) == 0 && (mod(current_year, 100) ~= 0 || mod(current_year, 400) == 0)
        days_in_year = 366;
    else
        days_in_year = 365;
    end
    
    for m = 1:12
        current_month = m;
        Mstr = sprintf('M%02d', current_month);
        
        % Get number of days in this month
        if current_month == 2
            if days_in_year == 366
                n_days = 29;
            else
                n_days = 28;
            end
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

[n_time, n_lat, n_lon] = size(UI_data);
fprintf('  Time steps per day: %d\n\n', n_time);

% ========================================================================
% REMOVE LEAP DAYS (FEB 29)
% ========================================================================

fprintf('Removing leap days (Feb 29)...\n');

leap_day_indices = [];
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    day_val = str2double(date_str(7:8));
    
    if month_val == 2 && day_val == 29
        leap_day_indices = [leap_day_indices, f];
    end
end

% Remove leap day files
all_data(leap_day_indices) = [];
dates_all(leap_day_indices) = [];
file_count = length(all_data);

fprintf('  Removed %d leap days (Feb 29)\n', length(leap_day_indices));
fprintf('  Files remaining after removal: %d\n\n', file_count);

% ========================================================================
% CALCULATE MONTHLY CLIMATOLOGY
% ========================================================================

fprintf('Calculating monthly climatology (month average)...\n');

% Month and year for each date
month_year_all = zeros(file_count, 2);  % [month, year]
for f = 1:file_count
    date_str = dates_all{f};
    year_val = str2double(date_str(1:4));
    month_val = str2double(date_str(5:6));
    
    month_year_all(f, 1) = month_val;
    month_year_all(f, 2) = year_val;
end

% Initialize climatology: [month, time, lat, lon]
climatology = zeros(12, n_time, n_lat, n_lon, 'single');
climatology_count = zeros(12, 1);

% Calculate climatology (mean for each month)
for f = 1:file_count
    month = month_year_all(f, 1);
    % Reshape to properly handle dimensions
    climatology(month, :, :, :) = climatology(month, :, :, :) + reshape(single(all_data{f}), 1, n_time, n_lat, n_lon);
    climatology_count(month) = climatology_count(month) + 1;
end

% Average
for month = 1:12
    if climatology_count(month) > 0
        climatology(month, :, :, :) = climatology(month, :, :, :) / climatology_count(month);
    else
        climatology(month, :, :, :) = NaN;
    end
end

fprintf('✓ Monthly climatology calculated for 12 months\n');
fprintf('  Months with data: %d\n\n', sum(climatology_count > 0));

% ========================================================================
% CALCULATE SPATIAL STATISTICS
% ========================================================================

fprintf('Calculating spatial statistics...\n');

% Time series of spatially-averaged UI, climatology, and anomalies
UI_timeseries = zeros(file_count, n_time);
climatology_timeseries = zeros(file_count, n_time);
anomaly_timeseries = zeros(file_count, n_time);
dates_datetime = zeros(file_count, 1);

for f = 1:file_count
    % Spatial mean
    for t = 1:n_time
        % Valid pixels (not NaN)
        UI_valid = squeeze(all_data{f}(t, :, :));
        UI_valid = UI_valid(~isnan(UI_valid));
        
        if ~isempty(UI_valid)
            UI_timeseries(f, t) = mean(UI_valid);
        else
            UI_timeseries(f, t) = NaN;
        end
        
        % Climatology
        month = month_year_all(f, 1);
        clim_valid = squeeze(climatology(month, t, :, :));
        clim_valid = clim_valid(~isnan(clim_valid));
        if ~isempty(clim_valid)
            climatology_timeseries(f, t) = mean(clim_valid);
        else
            climatology_timeseries(f, t) = NaN;
        end
        
        % Anomaly
        anomaly_timeseries(f, t) = UI_timeseries(f, t) - climatology_timeseries(f, t);
    end
    
    % Convert date string to datetime
    date_str = dates_all{f};
    dates_datetime(f) = datenum(date_str, 'yyyymmdd');
end

% Convert to datetime
dates_datetime = datetime(dates_datetime, 'ConvertFrom', 'datenum');

fprintf('✓ Spatial statistics calculated\n\n');

% ========================================================================
% ADDITIONAL STATISTICS
% ========================================================================

fprintf('Computing additional statistics...\n');

% Extract year and month as vectors (not using datetime indexing directly)
years_vec = year(dates_datetime);
months_vec = month(dates_datetime);

% Statistics by year
yearly_mean = zeros(n_years, 1);
yearly_std = zeros(n_years, 1);
yearly_max = zeros(n_years, 1);
yearly_min = zeros(n_years, 1);

for y = 1:n_years
    current_year = years(y);
    % Find indices for this year (use logical indexing, not datetime)
    year_mask = (years_vec == current_year);
    year_data = UI_timeseries(year_mask, :);
    year_data_valid = year_data(~isnan(year_data));
    
    if ~isempty(year_data_valid)
        yearly_mean(y) = mean(year_data_valid);
        yearly_std(y) = std(year_data_valid);
        yearly_max(y) = max(year_data_valid);
        yearly_min(y) = min(year_data_valid);
    else
        yearly_mean(y) = NaN;
        yearly_std(y) = NaN;
        yearly_max(y) = NaN;
        yearly_min(y) = NaN;
    end
end

% Overall statistics
overall_mean = mean(UI_timeseries(~isnan(UI_timeseries)));
overall_std = std(UI_timeseries(~isnan(UI_timeseries)));
overall_max = max(UI_timeseries(~isnan(UI_timeseries)));
overall_min = min(UI_timeseries(~isnan(UI_timeseries)));

fprintf('✓ Additional statistics calculated\n');
fprintf('  Overall mean UI: %.4f m²/s\n', overall_mean);
fprintf('  Overall std UI: %.4f m²/s\n', overall_std);
fprintf('  Overall max UI: %.4f m²/s\n', overall_max);
fprintf('  Overall min UI: %.4f m²/s\n\n', overall_min);

% ========================================================================
% CREATE COMPREHENSIVE PLOTS
% ========================================================================

fprintf('Creating plots...\n');

% Set figure visibility to off for all plots
set(0, 'DefaultFigureVisible', 'off');

% Plot 1: Time series of monthly-mean UI (all years)
fig1 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
plot(dates_datetime, mean(UI_timeseries, 2), 'b-', 'LineWidth', 1);
hold on
plot(dates_datetime, mean(climatology_timeseries, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'Climatology');
xlabel('Time (2005-2024)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ekman Upwelling Index Time Series (50km mask) - 2005 to 2024', 'FontSize', 14, 'FontWeight', 'bold');
legend('Daily mean UI', 'Monthly climatology', 'Location', 'best', 'FontSize', 11);
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XAxis.TickLabelFormat = 'yyyy';
xtickangle(45);
% Save as PNG and FIG
saveas(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.png'), 'png');
saveas(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.fig'), 'fig');
% Save as JPEG
print(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.jpg'), '-djpeg', '-r150');
close(fig1);
fprintf('  ✓ Saved: UI_timeseries_2005_2024 (png, fig, jpg)\n');

% Plot 2: Monthly climatology (mean for each month)
fig2 = figure('Position', [100 100 1200 600], 'Color', 'w', 'Visible', 'off');
climatology_mean_month = squeeze(mean(climatology, [3, 4], 'omitnan'));  % Average over lat/lon
climatology_mean_month_daily = mean(climatology_mean_month, 2);  % Average over time steps
bar(1:12, climatology_mean_month_daily, 'FaceColor', [0.8, 0.2, 0.2], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
plot(1:12, climatology_mean_month_daily, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Month', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Monthly Climatology of UI_{Ekman} (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
% Set month labels
month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
set(gca, 'XTick', 1:12, 'XTickLabel', month_labels);
xlim([0.5, 12.5]);
% Save as PNG and FIG
saveas(fig2, fullfile(plot_path, 'UI_climatology_monthly.png'), 'png');
saveas(fig2, fullfile(plot_path, 'UI_climatology_monthly.fig'), 'fig');
% Save as JPEG
print(fig2, fullfile(plot_path, 'UI_climatology_monthly.jpg'), '-djpeg', '-r150');
close(fig2);
fprintf('  ✓ Saved: UI_climatology_monthly (png, fig, jpg)\n');

% Plot 3: Anomalies time series
fig3 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
anomalies_mean = mean(anomaly_timeseries, 2);
plot(dates_datetime, anomalies_mean, 'Color', [0.2, 0.6, 1], 'LineWidth', 1);
hold on
% Add zero line
plot(dates_datetime([1, end]), [0, 0], 'k--', 'LineWidth', 1.5);
% Positive anomalies in blue, negative in red
idx_pos = anomalies_mean >= 0;
idx_neg = anomalies_mean < 0;
scatter(dates_datetime(idx_pos), anomalies_mean(idx_pos), 20, 'b', 'filled', 'DisplayName', 'Positive anomalies');
scatter(dates_datetime(idx_neg), anomalies_mean(idx_neg), 20, 'r', 'filled', 'DisplayName', 'Negative anomalies');
xlabel('Time (2005-2024)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} Anomaly (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ekman Upwelling Index Anomalies (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
ax.XAxis.TickLabelFormat = 'yyyy';
xtickangle(45);
% Save as PNG and FIG
saveas(fig3, fullfile(plot_path, 'UI_anomalies_timeseries.png'), 'png');
saveas(fig3, fullfile(plot_path, 'UI_anomalies_timeseries.fig'), 'fig');
% Save as JPEG
print(fig3, fullfile(plot_path, 'UI_anomalies_timeseries.jpg'), '-djpeg', '-r150');
close(fig3);
fprintf('  ✓ Saved: UI_anomalies_timeseries (png, fig, jpg)\n');

% Plot 4: Yearly statistics box plot
fig4 = figure('Position', [100 100 1200 600], 'Color', 'w', 'Visible', 'off');
% Prepare data for boxplot - flatten all data into single array
all_year_data = [];
group_indices = [];

for y = 1:n_years
    current_year = years(y);
    year_mask = (years_vec == current_year);
    year_data_flat = UI_timeseries(year_mask, :);
    year_data_flat = year_data_flat(~isnan(year_data_flat));
    year_data_flat = year_data_flat(:);  % Ensure column vector
    
    all_year_data = [all_year_data; year_data_flat];
    group_indices = [group_indices; repmat(y, length(year_data_flat), 1)];
end

boxplot(all_year_data, group_indices, 'Labels', string(years));
xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Annual Distribution of UI_{Ekman} (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
xtickangle(45);
% Save as PNG and FIG
saveas(fig4, fullfile(plot_path, 'UI_yearly_boxplot.png'), 'png');
saveas(fig4, fullfile(plot_path, 'UI_yearly_boxplot.fig'), 'fig');
% Save as JPEG
print(fig4, fullfile(plot_path, 'UI_yearly_boxplot.jpg'), '-djpeg', '-r150');
close(fig4);
fprintf('  ✓ Saved: UI_yearly_boxplot (png, fig, jpg)\n');

% Plot 5: Yearly mean with error bars
fig5 = figure('Position', [100 100 1200 600], 'Color', 'w', 'Visible', 'off');
errorbar(years, yearly_mean, yearly_std, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'CapSize', 5);
xlabel('Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Annual Mean UI_{Ekman} with Standard Deviation (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
% Save as PNG and FIG
saveas(fig5, fullfile(plot_path, 'UI_yearly_mean_std.png'), 'png');
saveas(fig5, fullfile(plot_path, 'UI_yearly_mean_std.fig'), 'fig');
% Save as JPEG
print(fig5, fullfile(plot_path, 'UI_yearly_mean_std.jpg'), '-djpeg', '-r150');
close(fig5);
fprintf('  ✓ Saved: UI_yearly_mean_std (png, fig, jpg)\n');

% Plot 6: Monthly box plot for each month across all years
fig6 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
% Prepare data for monthly boxplot
all_month_data = [];
month_group_indices = [];

month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

for month_idx = 1:12
    month_mask = (months_vec == month_idx);
    month_data_flat = UI_timeseries(month_mask, :);
    month_data_flat = month_data_flat(~isnan(month_data_flat));
    month_data_flat = month_data_flat(:);  % Ensure column vector
    
    all_month_data = [all_month_data; month_data_flat];
    month_group_indices = [month_group_indices; repmat(month_idx, length(month_data_flat), 1)];
end

boxplot(all_month_data, month_group_indices, 'Labels', month_labels);
xlabel('Month', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Monthly Distribution of UI_{Ekman} (50km mask, 2005-2024)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
xtickangle(45);
% Save as PNG and FIG
saveas(fig6, fullfile(plot_path, 'UI_monthly_boxplot.png'), 'png');
saveas(fig6, fullfile(plot_path, 'UI_monthly_boxplot.fig'), 'fig');
% Save as JPEG
print(fig6, fullfile(plot_path, 'UI_monthly_boxplot.jpg'), '-djpeg', '-r150');
close(fig6);
fprintf('  ✓ Saved: UI_monthly_boxplot (png, fig, jpg)\n');

% ========================================================================
% SAVE CLIMATOLOGY AND ANOMALIES TO NETCDF
% ========================================================================

fprintf('Saving climatology and anomalies to NetCDF files...\n');

% Save climatology
climatology_filepath = fullfile(output_path, sprintf('UI_ekman_monthly_climatology_%s.nc', mask_name));
if isfile(climatology_filepath)
    delete(climatology_filepath);
end

mode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(climatology_filepath, mode);

% Define dimensions
month_dimid = netcdf.defDim(ncid, 'month', 12);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);

% Define variables
month_varid = netcdf.defVar(ncid, 'month', 'NC_INT', month_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
count_varid = netcdf.defVar(ncid, 'count', 'NC_INT', month_dimid);
climatology_varid = netcdf.defVar(ncid, 'UI_ekman_climatology', 'NC_FLOAT', ...
    [month_dimid, time_dimid, lat_dimid, lon_dimid]);

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid, month_varid, int32(1:12));
netcdf.putVar(ncid, time_varid, int32(1:n_time));
netcdf.putVar(ncid, lat_varid, double(lat));
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, count_varid, int32(climatology_count));
netcdf.putVar(ncid, climatology_varid, climatology);

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, climatology_varid, 'long_name', 'Monthly climatology of UI_Ekman');
netcdf.putAtt(ncid, climatology_varid, 'units', 'm^2 s^-1');
netcdf.putAtt(ncid, climatology_varid, 'description', 'Mean UI for each month (2005-2024, leap days removed)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Monthly Climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');

netcdf.close(ncid);

fprintf('✓ Saved climatology: %s\n', climatology_filepath);

% Save time series and statistics
timeseries_filepath = fullfile(output_path, sprintf('UI_ekman_timeseries_statistics_monthly_%s.nc', mask_name));
if isfile(timeseries_filepath)
    delete(timeseries_filepath);
end

ncid = netcdf.create(timeseries_filepath, mode);

% Define dimensions
file_dimid = netcdf.defDim(ncid, 'file', file_count);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
year_dimid = netcdf.defDim(ncid, 'year', n_years);
month_dimid = netcdf.defDim(ncid, 'month', 12);

% Define variables
month_year_varid = netcdf.defVar(ncid, 'month_year', 'NC_INT', [file_dimid, 2]);
ui_ts_varid = netcdf.defVar(ncid, 'UI_timeseries', 'NC_FLOAT', ...
    [file_dimid, time_dimid]);
anom_ts_varid = netcdf.defVar(ncid, 'anomaly_timeseries', 'NC_FLOAT', ...
    [file_dimid, time_dimid]);
year_varid = netcdf.defVar(ncid, 'years', 'NC_INT', year_dimid);
yearly_mean_varid = netcdf.defVar(ncid, 'yearly_mean', 'NC_FLOAT', year_dimid);
yearly_std_varid = netcdf.defVar(ncid, 'yearly_std', 'NC_FLOAT', year_dimid);

netcdf.endDef(ncid);

% Write data
netcdf.putVar(ncid, month_year_varid, int32(month_year_all));
netcdf.putVar(ncid, ui_ts_varid, single(UI_timeseries));
netcdf.putVar(ncid, anom_ts_varid, single(anomaly_timeseries));
netcdf.putVar(ncid, year_varid, int32(years));
netcdf.putVar(ncid, yearly_mean_varid, single(yearly_mean));
netcdf.putVar(ncid, yearly_std_varid, single(yearly_std));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, ui_ts_varid, 'long_name', 'Spatially-averaged UI_Ekman time series');
netcdf.putAtt(ncid, anom_ts_varid, 'long_name', 'UI_Ekman anomalies (daily - monthly climatology)');
netcdf.putAtt(ncid, yearly_mean_varid, 'long_name', 'Annual mean UI_Ekman');
netcdf.putAtt(ncid, yearly_std_varid, 'long_name', 'Annual standard deviation of UI_Ekman');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Time Series and Statistics (Monthly Climatology)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_mean', overall_mean);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_std', overall_std);

netcdf.close(ncid);

fprintf('✓ Saved time series: %s\n\n', timeseries_filepath);

% ========================================================================
% SUMMARY
% ========================================================================

fprintf('\n========================================\n');
fprintf('✓ ANALYSIS COMPLETE!\n');
fprintf('========================================\n');
fprintf('Results saved to: %s\n', output_path);
fprintf('Plots saved to: %s\n', plot_path);
fprintf('\nStatistics Summary (2005-2024, leap days removed):\n');
fprintf('  Overall mean UI: %.4f m²/s\n', overall_mean);
fprintf('  Overall std: %.4f m²/s\n', overall_std);
fprintf('  Overall max: %.4f m²/s (upwelling favorable)\n', overall_max);
fprintf('  Overall min: %.4f m²/s (downwelling favorable)\n', overall_min);
fprintf('  Files processed: %d\n', file_count);
fprintf('========================================\n');