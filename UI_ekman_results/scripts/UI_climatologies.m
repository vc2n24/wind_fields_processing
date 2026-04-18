% ========================================================================
% CALCULATE UI_EKMAN CLIMATOLOGY, ANOMALIES & STATISTICS (2005-2024)
% ========================================================================
clear all; close all; clc;

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman_results';
output_path = fullfile(ui_results_path, 'analysis');
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
years = 2005:2006;  % All years
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
% CALCULATE DAILY CLIMATOLOGY
% ========================================================================

fprintf('Calculating daily climatology (day of year average)...\n');

% Day of year for each date (adjust for removed leap days)
day_of_year_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    year_val = str2double(date_str(1:4));
    month_val = str2double(date_str(5:6));
    day_val = str2double(date_str(7:8));
    
    % Calculate day of year (without Feb 29)
    days_in_months = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
    day_of_year = sum(days_in_months(1:month_val-1)) + day_val;
    day_of_year_all(f) = day_of_year;
end

% Initialize climatology: [day_of_year, time, lat, lon] - 365 days only
climatology = zeros(365, n_time, n_lat, n_lon, 'single');
climatology_count = zeros(365, 1);

% Calculate climatology (mean for each day of year)
for f = 1:file_count
    doy = day_of_year_all(f);
    % Reshape to properly handle dimensions
    climatology(doy, :, :, :) = climatology(doy, :, :, :) + reshape(single(all_data{f}), 1, n_time, n_lat, n_lon);
    climatology_count(doy) = climatology_count(doy) + 1;
end

% Average
for doy = 1:365
    if climatology_count(doy) > 0
        climatology(doy, :, :, :) = climatology(doy, :, :, :) / climatology_count(doy);
    else
        climatology(doy, :, :, :) = NaN;
    end
end

fprintf('✓ Climatology calculated for 365 days\n');
fprintf('  Days with data: %d\n\n', sum(climatology_count > 0));

% ========================================================================
% CALCULATE ANOMALIES (FROM DAILY CLIMATOLOGY)
% ========================================================================

fprintf('Calculating anomalies from daily climatology...\n');

anomalies = cell(file_count, 1);
for f = 1:file_count
    doy = day_of_year_all(f);
    % Proper dimension handling
    climatology_field = squeeze(climatology(doy, :, :, :));
    anomalies{f} = single(all_data{f}) - climatology_field;
end

fprintf('✓ Anomalies from daily climatology calculated\n\n');

% ========================================================================
% CALCULATE SPATIAL STATISTICS (DAILY)
% ========================================================================

fprintf('Calculating spatial statistics from daily climatology...\n');

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
        doy = day_of_year_all(f);
        clim_valid = squeeze(climatology(doy, t, :, :));
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

fprintf('✓ Spatial statistics from daily climatology calculated\n\n');

% ========================================================================
% CALCULATE ADDITIONAL STATISTICS
% ========================================================================

fprintf('Computing additional statistics...\n');

% Statistics by year
yearly_mean = zeros(n_years, 1);
yearly_std = zeros(n_years, 1);
yearly_max = zeros(n_years, 1);
yearly_min = zeros(n_years, 1);

for y = 1:n_years
    current_year = years(y);
    % Find indices for this year
    year_mask = year(dates_datetime) == current_year;
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
% CREATE COMPREHENSIVE PLOTS (DAILY CLIMATOLOGY)
% ========================================================================

fprintf('Creating plots from daily climatology...\n');

% Set figure visibility to off for all plots
set(0, 'DefaultFigureVisible', 'off');

% Plot 1: Time series of daily-mean UI (all years)
fig1 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
plot(dates_datetime, mean(UI_timeseries, 2), 'b-', 'LineWidth', 1);
hold on
plot(dates_datetime, mean(climatology_timeseries, 2), 'r-', 'LineWidth', 2, 'DisplayName', 'Climatology');
xlabel('Time (2005-2024)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Ekman Upwelling Index Time Series (50km mask) - 2005 to 2024', 'FontSize', 14, 'FontWeight', 'bold');
legend('Daily mean UI', 'Daily climatology', 'Location', 'best', 'FontSize', 11);
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
% Set x-axis to show year labels at regular intervals
ax.XAxis.TickLabelFormat = 'yyyy';
xtickangle(45);
% Save as PNG and FIG
saveas(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.png'), 'png');
saveas(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.fig'), 'fig');
% Save as JPEG
print(fig1, fullfile(plot_path, 'UI_timeseries_2005_2024.jpg'), '-djpeg', '-r150');
close(fig1);
fprintf('  ✓ Saved: UI_timeseries_2005_2024 (png, fig, jpg)\n');

% Plot 2: Climatology (mean day of year)
fig2 = figure('Position', [100 100 1200 600], 'Color', 'w', 'Visible', 'off');
climatology_mean_doy = squeeze(mean(climatology, [3, 4], 'omitnan'));  % Average over lat/lon
climatology_mean_doy_daily = mean(climatology_mean_doy, 2);  % Average over time steps
plot(1:365, climatology_mean_doy_daily, 'r-', 'LineWidth', 2);
xlabel('Day of Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Daily Climatology of UI_{Ekman} (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
% Add month labels (365 days, no leap day)
month_days = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
set(gca, 'XTick', month_days, 'XTickLabel', month_labels);
xlim([1, 365]);
% Save as PNG and FIG
saveas(fig2, fullfile(plot_path, 'UI_climatology_daily.png'), 'png');
saveas(fig2, fullfile(plot_path, 'UI_climatology_daily.fig'), 'fig');
% Save as JPEG
print(fig2, fullfile(plot_path, 'UI_climatology_daily.jpg'), '-djpeg', '-r150');
close(fig2);
fprintf('  ✓ Saved: UI_climatology_daily (png, fig, jpg)\n');

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
    year_mask = year(dates_datetime) == current_year;
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

fprintf('✓ All daily climatology plots created\n\n');

% ========================================================================
% ========================================================================
% SAVE DAILY CLIMATOLOGY AND ANOMALIES TO NETCDF
% ========================================================================

fprintf('Saving daily climatology and anomalies to NetCDF files...\n');

% Save daily climatology
climatology_filepath = fullfile(output_path, sprintf('UI_ekman_daily_climatology_%s.nc', mask_name));
if isfile(climatology_filepath)
    delete(climatology_filepath);
end

mode = bitor(netcdf.getConstant('CLOBBER'), netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(climatology_filepath, mode);

% Define dimensions
doy_dimid = netcdf.defDim(ncid, 'day_of_year', 365);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);

% Define variables
doy_varid = netcdf.defVar(ncid, 'day_of_year', 'NC_INT', doy_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
count_varid = netcdf.defVar(ncid, 'count', 'NC_INT', doy_dimid);
% Define with dimensions in NetCDF order: [doy, time, lat, lon]
climatology_varid = netcdf.defVar(ncid, 'UI_ekman_climatology', 'NC_FLOAT', ...
    [doy_dimid, time_dimid, lat_dimid, lon_dimid]);

netcdf.endDef(ncid);

% Write data - NO permute, write as MATLAB native order [doy, time, lat, lon]
netcdf.putVar(ncid, doy_varid, int32(1:365));
netcdf.putVar(ncid, time_varid, int32(1:n_time));
netcdf.putVar(ncid, lat_varid, double(lat));
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, count_varid, int32(climatology_count));
netcdf.putVar(ncid, climatology_varid, single(climatology));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, climatology_varid, 'long_name', 'Daily climatology of UI_Ekman');
netcdf.putAtt(ncid, climatology_varid, 'units', 'm^2 s^-1');
netcdf.putAtt(ncid, climatology_varid, 'description', 'Mean UI for each day of year (2005-2024, leap days removed). Note: MATLAB stores in column-major order, so when read externally dimensions are [time, lat, lon, doy]');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Daily Climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');

netcdf.close(ncid);

fprintf('✓ Saved daily climatology: %s\n', climatology_filepath);

% Save time series and statistics
timeseries_filepath = fullfile(output_path, sprintf('UI_ekman_timeseries_statistics_%s.nc', mask_name));
if isfile(timeseries_filepath)
    delete(timeseries_filepath);
end

ncid = netcdf.create(timeseries_filepath, mode);

% Define dimensions
file_dimid = netcdf.defDim(ncid, 'file', file_count);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
year_dimid = netcdf.defDim(ncid, 'year', n_years);

% Define variables
doy_varid = netcdf.defVar(ncid, 'day_of_year', 'NC_INT', file_dimid);
ui_ts_varid = netcdf.defVar(ncid, 'UI_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
anom_ts_varid = netcdf.defVar(ncid, 'anomaly_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
year_varid = netcdf.defVar(ncid, 'years', 'NC_INT', year_dimid);
yearly_mean_varid = netcdf.defVar(ncid, 'yearly_mean', 'NC_FLOAT', year_dimid);
yearly_std_varid = netcdf.defVar(ncid, 'yearly_std', 'NC_FLOAT', year_dimid);

netcdf.endDef(ncid);

% Write data - NO permute, write as MATLAB native order [file, time]
netcdf.putVar(ncid, doy_varid, int32(day_of_year_all));
netcdf.putVar(ncid, ui_ts_varid, single(UI_timeseries));
netcdf.putVar(ncid, anom_ts_varid, single(anomaly_timeseries));
netcdf.putVar(ncid, year_varid, int32(years));
netcdf.putVar(ncid, yearly_mean_varid, single(yearly_mean));
netcdf.putVar(ncid, yearly_std_varid, single(yearly_std));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, ui_ts_varid, 'long_name', 'Spatially-averaged UI_Ekman time series');
netcdf.putAtt(ncid, anom_ts_varid, 'long_name', 'UI_Ekman anomalies (daily - daily climatology)');
netcdf.putAtt(ncid, yearly_mean_varid, 'long_name', 'Annual mean UI_Ekman');
netcdf.putAtt(ncid, yearly_std_varid, 'long_name', 'Annual standard deviation of UI_Ekman');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Time Series and Statistics');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_mean', overall_mean);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_std', overall_std);

netcdf.close(ncid);

fprintf('✓ Saved daily time series: %s\n\n', timeseries_filepath);

% Save time series and statistics
timeseries_filepath = fullfile(output_path, sprintf('UI_ekman_timeseries_statistics_%s.nc', mask_name));
if isfile(timeseries_filepath)
    delete(timeseries_filepath);
end

ncid = netcdf.create(timeseries_filepath, mode);

% Define dimensions
file_dimid = netcdf.defDim(ncid, 'file', file_count);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
year_dimid = netcdf.defDim(ncid, 'year', n_years);

% Define variables - permute to [time, file] for MATLAB
doy_varid = netcdf.defVar(ncid, 'day_of_year', 'NC_INT', file_dimid);
ui_ts_varid = netcdf.defVar(ncid, 'UI_timeseries', 'NC_FLOAT', [time_dimid, file_dimid]);
anom_ts_varid = netcdf.defVar(ncid, 'anomaly_timeseries', 'NC_FLOAT', [time_dimid, file_dimid]);
year_varid = netcdf.defVar(ncid, 'years', 'NC_INT', year_dimid);
yearly_mean_varid = netcdf.defVar(ncid, 'yearly_mean', 'NC_FLOAT', year_dimid);
yearly_std_varid = netcdf.defVar(ncid, 'yearly_std', 'NC_FLOAT', year_dimid);

netcdf.endDef(ncid);

% Write data - permute from [file, time] to [time, file]
netcdf.putVar(ncid, doy_varid, int32(day_of_year_all));
netcdf.putVar(ncid, ui_ts_varid, single(permute(UI_timeseries, [2, 1])));
netcdf.putVar(ncid, anom_ts_varid, single(permute(anomaly_timeseries, [2, 1])));
netcdf.putVar(ncid, year_varid, int32(years));
netcdf.putVar(ncid, yearly_mean_varid, single(yearly_mean));
netcdf.putVar(ncid, yearly_std_varid, single(yearly_std));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, ui_ts_varid, 'long_name', 'Spatially-averaged UI_Ekman time series');
netcdf.putAtt(ncid, anom_ts_varid, 'long_name', 'UI_Ekman anomalies (daily - daily climatology)');
netcdf.putAtt(ncid, yearly_mean_varid, 'long_name', 'Annual mean UI_Ekman');
netcdf.putAtt(ncid, yearly_std_varid, 'long_name', 'Annual standard deviation of UI_Ekman');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Time Series and Statistics');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_mean', overall_mean);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'overall_std', overall_std);

netcdf.close(ncid);

fprintf('✓ Saved daily time series: %s\n\n', timeseries_filepath);

% ========================================================================
% CALCULATE MONTHLY CLIMATOLOGY
% ========================================================================

fprintf('Calculating monthly climatology (calendar month average)...\n');

% Month values for each file
month_all = zeros(file_count, 1);
for f = 1:file_count
    date_str = dates_all{f};
    month_val = str2double(date_str(5:6));
    month_all(f) = month_val;
end

% Initialize monthly climatology: [month, time, lat, lon]
monthly_climatology_raw = zeros(12, n_time, n_lat, n_lon, 'single');
monthly_count = zeros(12, 1);

% Calculate monthly climatology (mean for each calendar month)
for f = 1:file_count
    m = month_all(f);
    % Accumulate
    monthly_climatology_raw(m, :, :, :) = monthly_climatology_raw(m, :, :, :) + reshape(single(all_data{f}), 1, n_time, n_lat, n_lon);
    monthly_count(m) = monthly_count(m) + 1;
end

% Average
for m = 1:12
    if monthly_count(m) > 0
        monthly_climatology_raw(m, :, :, :) = monthly_climatology_raw(m, :, :, :) / monthly_count(m);
    else
        monthly_climatology_raw(m, :, :, :) = NaN;
    end
end

fprintf('✓ Monthly climatology calculated for 12 months\n');
fprintf('  Months with data: %d\n\n', sum(monthly_count > 0));

% ========================================================================
% CALCULATE MONTHLY ANOMALIES
% ========================================================================

fprintf('Calculating anomalies from monthly climatology...\n');

anomalies_monthly = cell(file_count, 1);
for f = 1:file_count
    m = month_all(f);
    % Proper dimension handling
    climatology_field = squeeze(monthly_climatology_raw(m, :, :, :));
    anomalies_monthly{f} = single(all_data{f}) - climatology_field;
end

fprintf('✓ Anomalies from monthly climatology calculated\n\n');

% ========================================================================
% CALCULATE SPATIAL STATISTICS (MONTHLY)
% ========================================================================

fprintf('Calculating spatial statistics from monthly climatology...\n');

% Time series of spatially-averaged UI and anomalies (monthly)
monthly_climatology_timeseries = zeros(file_count, n_time);
anomaly_timeseries_monthly = zeros(file_count, n_time);

for f = 1:file_count
    % Spatial mean
    for t = 1:n_time
        % Monthly climatology
        m = month_all(f);
        clim_valid = squeeze(monthly_climatology_raw(m, t, :, :));
        clim_valid = clim_valid(~isnan(clim_valid));
        if ~isempty(clim_valid)
            monthly_climatology_timeseries(f, t) = mean(clim_valid);
        else
            monthly_climatology_timeseries(f, t) = NaN;
        end
        
        % Anomaly
        anomaly_timeseries_monthly(f, t) = UI_timeseries(f, t) - monthly_climatology_timeseries(f, t);
    end
end

fprintf('✓ Spatial statistics from monthly climatology calculated\n\n');

% ========================================================================
% ========================================================================
% SAVE MONTHLY CLIMATOLOGY TO NETCDF
% ========================================================================

fprintf('Saving monthly climatology to NetCDF file...\n');

monthly_climatology_filepath = fullfile(output_path, sprintf('UI_ekman_monthly_climatology_%s.nc', mask_name));
if isfile(monthly_climatology_filepath)
    delete(monthly_climatology_filepath);
end

ncid = netcdf.create(monthly_climatology_filepath, mode);

% Define dimensions
month_dimid = netcdf.defDim(ncid, 'month', 12);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);

% Define variables - NO permute, write as MATLAB native order [month, time, lat, lon]
month_varid = netcdf.defVar(ncid, 'month', 'NC_INT', month_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
count_varid = netcdf.defVar(ncid, 'count', 'NC_INT', month_dimid);
monthly_climatology_varid = netcdf.defVar(ncid, 'UI_ekman_monthly_climatology', 'NC_FLOAT', ...
    [month_dimid, time_dimid, lat_dimid, lon_dimid]);

netcdf.endDef(ncid);

% Write data - NO permute
netcdf.putVar(ncid, month_varid, int32(1:12));
netcdf.putVar(ncid, time_varid, int32(1:n_time));
netcdf.putVar(ncid, lat_varid, double(lat));
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, count_varid, int32(monthly_count));
netcdf.putVar(ncid, monthly_climatology_varid, single(monthly_climatology_raw));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, monthly_climatology_varid, 'long_name', 'Monthly climatology of UI_Ekman');
netcdf.putAtt(ncid, monthly_climatology_varid, 'units', 'm^2 s^-1');
netcdf.putAtt(ncid, monthly_climatology_varid, 'description', 'Mean UI for each calendar month (2005-2024). Note: MATLAB stores in column-major order, so when read externally dimensions are [time, lat, lon, month]');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Monthly Climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');

netcdf.close(ncid);

fprintf('✓ Saved monthly climatology: %s\n\n', monthly_climatology_filepath);

% ========================================================================
% SAVE MONTHLY ANOMALIES TO NETCDF
% ========================================================================

fprintf('Saving monthly anomalies to NetCDF file...\n');

monthly_timeseries_filepath = fullfile(output_path, sprintf('UI_ekman_monthly_timeseries_statistics_%s.nc', mask_name));
if isfile(monthly_timeseries_filepath)
    delete(monthly_timeseries_filepath);
end

ncid = netcdf.create(monthly_timeseries_filepath, mode);

% Define dimensions
file_dimid = netcdf.defDim(ncid, 'file', file_count);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
year_dimid = netcdf.defDim(ncid, 'year', n_years);

% Define variables - NO permute
ui_ts_varid = netcdf.defVar(ncid, 'UI_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
clim_ts_varid = netcdf.defVar(ncid, 'monthly_climatology_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
anom_ts_varid = netcdf.defVar(ncid, 'anomaly_timeseries_monthly', 'NC_FLOAT', [file_dimid, time_dimid]);
year_varid = netcdf.defVar(ncid, 'years', 'NC_INT', year_dimid);

netcdf.endDef(ncid);

% Write data - NO permute
netcdf.putVar(ncid, ui_ts_varid, single(UI_timeseries));
netcdf.putVar(ncid, clim_ts_varid, single(monthly_climatology_timeseries));
netcdf.putVar(ncid, anom_ts_varid, single(anomaly_timeseries_monthly));
netcdf.putVar(ncid, year_varid, int32(years));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, anom_ts_varid, 'long_name', 'UI_Ekman anomalies from monthly climatology');
netcdf.putAtt(ncid, anom_ts_varid, 'description', 'Daily UI - calendar month climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Monthly Climatology Time Series and Anomalies');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');

netcdf.close(ncid);

fprintf('✓ Saved monthly time series: %s\n\n', monthly_timeseries_filepath);

% ========================================================================
% CALCULATE SMOOTH DAILY CLIMATOLOGY WITH SPLINE SMOOTHING
% ========================================================================

fprintf('Calculating smooth daily climatology using spline interpolation...\n');

% Initialize: smooth climatology [365, n_time, n_lat, n_lon]
smooth_climatology = zeros(365, n_time, n_lat, n_lon, 'single');

% ========================================================================
% STEP 1: Apply spline interpolation
% ========================================================================

fprintf('  Step 1: Applying spline interpolation...\n');

% Day-of-year corresponding to middle of each month (365-day calendar, no leap days)
month_days = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349];  % Middle of each month

% Create daily values 1:365
days_doy = 1:365;

% For each time step and spatial point, fit spline through 12 monthly points
for t = 1:n_time
    if mod(t, 2) == 0
        fprintf('    Time step %d/%d\n', t, n_time);
    end
    
    for i = 1:n_lat
        for j = 1:n_lon
            % Extract monthly climatology values for this location and time
            monthly_values = squeeze(monthly_climatology_raw(:, t, i, j));
            
            % Check if all values are valid (not NaN)
            if ~any(isnan(monthly_values))
                % Add periodic boundary conditions for smooth wrap-around
                % (December to January smooth transition)
                x_months = [month_days(12)-365, month_days, month_days(1)+365];  % Wrap around
                y_values = [monthly_values(12); monthly_values; monthly_values(1)];
                
                % Fit cubic spline with periodic boundary conditions
                pp = spline(x_months, y_values);
                
                % Evaluate spline at all days of year
                smooth_climatology(t, i, j, :) = ppval(pp, days_doy);
                
            else
                % If any monthly value is NaN, set entire day-of-year to NaN
                smooth_climatology(t, i, j, :) = NaN;
            end
        end
    end
end

fprintf('  ✓ Spline interpolation complete\n\n');

% ========================================================================
% CALCULATE ANOMALIES FROM SMOOTH CLIMATOLOGY
% ========================================================================

fprintf('Calculating anomalies from smooth climatology...\n');

anomalies_smooth = cell(file_count, 1);
for f = 1:file_count
    doy = day_of_year_all(f);
    % Proper dimension handling
    climatology_field = squeeze(smooth_climatology(doy, :, :, :));
    anomalies_smooth{f} = single(all_data{f}) - climatology_field;
end

fprintf('✓ Anomalies from smooth climatology calculated\n\n');

% ========================================================================
% CALCULATE SPATIAL STATISTICS FROM SMOOTH CLIMATOLOGY
% ========================================================================

fprintf('Calculating spatial statistics from smooth climatology...\n');

% Time series of spatially-averaged UI and anomalies
UI_timeseries_smooth = zeros(file_count, n_time);
smooth_climatology_timeseries = zeros(file_count, n_time);
anomaly_timeseries_smooth = zeros(file_count, n_time);

for f = 1:file_count
    for t = 1:n_time
        % Spatial mean from raw data (same as before)
        UI_valid = squeeze(all_data{f}(t, :, :));
        UI_valid = UI_valid(~isnan(UI_valid));
        
        if ~isempty(UI_valid)
            UI_timeseries_smooth(f, t) = mean(UI_valid);
        else
            UI_timeseries_smooth(f, t) = NaN;
        end
        
        % Smooth climatology
        doy = day_of_year_all(f);
        clim_valid = squeeze(smooth_climatology(doy, t, :, :));
        clim_valid = clim_valid(~isnan(clim_valid));
        if ~isempty(clim_valid)
            smooth_climatology_timeseries(f, t) = mean(clim_valid);
        else
            smooth_climatology_timeseries(f, t) = NaN;
        end
        
        % Anomaly
        anomaly_timeseries_smooth(f, t) = UI_timeseries_smooth(f, t) - smooth_climatology_timeseries(f, t);
    end
end

fprintf('✓ Spatial statistics from smooth climatology calculated\n\n');

% ========================================================================
% SAVE SMOOTH CLIMATOLOGY TO NETCDF
% ========================================================================

fprintf('Saving smooth climatology to NetCDF file...\n');

smooth_climatology_filepath = fullfile(output_path, sprintf('UI_ekman_smooth_daily_climatology_%s.nc', mask_name));
if isfile(smooth_climatology_filepath)
    delete(smooth_climatology_filepath);
end

ncid = netcdf.create(smooth_climatology_filepath, mode);

% Define dimensions
doy_dimid = netcdf.defDim(ncid, 'day_of_year', 365);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
lat_dimid = netcdf.defDim(ncid, 'latitude', n_lat);
lon_dimid = netcdf.defDim(ncid, 'longitude', n_lon);

% Define variables - NO permute, write as MATLAB native order [doy, time, lat, lon]
doy_varid = netcdf.defVar(ncid, 'day_of_year', 'NC_INT', doy_dimid);
time_varid = netcdf.defVar(ncid, 'time', 'NC_INT', time_dimid);
lat_varid = netcdf.defVar(ncid, 'latitude', 'NC_DOUBLE', lat_dimid);
lon_varid = netcdf.defVar(ncid, 'longitude', 'NC_DOUBLE', lon_dimid);
smooth_climatology_varid = netcdf.defVar(ncid, 'UI_ekman_smooth_climatology', 'NC_FLOAT', ...
    [doy_dimid, time_dimid, lat_dimid, lon_dimid]);

netcdf.endDef(ncid);

% Write data - NO permute
netcdf.putVar(ncid, doy_varid, int32(1:365));
netcdf.putVar(ncid, time_varid, int32(1:n_time));
netcdf.putVar(ncid, lat_varid, double(lat));
netcdf.putVar(ncid, lon_varid, double(lon));
netcdf.putVar(ncid, smooth_climatology_varid, single(smooth_climatology));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, smooth_climatology_varid, 'long_name', 'Smooth daily climatology of UI_Ekman via spline interpolation');
netcdf.putAtt(ncid, smooth_climatology_varid, 'units', 'm^2 s^-1');
netcdf.putAtt(ncid, smooth_climatology_varid, 'description', 'Cubic spline interpolation through 12 monthly climatologies. Note: MATLAB stores in column-major order, so when read externally dimensions are [time, lat, lon, doy]');
netcdf.putAtt(ncid, smooth_climatology_varid, 'interpolation_method', 'cubic_spline_with_periodic_boundaries');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Smooth Daily Climatology (Spline)');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'leap_days', 'removed');

netcdf.close(ncid);

fprintf('✓ Saved smooth climatology: %s\n', smooth_climatology_filepath);

% Save smooth climatology time series and anomalies
timeseries_smooth_filepath = fullfile(output_path, sprintf('UI_ekman_smooth_timeseries_statistics_%s.nc', mask_name));
if isfile(timeseries_smooth_filepath)
    delete(timeseries_smooth_filepath);
end

ncid = netcdf.create(timeseries_smooth_filepath, mode);

% Define dimensions
file_dimid = netcdf.defDim(ncid, 'file', file_count);
time_dimid = netcdf.defDim(ncid, 'time', n_time);
year_dimid = netcdf.defDim(ncid, 'year', n_years);

% Define variables - NO permute
ui_ts_varid = netcdf.defVar(ncid, 'UI_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
clim_ts_varid = netcdf.defVar(ncid, 'smooth_climatology_timeseries', 'NC_FLOAT', [file_dimid, time_dimid]);
anom_ts_varid = netcdf.defVar(ncid, 'anomaly_timeseries_smooth', 'NC_FLOAT', [file_dimid, time_dimid]);
year_varid = netcdf.defVar(ncid, 'years', 'NC_INT', year_dimid);

netcdf.endDef(ncid);

% Write data - NO permute
netcdf.putVar(ncid, ui_ts_varid, single(UI_timeseries_smooth));
netcdf.putVar(ncid, clim_ts_varid, single(smooth_climatology_timeseries));
netcdf.putVar(ncid, anom_ts_varid, single(anomaly_timeseries_smooth));
netcdf.putVar(ncid, year_varid, int32(years));

% Add attributes
netcdf.reDef(ncid);
netcdf.putAtt(ncid, anom_ts_varid, 'long_name', 'UI_Ekman anomalies from smooth spline climatology');
netcdf.putAtt(ncid, anom_ts_varid, 'description', 'Daily UI - smooth daily climatology');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'title', 'UI_Ekman Smooth Climatology Time Series and Anomalies');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'mask', mask_name);
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'period', '2005-2024');
netcdf.putAtt(ncid, netcdf.getConstant('NC_GLOBAL'), 'interpolation_method', 'cubic_spline_with_periodic_boundaries');

netcdf.close(ncid);

fprintf('✓ Saved smooth climatology time series: %s\n\n', timeseries_smooth_filepath);
% ========================================================================
% CREATE COMPREHENSIVE COMPARISON PLOTS
% ========================================================================

fprintf('Creating comparison plots...\n');

% Set figure visibility to off for all plots
set(0, 'DefaultFigureVisible', 'off');

% Plot 1: Comparison - Original daily vs smooth climatology
fig_comp1 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
climatology_mean_doy = squeeze(mean(climatology, [3, 4], 'omitnan'));
climatology_mean_doy_daily = mean(climatology_mean_doy, 2);
smooth_climatology_mean = squeeze(mean(smooth_climatology, [3, 4], 'omitnan'));
smooth_climatology_mean_daily = mean(smooth_climatology_mean, 2);

plot(1:365, climatology_mean_doy_daily, 'b-', 'LineWidth', 2, 'DisplayName', 'Original 365-day');
hold on
plot(1:365, smooth_climatology_mean_daily, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Smooth (spline)');
xlabel('Day of Year', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparison: Original vs Smooth Daily Climatology (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
month_days = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
set(gca, 'XTick', month_days, 'XTickLabel', month_labels);
xlim([1, 365]);

saveas(fig_comp1, fullfile(plot_path, 'UI_climatology_original_vs_smooth.png'), 'png');
saveas(fig_comp1, fullfile(plot_path, 'UI_climatology_original_vs_smooth.fig'), 'fig');
print(fig_comp1, fullfile(plot_path, 'UI_climatology_original_vs_smooth.jpg'), '-djpeg', '-r150');
close(fig_comp1);
fprintf('  ✓ Saved: UI_climatology_original_vs_smooth (png, fig, jpg)\n');

% Plot 2: Comparison of all three anomalies (daily, monthly, smooth)
fig_comp2 = figure('Position', [100 100 1400 600], 'Color', 'w', 'Visible', 'off');
anomalies_daily_mean = mean(anomaly_timeseries, 2);
anomalies_monthly_mean = mean(anomaly_timeseries_monthly, 2);
anomalies_smooth_mean = mean(anomaly_timeseries_smooth, 2);

plot(dates_datetime, anomalies_daily_mean, 'Color', [0.5, 0.5, 1], 'LineWidth', 1, 'DisplayName', 'From daily climatology');
hold on
plot(dates_datetime, anomalies_monthly_mean, 'Color', [1, 0.5, 0.5], 'LineWidth', 1.5, 'DisplayName', 'From monthly climatology');
plot(dates_datetime, anomalies_smooth_mean, 'g-', 'LineWidth', 1.5, 'DisplayName', 'From smooth climatology');
plot(dates_datetime([1, end]), [0, 0], 'k--', 'LineWidth', 1);

xlabel('Time (2005-2024)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} Anomaly (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Comparison: Anomalies from Different Climatologies (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';

saveas(fig_comp2, fullfile(plot_path, 'UI_anomalies_all_comparison.png'), 'png');
saveas(fig_comp2, fullfile(plot_path, 'UI_anomalies_all_comparison.fig'), 'fig');
print(fig_comp2, fullfile(plot_path, 'UI_anomalies_all_comparison.jpg'), '-djpeg', '-r150');
close(fig_comp2);
fprintf('  ✓ Saved: UI_anomalies_all_comparison (png, fig, jpg)\n');

% Plot 3: Monthly climatology (step function)
fig_monthly = figure('Position', [100 100 1200 600], 'Color', 'w', 'Visible', 'off');
monthly_climatology_mean = squeeze(mean(monthly_climatology_raw, [3, 4], 'omitnan'));
monthly_climatology_mean_daily = mean(monthly_climatology_mean, 2);

bar(1:12, monthly_climatology_mean_daily, 'FaceColor', [1, 0.7, 0.5], 'EdgeColor', 'k', 'LineWidth', 1.5);
hold on
plot(1:12, monthly_climatology_mean_daily, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);

xlabel('Month', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('UI_{Ekman} (m²/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Monthly Climatology of UI_{Ekman} (50km mask)', 'FontSize', 14, 'FontWeight', 'bold');
grid on
ax = gca;
ax.FontSize = 11;
ax.FontWeight = 'bold';
month_labels = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};
set(gca, 'XTick', 1:12, 'XTickLabel', month_labels);
xlim([0.5, 12.5]);

saveas(fig_monthly, fullfile(plot_path, 'UI_monthly_climatology.png'), 'png');
saveas(fig_monthly, fullfile(plot_path, 'UI_monthly_climatology.fig'), 'fig');
print(fig_monthly, fullfile(plot_path, 'UI_monthly_climatology.jpg'), '-djpeg', '-r150');
close(fig_monthly);
fprintf('  ✓ Saved: UI_monthly_climatology (png, fig, jpg)\n');

fprintf('✓ All comparison plots created\n\n');

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
fprintf('  Files processed: %d\n\n', file_count);

fprintf('Output Files:\n');
fprintf('  1. UI_ekman_daily_climatology_%s.nc\n', mask_name);
fprintf('  2. UI_ekman_timeseries_statistics_%s.nc\n', mask_name);
fprintf('  3. UI_ekman_monthly_climatology_%s.nc\n', mask_name);
fprintf('  4. UI_ekman_monthly_timeseries_statistics_%s.nc\n', mask_name);
fprintf('  5. UI_ekman_smooth_daily_climatology_%s.nc\n', mask_name);
fprintf('  6. UI_ekman_smooth_timeseries_statistics_%s.nc\n\n', mask_name);

fprintf('Dimension Order (in NetCDF files after MATLAB flipping):\n');
fprintf('  Climatology variables: [time, lat, lon, doy] or [time, lat, lon, month]\n');
fprintf('  Time series variables: [time, file] or [time, year]\n');

fprintf('Comparison Options:\n');
fprintf('  - Daily anomalies: deviations from daily climatology (synoptic scale)\n');
fprintf('  - Monthly anomalies: deviations from calendar month climatology (seasonal scale)\n');
fprintf('  - Smooth anomalies: deviations from spline-smoothed climatology (sub-seasonal scale)\n');
fprintf('========================================\n');