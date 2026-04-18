clear; clc; close all;

% ========================================================================
% PLOT UI_EKMAN AS SPATIAL MAPS - LOOP OVER ALL DATA
% ========================================================================

% Define paths
ui_results_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/UI_ekman';
plot_path = fullfile(ui_results_path, 'spatial_maps');

% Create output directory for plots
if ~isfolder(plot_path)
    mkdir(plot_path);
end

% ========================================================================
% PROCESSING PARAMETERS
% ========================================================================

days = 1:31;                    % All days in month
months = 1:12;                  % All months
years_numbers = 2005:2024;      % All years
mask_distances = {'50km', '150km', '300km'};  % All masks

% ========================================================================
% MAIN LOOP OVER FILES AND MASKS
% ========================================================================

total_files = length(years_numbers) * length(months) * length(days) * length(mask_distances);
file_count = 0;

for y = 1:length(years_numbers)
    
    year = years_numbers(y);
    Ystr = sprintf('Y%04d', year);
    
    for m = 1:length(months)
        
        month = months(m);
        Mstr = sprintf('M%02d', month);
        
        % Determine number of days in this month
        if month == 2
            if mod(year, 4) == 0 && (mod(year, 100) ~= 0 || mod(year, 400) == 0)
                n_days = 29;  % Leap year
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
            Dnumber = sprintf('%02d', day);
            
            for mask_idx = 1:length(mask_distances)
                
                mask_name = mask_distances{mask_idx};
                file_count = file_count + 1;
                
                % ============================================================
                % READ UI_EKMAN NETCDF FILE
                % ============================================================
                
                % Build filename safely
                year_str = num2str(year);
                month_str = sprintf('%02d', month);
                day_str = sprintf('%02d', day);
                input_filename = ['UI_ekman_' year_str month_str day_str '_' mask_name '.nc'];
                
                % Build full path
                mask_dir = ['UI_ekman_' mask_name];
                input_filepath = fullfile(ui_results_path, mask_dir, Ystr, Mstr, input_filename);
                
                % Check if file exists
                if ~isfile(input_filepath)
                    if mod(file_count, 500) == 0
                        fprintf('File %d/%d not found: %s%s%s_%s\n', file_count, total_files, year_str, month_str, day_str, mask_name);
                    end
                    continue;
                end
                
                if mod(file_count, 500) == 0 || file_count == 1
                    fprintf('Processing file %d/%d: %s%s%s_%s\n', file_count, total_files, year_str, month_str, day_str, mask_name);
                end
                
                try
                    % Read coordinates and data
                    lon = ncread(input_filepath, 'longitude');
                    lat = ncread(input_filepath, 'latitude');
                    UI_ekman = ncread(input_filepath, 'UI_ekman');  % [time, lat, lon]
                    
                    [n_time, n_lat, n_lon] = size(UI_ekman);
                    
                catch ME
                    fprintf('  Error reading %s: %s\n', input_filename, ME.message);
                    continue;
                end
                
                % ============================================================
                % CREATE COMPOSITE MAP (AVERAGE OF ALL TIME STEPS)
                % ============================================================
                
                % Calculate mean UI across time
                UI_mean_time = mean(UI_ekman, 1, 'omitnan');
                UI_mean_time = squeeze(UI_mean_time);  % [lat, lon]
                
                % Create figure
                fig = figure('Position', [100 100 1200 800], 'Color', 'w');
                
                % Create meshgrid for plotting
                [LON_GRID, LAT_GRID] = meshgrid(lon, lat);
                
                % Plot mean UI_Ekman with colormap
                h = pcolor(LON_GRID, LAT_GRID, UI_mean_time);
                set(h, 'EdgeColor', 'none');
                
                % Colormap settings
                colormap(jet);
                cbar = colorbar;
                cbar.Label.String = 'UI_{Ekman} (m^2 s^{-1})';
                cbar.Label.FontSize = 12;
                cbar.Label.FontWeight = 'bold';
                
                % Add contours
                hold on
                contour(LON_GRID, LAT_GRID, UI_mean_time, 10, 'k-', 'LineWidth', 0.5);
                
                % Formatting
                xlabel('Longitude (°E)', 'FontSize', 14, 'FontWeight', 'bold');
                ylabel('Latitude (°N)', 'FontSize', 14, 'FontWeight', 'bold');
                title(sprintf('Mean UI_{Ekman} - %s mask\n%04d-%02d-%02d', ...
                    mask_name, year, month, day), ...
                    'FontSize', 16, 'FontWeight', 'bold');
                
                set(gca, 'YDir', 'normal');
                grid on
                ax = gca;
                ax.FontSize = 12;
                ax.FontWeight = 'bold';
                
                % Statistics
                UI_valid = UI_mean_time(~isnan(UI_mean_time));
                if ~isempty(UI_valid)
                    mean_ui = mean(UI_valid);
                    min_ui = min(UI_valid);
                    max_ui = max(UI_valid);
                    std_ui = std(UI_valid);
                    
                    stats_text = sprintf('Mean: %.3e m²/s\nStd: %.3e m²/s\nMin: %.3e m²/s\nMax: %.3e m²/s', ...
                        mean_ui, std_ui, min_ui, max_ui);
                    annotation('textbox', [0.02 0.7 0.15 0.25], 'String', stats_text, ...
                        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10, ...
                        'FontWeight', 'bold');
                end
                
                % Save figure
                output_filename = ['UI_ekman_mean_' year_str month_str day_str '_' mask_name '.png'];
                output_filepath = fullfile(plot_path, output_filename);
                saveas(fig, output_filepath, 'png');
                
                close(fig);
                
                % ============================================================
                % CREATE STANDARD DEVIATION MAP
                % ============================================================
                
                % Calculate std UI across time
                UI_std_time = std(UI_ekman, 1, 'omitnan');
                UI_std_time = squeeze(UI_std_time);  % [lat, lon]
                
                % Create figure
                fig = figure('Position', [100 100 1200 800], 'Color', 'w');
                
                % Plot std
                h = pcolor(LON_GRID, LAT_GRID, UI_std_time);
                set(h, 'EdgeColor', 'none');
                
                colormap(hot);
                cbar = colorbar;
                cbar.Label.String = 'Std Dev UI_{Ekman} (m^2 s^{-1})';
                cbar.Label.FontSize = 12;
                cbar.Label.FontWeight = 'bold';
                
                hold on
                contour(LON_GRID, LAT_GRID, UI_std_time, 10, 'k-', 'LineWidth', 0.5);
                
                xlabel('Longitude (°E)', 'FontSize', 14, 'FontWeight', 'bold');
                ylabel('Latitude (°N)', 'FontSize', 14, 'FontWeight', 'bold');
                title(sprintf('Std Dev UI_{Ekman} - %s mask\n%04d-%02d-%02d', ...
                    mask_name, year, month, day), ...
                    'FontSize', 16, 'FontWeight', 'bold');
                
                set(gca, 'YDir', 'normal');
                grid on
                ax = gca;
                ax.FontSize = 12;
                ax.FontWeight = 'bold';
                
                % Statistics
                UI_std_valid = UI_std_time(~isnan(UI_std_time));
                if ~isempty(UI_std_valid)
                    mean_std = mean(UI_std_valid);
                    min_std = min(UI_std_valid);
                    max_std = max(UI_std_valid);
                    
                    stats_text = sprintf('Mean Std: %.3e m²/s\nMin Std: %.3e m²/s\nMax Std: %.3e m²/s', ...
                        mean_std, min_std, max_std);
                    annotation('textbox', [0.02 0.7 0.15 0.2], 'String', stats_text, ...
                        'BackgroundColor', 'white', 'EdgeColor', 'black', 'FontSize', 10, ...
                        'FontWeight', 'bold');
                end
                
                % Save figure
                output_filename = ['UI_ekman_std_' year_str month_str day_str '_' mask_name '.png'];
                output_filepath = fullfile(plot_path, output_filename);
                saveas(fig, output_filepath, 'png');
                
                close(fig);
                
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('✓ ALL SPATIAL MAPS CREATED!\n');
fprintf('Total files processed: %d\n', file_count);
fprintf('Saved to: %s\n', plot_path);
fprintf('========================================\n');