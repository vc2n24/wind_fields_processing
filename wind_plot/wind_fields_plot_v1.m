%Script created in date 11/11/2025 by Vanessa Cavucci
% plot CCMP data wind fields of the CUS, in 12 plots per year 2024
% only the monthly average is plotted

clear; clc; close all;

% Base directory (update path style for Windows)
base_path = '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/Y2024';
save_path = '/noc/users/vc2n24/wind_fields/plot';
% Define years and months

months = 1:12;

% Loop over MONTHS
%figure for all months 

figure('Position',[100 100 1200 900], 'Color','w' );

    for m = 1:length(months)
        
        Mstr = sprintf('M%02d', months(m));
        Mnumber = sprintf ('%02d', months(m));

        % Construct full NetCDF file path
        % Example filename: CCMP_Wind_Analysis_202412_monthly_mean_V03.1_L4.nc
        file_name = sprintf(['CCMP_Wind_Analysis_2024' Mnumber '_monthly_mean_V03.1_L4.nc']);
        file_path = fullfile(base_path, Mstr, file_name);

        % Check if file exists
        if ~isfile(file_path)
            warning('File not found: %s', file_path);
            continue;
        end

        fprintf('Processing %s\n', file_path);

        % Read variables (use actual variable names in your file)
        lon = ncread(file_path, 'longitude');
        lat = ncread(file_path, 'latitude');
        uwnd = ncread(file_path, 'u');
        vwnd = ncread(file_path, 'v');

        % Compute wind speed magnitude
        wspd = sqrt(uwnd.^2 + vwnd.^2);

        % Plot — using a subset for quiver (to reduce clutter)
        subplot(4, 3, m);
        contourf(lon, lat, wspd', 20, 'LineColor', 'none');
       hold on;
       

    skip = 5;
    quiver(lon(1:skip:end), lat(1:skip:end), ...
           uwnd(1:skip:end, 1:skip:end)', ...
           vwnd(1:skip:end, 1:skip:end)', ...
           'k', 'AutoScale', 'on', 'AutoScaleFactor', 1.2);
    hold off;

    % Label subplot
    title(sprintf('Month %02d', months(m)));
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar;

    end



% Add an overall title for the entire figure
sgtitle('2024 Monthly Mean Wind Speed');

% Save the combined figure
saveas(gcf, fullfile(save_path, 'WindMap_2024_AllMonths_v1.png'));

