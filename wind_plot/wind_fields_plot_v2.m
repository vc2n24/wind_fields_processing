% Plotting netCDF files of wind fields with vectors, months in 12 subplots


clear; clc; close all;

% Base directory (update path style for Windows)
base_path = '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/Y2024';
save_path = '/noc/users/vc2n24/wind_fields/plot';
% Define years and months

months = 1:12;

%define region of interest

lonlim = [-30 -5]
latlim = [ 20 30]

% Ensure numeric, ascending, and column vectors if needed
lonlim = sort(double(lonlim(:)))';
latlim = sort(double(latlim(:)))';


coast = load('coastlines');

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
        wspd = ncread(file_path, 'w');

   

        % Plot — using a subset for quiver (to reduce clutter)
        subplot(4, 3, m);
        contourf(lon, lat, wspd', 40, 'LineColor', 'none');
       hold on;
       

    skip = 10;
    quiver(lon(1:skip:end), lat(1:skip:end), ...
           uwnd(1:skip:end, 1:skip:end)', ...
           vwnd(1:skip:end, 1:skip:end)', ...
           'k', 'AutoScale', 'on', 'AutoScaleFactor', 1.2);

        % Overlay land mask (coastlines)
    plot(coast.coastlon, coast.coastlat, 'k', 'LineWidth', 1);

        % Zoom to region
    xlim(lonlim);
    ylim(latlim);
    set(gca,'YDir','normal');  % keep north upwards

    % Label subplot
    title(sprintf('Month %02d', months(m)));
    xlabel('Longitude');
    ylabel('Latitude');
    colorbar;
   
    colormap(cmocean('speed'));
    caxis([0 12]);
    cb = colorbar;
    ylabel(cb, 'Wind speed (m s^{-1})', 'FontSize', 9);
   

    end


% Add an overall title for the entire figure
sgtitle('2024 Monthly Mean Wind Speed');

% Save the combined figure
saveas(gcf, fullfile(save_path, 'WindMap_2024_AllMonths_v2.png'));

