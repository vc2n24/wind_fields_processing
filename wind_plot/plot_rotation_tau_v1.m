%Plot the wind stress along shore intensity on the background and the wind
%speed vectors of the monthly mean data of ccmp


clear; clc; close all;

%file_paths
base_path= '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/';
save_path = '/noc/users/vc2n24/wind_fields/plot';
mkdir(save_path, 'wind_stress');
save_path = fullfile(save_path, 'wind_stress');


%Define region of interest
lonlim = [-30 -5];
latlim = [20 30];

%load gebco subset with lon_sub; lat_sub and elev_sub
[lon_sub, lat_sub, elev_sub] = open_gebco;

%% upload wind fields data, wind stress and components

months = 1:12;
years_numbers= 2024:2024;

% Loop over all files


for y= 1:length(years_numbers)
    
    Ystr= sprintf('Y%04d', years_numbers(y));
    Ynumber= sprintf('%04d', years_numbers(y));
    figure('Position',[100 100 1200 900], 'Color','w' );

    for m = 1:length(months)
        
        Mstr = sprintf('M%02d', months(m));
        Mnumber = sprintf ('%02d', months(m));

            
        % Construct full NetCDF file path
        % Example filename: CCMP_Wind_Analysis_202412_monthly_mean_V03.1_L4.nc
        file_name = sprintf(['CCMP_Wind_Analysis_' Ynumber Mnumber '_monthly_mean_V03.1_L4.nc']);
        file_path = fullfile(base_path, Ystr, Mstr, file_name);

        % Check if file exists
        
        if ~isfile(file_path)
            warning('File not found: %s', file_path);
            continue;
        end

        fprintf('Processing %s\n', file_path);

        % Read variables 
        lon = ncread(file_path, 'longitude');
        lat = ncread(file_path, 'latitude');
        uwnd = ncread(file_path, 'u');
        vwnd = ncread(file_path, 'v');
        wspd = ncread(file_path, 'w');
        u_along = ncread(file_path, 'u_along');
        u_cross = ncread(file_path, 'u_cross');
        taux = ncread(file_path, 'taux');
        tauy = ncread(file_path, 'tauy');
        tau_along = ncread(file_path, 'tau_along');
        tau_cross = ncread(file_path, 'tau_cross');
        

        % Plot — using a subset for quiver 
        subplot(4, 3, m);
        contourf(lon_sub, lat_sub, elev_sub, 8, 'LineColor', 'none');
       % colormap(gca, flipud(gray));  % gray background
        hold on;

        % Plot wind stress along-shore on top
        %need to extract the taux 2D variable cause right now is
        %timexlatxlon but contourf is only 2D lat lon
        %problem! the direction lonxlatxtime are at the reverse! to change

        contourf(Lon, Lat, taux', 12, 'LineColor', 'none', 'LineStyle', 'none');
        skip = 10;
        quiver(lon(1:skip:end), lat(1:skip:end), ...
           uwnd(1:skip:end, 1:skip:end)', ...
           vwnd(1:skip:end, 1:skip:end)', ...
           'k', 'AutoScale', 'on', 'AutoScaleFactor', 1.2);

%coastline and 200 m isobath from GEBCO
contour(lon_sub, lat_sub, elev_sub, [0 0], 'k', 'LineWidth', 1.2);      % 0 m contour = coast
contour(lon_sub, lat_sub, elev_sub, [-200 -200], 'b', 'LineWidth', 1);  % 200 m isobath

    % Zoom to region
    xlim(lonlim);
    ylim(latlim);
    set(gca,'YDir','normal');  % keep north upwards

    % Label subplot
    title(sprintf('Month %02d', months(m)));
    xlabel('Longitude');
    ylabel('Latitude');

        
    end

h = colorbar('Position', [0.92 0.11 0.02 0.78]);  
colormap(cmocean('balance'));
clim([37 40]);
ylabel(h, 'Wind Stress Along-shore (N m^{-2})')

% Add an overall title for the entire figure
sgtitle([Ynumber ' Monthly Mean Wind Stress Along-shore']);

% Save the combined figure

saveas(gcf, fullfile(save_path, ['WindStress_' num2str(Ynumber) '_AllMonths.png']));

close all
end
