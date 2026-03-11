%Plot the wind stress along shore intensity on the background and the wind
%speed vectors of the monthly mean data of ccmp


clear; clc; close all;

%file_paths
base_path_products= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products';
base_path_ccmp = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/ccmp_data';
save_path = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/plot';
mkdir(save_path, 'plot_rotation_tau_v3');
save_path = fullfile(save_path, 'plot_rotation_tau_v3');


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
        file_name_ccmp = sprintf(['CCMP_Wind_Analysis_' Ynumber Mnumber '_monthly_mean_V03.1_L4.nc']);
        file_path_ccmp = fullfile(base_path_ccmp, Ystr, Mstr, file_name_ccmp);
        
        file_name_products = sprintf(['CCMP_products_monthly_' Ynumber Mnumber '_V03.1_L4.nc']);
        file_path_products = fullfile(base_path_products, Ystr, Mstr, file_name_products);
        % Check if file exists
        
        if ~isfile(file_path_ccmp)
            warning('File not found: %s', file_path_ccmp);
            continue;
        end

        fprintf('Processing %s\n', file_path_ccmp);

        % Read variables 
        lon = ncread(file_path_ccmp, 'longitude');
        lat = ncread(file_path_ccmp, 'latitude');
        lon = lon(:);
        lat = lat(:);
         
        [Lon, Lat] = meshgrid(lon, lat);

        uwnd = ncread(file_path_ccmp, 'u');
        vwnd = ncread(file_path_ccmp, 'v');
        wspd = ncread(file_path_ccmp, 'w');
        
        %transpose
        % to [lat x lon]
        uwnd = uwnd';   
        vwnd = vwnd';
        wspd = wspd';
        
        %read tau_along already in lat, lon
        tau_along = ncread(file_path_products, 'tau_along');
        tau_along = tau_along'; 

        %create meshgrid
        [Lon, Lat] = meshgrid(lon, lat);

        % Plot subset for quiver 
        subplot(4, 3, m);
        contourf(lon_sub, lat_sub, elev_sub, 8, 'LineColor', 'none');
         % colormap(gca, flipud(gray));  % gray background
         hold on;
         
        contourf(Lon, Lat, tau_along, 15, 'LineColor', 'none', 'LineStyle', 'none');
        %subsample and plot quiver

        colormap(cmocean('balance'));
        clim([-0.2 0.2]);

        skip = 10;
        Lonq = Lon(1:skip:end, 1:skip:end);
        Latq = Lat(1:skip:end, 1:skip:end);
        uq   = uwnd(1:skip:end, 1:skip:end);
        vq   = vwnd(1:skip:end, 1:skip:end);
        quiver(Lonq, Latq, uq, vq, 'k', 'AutoScale', 'on', 'AutoScaleFactor', 1.2);


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
ylabel(h, 'Wind Stress Along-shore (N m^{-2})')

% Add an overall title for the entire figure
sgtitle([Ynumber ' Monthly Mean Wind Stress Along-shore']);

% Save the combined figure

saveas(gcf, fullfile(save_path, ['tau_3' num2str(Ynumber) '_AllMonths.png']));

close all
end
