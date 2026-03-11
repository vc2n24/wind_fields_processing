% Plotting netCDF files of wind fields with vectors, months in 12 subplots


clear; clc; close all;


% Base directory (update path style for Windows)
base_path = '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/';
save_path = '/noc/users/vc2n24/wind_fields/plot';
mkdir(save_path, 'v4');
save_path = fullfile(save_path, 'v4');

%read gebco data bathymetry
%read coordinates 
% Read coordinates

%dx
dx=5;


gebco_file = ['/noc/users/vc2n24/data/GEBCO_11_Nov_2025/GEBCO_11_Nov_2025/' ...
    'gebco_2025_n40.0_s8.0_w-40.0_e8.0.nc'];
lon_gebco = ncread(gebco_file, 'lon');
lat_gebco = ncread(gebco_file, 'lat');



%Define region of interest
lonlim = [-30 -5];
latlim = [20 30];


idxlon = find(lon_gebco >= lonlim(1) & lon_gebco <= lonlim(2));
idxlat = find(lat_gebco >= latlim(1) & lat_gebco <= latlim(2));


lon_sub = ncread(gebco_file,'lon', idxlon(1), length(idxlon));
lat_sub = ncread(gebco_file,'lat', idxlat(1), length(idxlat));
elev_sub = ncread(gebco_file,'elevation', [idxlon(1) idxlat(1)], [length(idxlon) length(idxlat)]);
elev_sub = elev_sub';   % transpose to [lat,lon]

lon_sub=lon_sub(1:dx:end);
lat_sub=lat_sub(1:dx:end);
elev_sub=elev_sub(1:dx:end,1:dx:end);

months = 1:12;
years_numbers= 2024:2024;


% Ensure numeric, ascending, and column vectors if needed
lonlim = sort(double(lonlim(:)))';
latlim = sort(double(latlim(:)))';


% Loop over MONTHS
%figure for all months 



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

   

        % Plot using a subset for quiver
        subplot(4, 3, m);
        contourf(lon_sub, lat_sub, elev_sub, 8, 'LineColor', 'none');
       % colormap(gca, flipud(gray));  % gray background
        hold on;

       % Plot wind speed on top
       contourf(lat, lon, wspd', 12, 'LineColor', 'none', 'LineStyle', 'none');
       
       

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
    cb = colorbar;
   
    colormap(cmocean('speed'));
    clim([0 12]);
    ylabel(cb, 'Wind speed (m s^{-1})', 'FontSize', 9);
   

    end


% Add an overall title for the entire figure
sgtitle([Ynumber ' Monthly Mean Wind Speed']);

% Save the combined figure

saveas(gcf, fullfile(save_path, ['WindMap_' num2str(Ynumber) '_AllMonths_v4.png']));

close all

end

