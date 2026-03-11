%script to compute rotate vectors in direction along shore Cape Juby to
%Cape Blanc and cross shore modified from rotation_daily_v1.m 
% ; also wind stress is computed and variables
%addedd to the original netcdf files of CCMP from Vanessa Cavucci (vc2n24@soton.ac.uk)
% input data ccmp from 'Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\ccmp_data';
% output data would be organized in fodlers by year and months
% with variables (u_along, u_cross, taux, tauy, tau_along, tau_cross)
% calculated form the ccmp.
% Modified script to save new variables/attributes to a *new* NetCDF file
% Files will be saved to Z:\mpoc\rpdmoc\users\vc2n24\wind_fields\data_products
% Organized in folders divided by year (Y2005-Y2025), month (M01-M12)
% File naming: CCMP_products_'YnumberMnumberDnumber'_V03.1_L4.nc
% General attributes added for provenance

clear; clc; close all;

base_path= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/ccmp_data';
save_base = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products';

%Two points to draw the coastal line from Cape Juby to Cape Blanc
cj_lat= 27.96; 
cj_lon= -13.20; 
cb_lat= 21.17; 
cb_lon= -17.20;

%transform them in radiants

rad= pi / 180; 

cj_lat_rad= cj_lat * rad; 
cj_lon_rad= cj_lon * rad;
cb_lat_rad= cb_lat * rad; 
cb_lon_rad= cb_lon * rad;

% calculate dx and dy being: dx (the difference between lon, and dy 
% the difference between lat of the north point cape juby (2) and south
% point cape blanc (1). 

dy= (cj_lat_rad - cb_lat_rad);
dx= (cj_lon_rad - cb_lon_rad) * cos ((cj_lat_rad+cb_lat_rad)/2);

theta= atan2(dx,dy);

% upload u, v component and calculate cross shore, along shore winds

days = 1:31;
months = 1:8;
years_numbers= 2025:2025;

% Loop over all files

for y= 1:length(years_numbers)

    Ystr= sprintf('Y%04d', years_numbers(y));
    Ynumber= sprintf('%04d', years_numbers(y));

    for m = 1:length(months)

        Mstr = sprintf('M%02d', months(m));
        Mnumber = sprintf ('%02d', months(m));

        for d = 1:length(days)

            Dnumber = sprintf ('%02d', days(d));
            file_name = sprintf(['CCMP_Wind_Analysis_' Ynumber Mnumber Dnumber '_V03.1_L4.nc']);
            file_path = fullfile(base_path, Ystr, Mstr, file_name);

            if ~isfile(file_path)
                warning('File not found: %s', file_path);
                continue;
            end

            fprintf('Processing %s\n', file_path);

            % Read variables
            uwnd = ncread(file_path, 'uwnd');
            vwnd = ncread(file_path, 'vwnd');
            wspd = ncread(file_path, 'ws');

            lat = ncread(file_path,'latitude');
            lon = ncread(file_path,'longitude');
            
            lat = double(lat(:));
            lon = double(lon(:));

            %vectors
            U = uwnd; 
            V = vwnd;
            
            % Permute to get the same dimensions as the original ccmp dataset
            %right now dimensions are lon, lat, time but I want them in lat,
            %lon , time
             
            uwnd = permute(uwnd, [2 1 3]);
            vwnd = permute(vwnd, [2 1 3]);
            wspd = permute(wspd, [2 1 3]);

            %scalars projections of wind vectors on the coast and offshore,
            %once transformed the wind vector rotated would be (u_along,
            %u_cross)

            u_along= -(U .* sin(theta) + V .* cos(theta));
            u_cross= - (U.* cos(theta) - V .*sin(theta));
            w_rotated = sqrt(u_along.^2 + u_cross.^2);
             
            % calculate wind stress tau from function windstress which is
            % using CD from Sylla and Santos
            [taux, tauy] = windstress(uwnd,vwnd);
            %rotate
            tau_along = - (taux .* sin(theta) + tauy .* cos(theta) );
            tau_cross = - (taux .* cos(theta) - tauy .* sin(theta) );

            %need to add the variables u_alon; u_cross; taux; tauy; tau_along;
            %tau_cross, ekman transport and specify in the attributes what they stay for etc..
            %nccreate ; ncwrite; ncwriteatt (size 3D: lat, lon, time)

            % Prepare output file path 
            out_dir = fullfile(save_base, Ystr, Mstr);
            if ~isfolder(out_dir)
                mkdir(out_dir);
            end
            out_file = sprintf('CCMP_products_%s%s%s_V03.1_L4.nc', Ynumber, Mnumber, Dnumber);
            out_path = fullfile(out_dir, out_file);

            % Create variables in the new file 
            nccreate(out_path, 'u_along', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'u_cross', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'taux', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'tauy', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'tau_along', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'tau_cross', 'Dimensions', {'latitude', 128, 'longitude', 192, 'time', 4}, 'Datatype','single', 'FillValue', single(-9999));
            nccreate(out_path, 'latitude',  'Dimensions', {'latitude', numel(lat)},  'Datatype','double');
            nccreate(out_path, 'longitude', 'Dimensions', {'longitude', numel(lon)}, 'Datatype','double');

            ncwrite(out_path, 'u_along', single(u_along));
            ncwrite(out_path, 'u_cross', single(u_cross));
            ncwrite(out_path, 'taux', single(taux));
            ncwrite(out_path, 'tauy', single(tauy));
            ncwrite(out_path, 'tau_along', single(tau_along));
            ncwrite(out_path, 'tau_cross', single(tau_cross));
            ncwrite(out_path,  'latitude',  lat);
            ncwrite(out_path,  'longitude', lon);

            %Variable attributes 
            ncwriteatt(out_path, 'u_along', 'standard_name', 'wind_speed_along_coast');
            ncwriteatt(out_path, 'u_along', 'long_name', 'Along-shore wind component');
            ncwriteatt(out_path, 'u_along', 'units', 'm s-1');
            ncwriteatt(out_path, 'u_along', 'height', '10 meters above sea-surface');
            ncwriteatt(out_path, 'u_along', 'method', ['Computed from rotated wind ' ...
                'using coastline angle theta: dy = cj_lat_rad - cb_lat_rad; dx = ' ...
                '(cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad+cb_lat_rad)/2); theta = atan2(dx,dy); ' ...
                'u_along = -(U .* sin(theta) + V .* cos(theta)).']);
            ncwriteatt(out_path, 'u_along', 'comment', 'Positive from Cape Juby to Cape Blanc');
            
            ncwriteatt(out_path, 'u_cross', 'standard_name', 'wind_speed_across_coast');
            ncwriteatt(out_path, 'u_cross', 'long_name', 'Cross-shore wind component');
            ncwriteatt(out_path, 'u_cross', 'units', 'm s-1');
            ncwriteatt(out_path, 'u_cross', 'height', '10 meters above sea-surface');
            ncwriteatt(out_path, 'u_cross', 'method', 'Computed from rotated wind using coastline angle theta');
            ncwriteatt(out_path, 'u_cross', 'comment', 'Positive from on-shore to offshore');
            
            ncwriteatt(out_path, 'taux', 'standard_name', 'eastward_wind_stress');
            ncwriteatt(out_path, 'taux', 'long_name', 'Zonal (eastward) wind stress');
            ncwriteatt(out_path, 'taux', 'units', 'N m-2');
            ncwriteatt(out_path, 'taux', 'height', '10 m above sea surface');
            ncwriteatt(out_path, 'taux', 'method', ['Computed from uwnd using taux = rho_air * C_D * ' ...
                '|W| * uwnd. Drag coefficient from Kara et al. (2005).']);

            ncwriteatt(out_path, 'tauy', 'standard_name', 'northward_wind_stress');
            ncwriteatt(out_path, 'tauy', 'long_name', 'Meridional (northward) wind stress');
            ncwriteatt(out_path, 'tauy', 'units', 'N m-2');
            ncwriteatt(out_path, 'tauy', 'height', '10 m above sea surface');
            ncwriteatt(out_path, 'tauy', 'method', ['Computed from vwnd using tauy = rho_air * C_D * ' ...
                '|W| * vwnd. Drag coefficient from Kara et al. (2005).']);

            ncwriteatt(out_path, 'tau_along', 'standard_name', 'wind_stress_along_coast');
            ncwriteatt(out_path, 'tau_along', 'long_name', 'Along-shore wind stress component');
            ncwriteatt(out_path, 'tau_along', 'units', 'N m-2');
            ncwriteatt(out_path, 'tau_along', 'height', '10 m above sea surface');
            ncwriteatt(out_path, 'tau_along', 'method', 'Computed by rotating (taux, tauy) using coastline angle theta');
            ncwriteatt(out_path, 'tau_along', 'comment', 'Positive along-shore = direction Cape Juby to Cape Blanc.');

            ncwriteatt(out_path, 'tau_cross', 'standard_name', 'wind_stress_across_coast');
            ncwriteatt(out_path, 'tau_cross', 'long_name', 'Cross-shore (offshore) wind stress component');
            ncwriteatt(out_path, 'tau_cross', 'units', 'N m-2');
            ncwriteatt(out_path, 'tau_cross', 'height', '10 m above sea surface');
            ncwriteatt(out_path, 'tau_cross', 'method', 'Computed by rotating (taux, tauy) using coastline angle theta');
            ncwriteatt(out_path, 'tau_cross', 'comment', 'Positive cross-shore = offshore (east to west).');
            
            ncwriteatt(out_path, 'latitude', 'units', 'degrees_north');
            ncwriteatt(out_path, 'longitude','units', 'degrees_east');
            
            %File level attributes for provenance 
            ncwriteatt(out_path, '/', 'source_file', file_path);
            ncwriteatt(out_path, '/', 'generation_date', string(datetime("now"),'yyyy-MM-dd HH:mm:ss'));
            ncwriteatt(out_path, '/', 'description', ...
                'Wind-derived products (rotated components and wind stress) created from CCMP input files.');
            ncwriteatt(out_path, '/', 'script_used', 'ccmp_daily_products_script.m');
        end
    end
end

