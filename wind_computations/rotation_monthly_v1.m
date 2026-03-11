%script to compute rotate vectors in direction along shore Cape Juby to
%Cape Blanc and cross shore; also wind stress is computed and variables
%addedd to the original netcdf files of CCMP from Vanessa Cavucci (vc2n24@soton.ac.uk)
% input data ccmp from '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/';
% output data would be same but with addedd variables (u_along, u_cross, taux, tauy, tau_along, tau_cross)



clear; clc; close all;

base_path= '/noc/users/vc2n24/wind_fields/subset/ccmp/v03.1/';

%% find the tangent to the coast between Cape Juby and Cape Blanc

cj_lat= 27.96;
cj_lon= -13.20;

cb_lat= 21.17;
cb_lon= -17.20;

%transform them in radians

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


%% upload u, v component and calculate cross shore, along shore winds

months = 1:08;
years_numbers= 2025:2025;

% Loop over all files


for y= 1:length(years_numbers)
    
    Ystr= sprintf('Y%04d', years_numbers(y));
    Ynumber= sprintf('%04d', years_numbers(y));
    

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

        % Read variables (use actual variable names in your file)
        uwnd = ncread(file_path, 'u');
        vwnd = ncread(file_path, 'v');
        wspd = ncread(file_path, 'w');

        %% rotate u, v vector with theta, considering the final outcome
        %should be along positive from Cape Juby to Cape Blanc, and can add
        %the cross shore winds as well positive offshore
        
        %vectors
        U = uwnd;
        V = vwnd;

        %scalars projections of wind vectors on the coast and offshore,
        %once transformed the wind vector rotated would be (u_along,
        %u_cross)

        u_along= -(U .* sin(theta) + V .* cos(theta));
        u_cross= - (U.* cos(theta) - V .*sin(theta));

        %calculate w to double check
        w_rotated = sqrt(u_along.^2 + u_cross.^2);
        %Check on calculations with ws which should results equal from both
        diff_w= abs(w_rotated - wspd); %plot is overlapping is well


        %% calculate wind stress tau from function windstress
        

        [taux, tauy] = windstress(uwnd,vwnd);
        
        % rotate stress using the same convention as the wind
        tau_along = - (taux .* sin(theta) + tauy .* cos(theta) );
        tau_cross = - (taux .* cos(theta) - tauy .* sin(theta) );
        
        %need to add the variables u_alon; u_cross; taux; tauy; tau_along;
        %tau_cross, ekman transport and specify in the attributes what they stay for etc..
        %nccreate ; ncwrite; ncwriteatt (size 3D: lon, lat, time)


       
        u_along = permute(u_along, [3 2 1]);
        u_cross = permute(u_cross, [3 2 1]);
        taux = permute(taux, [3 2 1]);
        tauy = permute(tauy, [3 2 1]);
        tau_along = permute(tau_along,[3 2 1]);
        tau_cross = permute(tau_cross,[3 2 1]);
        
          
        nccreate(file_path, 'u_along', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));
        nccreate(file_path, 'u_cross', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));
        nccreate(file_path, 'taux', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));
        nccreate(file_path, 'tauy', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));
        nccreate(file_path, 'tau_along', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));
        nccreate(file_path, 'tau_cross', 'Dimensions', {'time', 4 , 'latitude', 128, 'longitude', 192}, 'Datatype','single', 'FillValue', single(-9999));

        
        ncwrite(file_path, 'u_along', single(u_along));
        ncwrite(file_path, 'u_cross', single(u_cross));
        ncwrite(file_path, 'taux', single(taux));
        ncwrite(file_path, 'tauy', single(tauy));
        ncwrite(file_path, 'tau_along', single(tau_along));
        ncwrite(file_path, 'tau_cross', single(tau_cross));

        ncwriteatt(file_path, 'u_along', 'standard_name', 'wind_speed_along_coast');
        ncwriteatt(file_path, 'u_along', 'long_name', 'Along-shore wind component');
        ncwriteatt(file_path, 'u_along', 'units', 'm s-1');
        ncwriteatt(file_path, 'u_along', 'height', '10 meters above sea-surface');
        ncwriteatt(file_path, 'u_along', 'method', ['Computed from rotated wind ' ...
            'using coastline angle theta: ', ...
     'dy = cj_lat_rad - cb_lat_rad; ', ...
     'dx = (cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad+cb_lat_rad)/2); ', ...
     'theta = atan2(dx,dy); ', ...
     'u_along = -(U .* sin(theta) + V .* cos(theta)).']);
        ncwriteatt(file_path, 'u_along', 'comment', 'Positive from Cape Juby to Cape Blanc');
          


        ncwriteatt(file_path, 'u_cross', 'standard_name', 'wind_speed_across_coast');
        ncwriteatt(file_path, 'u_cross', 'long_name', 'Cross-shore wind component');
        ncwriteatt(file_path, 'u_cross', 'units', 'm s-1');
        ncwriteatt(file_path, 'u_cross', 'height', '10 meters above sea-surface');
        ncwriteatt(file_path, 'u_cross', 'method', ['Computed from rotated wind ' ...
            'using coastline angle theta: ', ...
     'dy = cj_lat_rad - cb_lat_rad; ', ...
     'dx = (cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad+cb_lat_rad)/2); ', ...
     'theta = atan2(dx,dy); ', ...
     'u_cross = -( U*cos(theta) - V*sin(theta) ).']);
        ncwriteatt(file_path, 'u_cross', 'comment', 'Positive from on-shore to offshore');

        ncwriteatt(file_path, 'taux', 'standard_name', 'eastward_wind_stress');
        ncwriteatt(file_path, 'taux', 'long_name', 'Zonal (eastward) wind stress');
        ncwriteatt(file_path, 'taux', 'units', 'N m-2');
        ncwriteatt(file_path, 'taux', 'height', '10 m above sea surface');
        ncwriteatt(file_path, 'taux', 'method', ...
            ['Computed from uwnd using taux = rho_air * C_D * |W| * uwnd. ', ...
     'Drag coefficient from Kara et al. (2005).']);
        

        ncwriteatt(file_path, 'tauy', 'standard_name', 'northward_wind_stress');
        ncwriteatt(file_path, 'tauy', 'long_name',     'Meridional (northward) wind stress');
        ncwriteatt(file_path, 'tauy', 'units',         'N m-2');
        ncwriteatt(file_path, 'tauy', 'height',        '10 m above sea surface');
        ncwriteatt(file_path, 'tauy', 'method', ...
            ['Computed from vwnd using tauy = rho_air * C_D * |W| * vwnd. ', ...
     'Drag coefficient from Kara et al. (2005).']);

        ncwriteatt(file_path, 'tau_along', 'standard_name', 'wind_stress_along_coast');
        ncwriteatt(file_path, 'tau_along', 'long_name',     'Along-shore wind stress component');
        ncwriteatt(file_path, 'tau_along', 'units',         'N m-2');
        ncwriteatt(file_path, 'tau_along', 'height',        '10 m above sea surface');
        ncwriteatt(file_path, 'tau_along', 'method', ...
            ['Computed by rotating (taux, tauy) using coastline angle theta: ', ...
     'tau_along = -( taux*sin(theta) + tauy*cos(theta) ). ', ...
     'Positive from Cape Juby toward Cape Blanc.']);
        ncwriteatt(file_path, 'tau_along', 'comment', ['Positive along-shore ' ...
            '= direction Cape Juby to Cape Blanc.']);

       ncwriteatt(file_path, 'tau_cross', 'standard_name', 'wind_stress_across_coast');
       ncwriteatt(file_path, 'tau_cross', 'long_name',     'Cross-shore (offshore) wind stress component');
       ncwriteatt(file_path, 'tau_cross', 'units',         'N m-2');
       ncwriteatt(file_path, 'tau_cross', 'height',        '10 m above sea surface');
       ncwriteatt(file_path, 'tau_cross', 'method', ...
           ['Computed by rotating (taux, tauy) using coastline angle theta: ', ...
     'tau_cross = -( taux*cos(theta) - tauy*sin(theta) ). ', ...
     'Positive from east to west (offshore).']);
       ncwriteatt(file_path, 'tau_cross', 'comment', ...
           'Positive cross-shore = offshore (east to west).');
        

    end
end
