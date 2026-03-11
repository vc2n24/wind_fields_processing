%Eliminate varibales, and attributes from NetCDF files.
%input file:ccmp_wind_fields with wrong latxlonxtime orientation of taux,
%tauy, u_along, u_cross and tau_along, tau_cross. Need to just remove those
%variables and rename the file with same name.



clear; clc; close all;

base_path= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/ccmp_data/';


% Loop over all files
%% 
days = 1:31;
months = 1:8;
years_numbers= 2025:2025;


for y= 1:length(years_numbers)
    
    Ystr= sprintf('Y%04d', years_numbers(y));
    Ynumber= sprintf('%04d', years_numbers(y));
    figure('Position',[100 100 1200 900], 'Color','w' );

    for m = 1:length(months)
        
        Mstr = sprintf('M%02d', months(m));
        Mnumber = sprintf ('%02d', months(m));

        %make loop for days 

        for d = 1:length(days)
        
            Dnumber = sprintf ('%02d', days(d));
            
        % Construct full NetCDF file path
        % Example filename: CCMP_Wind_Analysis_202412_monthly_mean_V03.1_L4.nc
        file_name = sprintf(['CCMP_Wind_Analysis_' Ynumber Mnumber Dnumber '_V03.1_L4.nc']);
        file_path = fullfile(base_path, Ystr, Mstr, file_name);

        % Check if file exists
        
        if ~isfile(file_path)
            warning('File not found: %s', file_path);
            continue;
        end

        fprintf('Processing %s\n', file_path);

        %prepare a temporaneous file to use in the nco function
        
        temp_file = [file_path, '.tmp'];
        
        %make a list of variable to delete, because system want accept to 
        %make the string name inside the arguments of the function need to
        %create the argument before cmd
        vars= {'u_along' ,'u_cross', 'taux', 'tauy', 'tau_along', 'tau_cross'};
        vararg= strjoin(vars, ',');
        
        %-O overwrite authomatically
        cmd = sprintf('ncks -O -x -v %s "%s" "%s"', vararg, file_path, temp_file);
        system(cmd)

        %overwrite the original file safely, nco does not like to overwrite
        movefile(temp_file, file_path, 'f');

        end
    end
end
