%elimate the file of monthly products from all the folders

path= '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/data_products';


months = 1:12;
years_numbers= 2005:2025;

% Loop over all files


for y= 1:length(years_numbers)
    
    Ystr= sprintf('Y%04d', years_numbers(y));
    Ynumber= sprintf('%04d', years_numbers(y));
    

    for m = 1:length(months)
        
        Mstr = sprintf('M%02d', months(m));
        Mnumber = sprintf ('%02d', months(m));


        file_name= sprintf('CCMP_products_montly_%s%s%s', Ynumber, Mnumber);
        full_file= fullfile(path, Ystr, Mstr, file_name);
        
        delete(full_file);



    end
end
