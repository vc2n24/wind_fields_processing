%rename varibale v_cross as u_cross from file : CCMP_Wind_Analysis_20240101_V03.1_L4.nc

file_path= 'Z:\users\vc2n24\wind_fields\subset\ccmp\v03.1\Y2024\M01';
name= 'CCMP_Wind_Analysis_20240101_V03.1_L4.nc';

full_path= fullfile(file_path, name);

%open file

cd(file_path)

%open netcdf file
necid = netcdf.open(name, 'NC_WRITE');

%put file in define mode
netcdf.reDef(necid)

% Get name of variable to change
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(necid,8);

% Rename the variable, using a capital letter to start the name.
netcdf.renameVar(necid,8,'u_cross')

% Verify that the name of the variable changed.
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(necid,8);

%close netcdf properly to save changes
netcdf.close(necid)