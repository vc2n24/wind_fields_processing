Folder created by Vanessa Cavucci (v.cavucci@soton.ac.uk) in date: 11/11/2025
--------------------------------------------------------------------------------
wind_fields_plot:
- v1 (plot the wind field molthy averages of CCMP dataset, only run through months of one year, not ocean mask applied);
- v2 (plot the wind field molthy averages of CCMP dataset, only run through months of one year, matlab ocean mask applied);
- v3 (plot the wind field molthy averages of CCMP dataset, only run through months of one multiple year, matlab ocean mask applied, adds colormap label)
- v4 (v3 (plot the wind field molthy averages of CCMP dataset, only run through months of one multiple year, gbco ocean mask applied, adds colormap label, adds coastline 200m isobath, there is a stride applied to gebco files which make the map with less resolution but faster to plot)
--------------------------------------------------------------------
plot_rotation_tau_v1 : is a scirpt to plot the tau components from the ccmp files with new variables inside, but it did not work well cause of the different dimensions of the variables from the files.

plot_rotation_tau_v2 : is the script which plot in the correct one.