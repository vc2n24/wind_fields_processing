%script to check on the isobath 200m and the distance from the line between
%Cape Juby and Cape Blanc, offshore of 50km. Using the gebco data.
% File input gebco loaded with gebco function

clear; clc; close all;

save_path_plot = '/noc/mpoc/rpdmoc/users/vc2n24/wind_fields/plot/plot_mask';
if ~isfolder(save_path_plot)
    mkdir(save_path_plot);
end

% load data gebco (assumes open_gebco is on path and returns lon, lat, elevation)
[lon_sub, lat_sub, elev_sub] = open_gebco;

%% find the tangent to the coast between Cape Juby and Cape Blanc

cj_lat = 27.96;
cj_lon = -13.20;

cb_lat = 21.17;
cb_lon = -17.20;

% transform them in radians
rad = pi / 180;

cj_lat_rad = cj_lat * rad;
cj_lon_rad = cj_lon * rad;

cb_lat_rad = cb_lat * rad;
cb_lon_rad = cb_lon * rad;

% calculate dx and dy: dx (difference in lon corrected by cos(lat)), dy (difference in lat)
dy = (cj_lat_rad - cb_lat_rad); % delta latitude
dx = (cj_lon_rad - cb_lon_rad) * cos((cj_lat_rad + cb_lat_rad) / 2); % delta longitude

theta = atan2(dx, dy);

% convert 50 km to radians
R = 6371; % Earth radius in km
off_rad = 50 / R;

% linear interpolation between cj and cb points
N = 200;
line = linspace(0, 1, N);

lat_coast_line = cb_lat_rad + line .* (cj_lat_rad - cb_lat_rad);
lon_coast_line = cb_lon_rad + line .* (cj_lon_rad - cb_lon_rad);

% normalize the coastline vector (direction from Cape Blanc to Cape Juby)
len = sqrt(dx^2 + dy^2);
ux = dx / len;
uy = dy / len;

% normal pointing offshore (rotate 90 deg)
nx = -uy;
ny = ux;

% offshore coastline 50 km
lat_off_line = lat_coast_line + ny * off_rad;
lon_off_line = lon_coast_line + (nx * off_rad ./ cos(lat_coast_line));

% transform them back to degrees
lat_cos_deg = lat_coast_line / rad;
lon_cos_deg = lon_coast_line / rad;

lat_off_deg = lat_off_line / rad;
lon_off_deg = lon_off_line / rad;

% Plot the coastline and the offshore line, together with 200 m isobath
figure;
contourf(lon_sub, lat_sub, elev_sub, 8, 'LineColor', 'none');
colormap(gca, flipud(gray));
hold on;

% coastline and 200 m isobath from GEBCO
contour(lon_sub, lat_sub, elev_sub, [0 0], 'k', 'LineWidth', 1.2);      % 0 m contour = coast
contour(lon_sub, lat_sub, elev_sub, [-200 -200], 'b', 'LineWidth', 1);  % 200 m isobath

plot(lon_cos_deg, lat_cos_deg, 'k', 'LineWidth', 2);
plot(lon_off_deg, lat_off_deg, 'r', 'LineWidth', 2);

xlabel('Longitude', 'FontSize', 12);
ylabel('Latitude', 'FontSize', 12);

% save figure
saveas(gcf, fullfile(save_path_plot, 'ui_mask.png'));