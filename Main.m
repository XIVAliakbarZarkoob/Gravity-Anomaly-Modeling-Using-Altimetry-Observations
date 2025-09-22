%% Aliakbar Zarkoob, AKA "XIV"
%  Gmail: XIV.Aliakbar.Zarkoob@gmail.com
%  Telegram: @XIVAliakbar

clc, clear, close all, beep off, format long g
set(0,'defaultTextInterpreter','latex')

%#ok<*SAGROW>
%#ok<*AGROW>
%#ok<*MINV>

%% Load Data & Initializations

% [dot_file, dot_path] = uigetfile('Select DOT files','DOT', '*.nc','MultiSelect','on');
% dot_file = sort(dot_file); lat = []; lng = []; dot = [];
% for i = 1:length(dot_file)
%     data = rncdf([dot_path, char(dot_file(i))]);
%     lat = [lat; data.glat_00];
%     lng = [lng; data.glon_00]; 
%     dot = [dot; data.dot_33];
% end
% DOT = table(lat, lng, dot);
% 
% [ssh_file, ssh_path] = uigetfile('Select SSH files','SSH', '*.nc','MultiSelect','on');
% ssh_file = sort(ssh_file); lat = []; lng = []; ssh = [];
% for i = 1:length(ssh_file)
%     data = rncdf([ssh_path, char(ssh_file(i))]);
%     lat = [lat; data.glat_00];
%     lng = [lng; data.glon_00]; 
%     ssh = [ssh; data.ssh_40];
% end
% SSH = table(lat, lng, ssh);

load('Data.mat')
LAT_LIM = [round(min(SSH.lat)), round(max(SSH.lat))];
LNG_LIM = [round(min(SSH.lng)), round(max(SSH.lng))];
clear data dot lat lng ssh i

%% Calculating Geoid Height

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery, title('Boarder')
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(DOT,"lat","lng",'filled',ColorVariable='dot',SizeData=5)
title('Dynamic Ocean Topography ($DOT$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(SSH,"lat","lng",'filled',ColorVariable='ssh',SizeData=5)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)
title('Sea Surface Height ($SSH$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');

% geolimits(LAT_LIM, LNG_LIM)

% % 3D Scatter
% figure()
% subplot(1, 2, 1)
% scatter3(SSH.lng, SSH.lat, SSH.ssh, 10, SSH.ssh,'filled')
% title('Sea Surface Height', 'FontSize', 14)
% colormap('turbo'), axis equal tight, grid on
% hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
%     'FontSize',12,'Interpreter','latex');
% subplot(1, 2, 2)
% scatter3(DOT.lng, DOT.lat, DOT.dot, 10, DOT.dot,'filled')
% title('Dynamic Ocean Topography', 'FontSize', 14)
% colormap('turbo'), axis equal tight, grid on
% hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
%     'FontSize',12,'Interpreter','latex');

figure()
subplot(1, 2, 1)
scatter(SSH.lng, SSH.lat, 5, SSH.ssh,'filled')
title('Sea Surface Height ($SSH$)', 'FontSize', 14)
xlabel('Longitude'), ylabel('Latitude')
colormap('turbo'), axis equal tight, grid on
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');
subplot(1, 2, 2)
scatter(DOT.lng, DOT.lat, 5, DOT.dot,'filled')
title('Dynamic Ocean Topography ($DOT$)', 'FontSize', 14)
xlabel('Longitude'), ylabel('Latitude')
colormap('turbo'), axis equal tight, grid on
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');
sgtitle('Altimeter Data', 'FontSize', 16)

GRID_STEP = 0.25; % Degrees
lng = (LNG_LIM(1):GRID_STEP:LNG_LIM(2))';
lat = (LAT_LIM(1):GRID_STEP:LAT_LIM(2))';
[lng, lat] = meshgrid(lng, lat);
lng = lng(:); lat = lat(:);

% figure
% addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
% geobasemap usgsimagery
% geoplot(lat, lng, '.w', 'LineWidth', 2)
% title('Study Grid', 'FontSize', 14), colormap('turbo')

METHOD = 'natural';
dot = griddata(real(full(double(DOT.lng))), real(full(double(DOT.lat))), ...
    real(full(double(DOT.dot))), lng, lat, METHOD);
ssh = griddata(real(full(double(SSH.lng))), real(full(double(SSH.lat))), ...
    real(full(double(SSH.ssh))), lng, lat, METHOD);
N = ssh - dot;
GRID = table(lng, lat, dot, ssh, N);


figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(GRID,"lat","lng",'filled',ColorVariable='ssh',SizeData=10)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)
title('Interpolated Sea Surface Height ($SSH$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(GRID,"lat","lng",'filled',ColorVariable='dot',SizeData=10)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)
title('Interpolated Dynamic Ocean Topography ($DOT$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(GRID,"lat","lng",'filled',ColorVariable='N',SizeData=10)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)
title('Geoid Height ($N$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');

%% Calculating Gravity Anomaly

dla = deg2rad(GRID_STEP);
dlo = deg2rad(GRID_STEP);
average_G = 980665; % mGal
R = 6371000; % meters

global_lo = (0:GRID_STEP:360)';
global_la = (-90:GRID_STEP:90)';
[global_la, global_lo] = meshgrid(global_la, global_lo);
global_la = global_la(:); global_lo = global_lo(:);
global_N = egm96geoid(global_la, global_lo);

idx = ~isnan(GRID.N);
study_N = (GRID.N(idx));
study_la = (GRID.lat(idx));
study_lo = (GRID.lng(idx));
number_study = height(study_N);
study_Grav = zeros(number_study, 1);
pos = zeros(number_study, 1); % Index
for i = 1:number_study
     pos(i) = find(((global_lo == GRID.lng(i)) + (global_la == GRID.lat(i))) == 2);
end

for i = 1:number_study 
    Nq_Np = global_N-study_N(i);
    Z_phin = ( sind((study_la(i)-global_la)/2).^2 + ...
        sind((study_lo(i)-global_lo)/2).^2.*cosd(study_la(i)).*cosd(global_la) ).^(3/2);
    component_G = Nq_Np.*cosd(global_la)./Z_phin; 
    component_G(isinf(component_G)) = 0;
    component_G(pos(i)) = 0; 
    study_Grav(i) = sum(component_G);
end
study_Grav = average_G*(-1*study_N/R-dla*dlo*study_Grav/(16*pi*R)); 
result = table(study_lo, study_la, study_Grav); 
result.Properties.VariableNames = {'lng', 'lat', 'dg'};

figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(result,"lat","lng",'filled',ColorVariable='dg',SizeData=20)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 2)
title('Gravity Anomaly ($\Delta{g}$)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex');

eval_N = readtable('XGM2019_3f13b4a31dc0c7fe3e26487f32ded04803fdd213adbd7d2ed3a318f90a3754c8.dat', 'NumHeaderLines', 32);
eval_N.Properties.VariableNames = {'idx', 'lng', 'lat', 'h', 'dg'};
figure
addCustomBasemap("usgsimagery","usgsimagery.mbtiles")
geobasemap usgsimagery
hold on
geoscatter(eval_N,"lat","lng",'filled',ColorVariable='dg',SizeData=35)
geoplot([LAT_LIM(1),LAT_LIM(1),LAT_LIM(2),LAT_LIM(2),LAT_LIM(1),LAT_LIM(1)], ...
    [LNG_LIM(1),LNG_LIM(2),LNG_LIM(2),LNG_LIM(1),LNG_LIM(1),LNG_LIM(2)], '-w', 'LineWidth', 6)
title('Gravity Anomaly ($\Delta{g}$) (XGM2019)', 'FontSize', 14), colormap('turbo')
hhh = colorbar(); set(get(hhh,'ylabel'),'String','[meters]', ...
    'FontSize',12,'Interpreter','latex'); 

%% Extra

% fid = fopen('grid.dat','w');
% fprintf(fid,'%f  %f\n', [GRID.lng(1:2:end) GRID.lat(1:2:end)]');
% fclose(fid);

