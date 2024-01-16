%% Confrontiamo i vari modelli
clear all
clc

load("result_data_multivariato_meteo_NOX.mat")
load("result_data_multivariato_meteo_PM10.mat")

tot_t_stat = (mean(result_data_multivariato_meteo_PM10{1,9},2) + mean(X{1,9},2))/2;

%{
ORDINE NELLA STRUCT
beta_cv;
theta_z_cv;
v_z_cv;
sigma_eta_cv;
G_cv;
sigma_eps_cv;
diag_varcov_cv;
log_likelihood_cv;
t_stat;
%}

tot = zeros(6);

for c = 1:9
    tot = tot + result_data_multivariato_meteo_NOX{1,3}{1,c};
end

tot = tot./9;

d = sqrt(diag(tot).*eye(6));
R = inv(d)*tot*inv(d);  

load("result_data_pm10_univariate_selected.mat")
result_data_pm10_univariate_selected{1,7};

%{
    UNI_NOX 16.1382
    UNI_PM10 9.8716
    
    BIVAR_NOX 16.0290
    BIVAR_PM10 9.8983

    MULTI_NOX 16.19
    MULTI_PM10 9.9054

    % come vedi sono uguali prof
%}

% migliore NOX = bivariato
% migliore PM10 = univariato

%% Confronto dei modelli sulla deviazione standard media

load("result_krig_univariato_NOX.mat");
load("result_krig_univariato_PM10.mat");
load("result_krig_bivariato.mat");
load("result_krig_multivariato.mat");

% y_hat, diag_std_y_hat, coordinate, zk_s (prima NOX poi PM10)

std_hat_univariato_NOX = result_krig_univariato_NOX{2};
std_hat_univariato_PM10 = result_krig_univariato_PM10{2};

std_hat_bivariato_NOX = result_krig_bivariato{2};
std_hat_bivariato_PM10 = result_krig_bivariato{6};

std_hat_multivariato_NOX = result_krig_multivariato{2};
std_hat_multivariato_PM10 = result_krig_multivariato{6};

y_hat_univariato_NOX = result_krig_univariato_NOX{2};
y_hat_univariato_PM10 = result_krig_univariato_PM10{2};

y_hat_bivariato_NOX = result_krig_bivariato{2};
y_hat_bivariato_PM10 = result_krig_bivariato{6};

y_hat_multivariato_NOX = result_krig_multivariato{2};
y_hat_multivariato_PM10 = result_krig_multivariato{6};

mean(mean(mean(std_hat_univariato_NOX, 3),2))
mean(mean(mean(std_hat_univariato_PM10,3), 2))

mean(mean(mean(std_hat_bivariato_NOX,3),2))
mean(mean(mean(std_hat_bivariato_PM10,3),2))

mean(mean(mean(std_hat_multivariato_NOX,3),2))
mean(mean(mean(std_hat_multivariato_PM10,3),2))

%% Grafici nox

madrid = shaperead('madrid-districtsgeojson.shp');

figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_univariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');


figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_bivariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');


figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_multivariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');

%% grafico y_hat
figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_univariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');


figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_bivariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');


figure
colorbar
caxis([0 100])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_multivariato_NOX(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');

%% grafici pm10

figure
colorbar

hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_univariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');



figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_bivariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');



figure
colorbar
caxis([0 50])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), std_hat_multivariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');


%% grafico y_hat
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_univariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');



figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_bivariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');


figure
colorbar
caxis([0 50])
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_multivariato_PM10(:,:,9),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','o','MarkerEdgeColor','y');