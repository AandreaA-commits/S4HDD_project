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

nanmean(nanmean(nanmean(std_hat_univariato_NOX(val), 3),2))
nanmean(nanmean(nanmean(std_hat_univariato_PM10(val),3), 2))

nanmean(nanmean(nanmean(std_hat_bivariato_NOX(val),3),2))
nanmean(nanmean(nanmean(std_hat_bivariato_PM10(val),3),2))

nanmean(nanmean(nanmean(std_hat_multivariato_NOX(val),3),2))
nanmean(nanmean(nanmean(std_hat_multivariato_PM10(val),3),2))

%% Grafico delle stazioni di PM10 e NOX
madrid = shaperead('contorno_madrid.shp');



figure
grid on
geoscatter(PM10_coord(:,1), PM10_coord(:,2), 'filled', 'b');
hold on
for i=1:length(madrid)
    geoplot(madrid(i).Y, madrid(i).X, "k-");
end


NOX_coord = result_krig_univariato_NOX{3};

figure
grid on
geoscatter(NOX_coord(:,1), NOX_coord(:,2), 'filled', 'b');
hold on
for i=1:length(madrid)
    geoplot(madrid(i).Y, madrid(i).X, "k-");
end

%% stazioni del meteo (non eseguire)

all_lat =[TEMPERATURA_lat; UMIDITA_lat; VEL_VENTO_lat; PRESSIONE_lat];
all_long = [TEMPERATURA_long; UMIDITA_long; VEL_VENTO_long; PRESSIONE_long];

all = [all_lat, all_long];

c = unique(all,"rows", 'stable');
figure
grid on
geoscatter(c(:,1), c(:,2), 'filled', 'b');
hold on
for i=1:length(madrid)
    geoplot(madrid(i).Y, madrid(i).X, "k-");
end

%% Cerco il giorno più inquinato
[val_max, idx_max] = max(nanmean(dati_NOX,1));
[val_min, idx_min] = min(nanmean(dati_NOX,1));

[val_max2, idx_max2] = max(nanmean(dati_PM25,1));
[val_min2, idx_min2] = min(nanmean(dati_PM25,1));

%% plot del giorno più inquinato

load("result_krig_bivariato.mat")

PM10_coord = result_krig_univariato_PM10{3};

val = inpolygon(a(:,1), a(:,2), madrid.Y, madrid.X);
reshape(val, 56, 56);
t = result_krig_bivariato{1,1}(:,:,idx_max);
t(val == 0) = NaN;

for i=1:56
    for j=1:56
    end

end

%stampa solo all'interno
%y_hat del massimo NOX
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), t,'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');

 
%y_hat del massimo NOX
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,1}(:,:,idx_max),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');

%std_hat del massimo NOX
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,2}(:,:,idx_max),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');

%y_hat del PM10 MASSIMO
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,5}(:,:,idx_max2),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(result_krig_bivariato{1,7}(:,1), result_krig_bivariato{1,7}(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');



%std_hat del PM10 massimo
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,6}(:,:,idx_max2),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(result_krig_bivariato{1,7}(:,1), result_krig_bivariato{1,7}(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');



%stampa della z_k relativa al NOx nel giorno più inquinato
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,4}(:,:,idx_max,1),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');


%z_k del pm10 nel giono più inquinato
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), result_krig_bivariato{1,4}(:,:,idx_max2,2),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{2,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','r');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');

% verifica che le y siano stazionarie (Nsomma)
plot(reshape(result_krig_bivariato{1,1}(1,1, :), 1, 365))
adftest(reshape(result_krig_bivariato{1,1}(1,1, :), 1, 365))

%% plot of the levation in the grid
figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), reshape(a(1:end,3),56,56),'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
colormap('pink');
xlim([-3.9 -3.5]);
ylim([40.3 40.7]);
xlabel('Longitude');
ylabel('Latitude');

%% plot della media della distribuzione dell'NO_x durante l'anno
y_hat_medio = mean(result_krig_bivariato{1,1}, 3);

figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56), y_hat_medio,'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');

%% verifica della stazionarietà della z e media 0
z_media_grid = mean(result_krig_bivariato{1,4}(:,:,:,1),3);
z_media_punti = mean(z_media_grid, 1);
mean(z_media_punti) %non ha media 0

%stazionarietà
z_temp = mean(result_krig_bivariato{1,4}(:,:,:,1),1:2);
adftest(reshape(z_temp, 365,1)) %almeno è stazionario

%% tabella dei valori
load("result_data_nox_univariate.mat")
load("result_data_pm10_univariate.mat")
load("result_data_pm10_univariate_selected.mat")
load("result_data_nox_univariate_selected.mat")

load("result_data_biavariate_NOX.mat")
load("result_data_biavariate_NOX_selected.mat")
load("result_data_bivariate_PM10.mat")
load("result_data_bivariate_PM10_selected.mat")

load("result_data_multivariato_meteo_NOX.mat")
load("result_data_multivariato_meteo_PM10.mat")
load("result_data_multivariato_meteo_PM10_selected.mat")

% estrazione della log-likelihood
l1 = mean(result_data_nox_univariate{8})
l2 = mean(result_data_pm10_univariate{8})
l3 = mean(result_data_pm10_univariate_selected{8})
l4 = mean(result_data_nox_univariate_selected{8})


l5 = mean(result_data_biavariate_NOX{8})
l6 = mean(result_data_biavariate_NOX_selected{8})
l7 = mean(result_data_biavariate_PM10{8})
l8 = mean(result_data_biavariate_PM10_selected{8})

l9 = mean(result_data_multivariato_meteo_NOX{8})
l10 = mean(result_data_multivariato_meteo_PM10{8})
l11 = mean(result_data_multivariato_meteo_PM10_selected{8})



result_data_nox_univaria


%% Prova di ritagliamento

madrid = shaperead('../madrid-districtsgeojson.shp');
a = load('./../kriging_regulaGrid_elevations.csv');
X = reshape(a(1:end,2),56,56);
Y = reshape(a(1:end,1),56,56);

inpolygon(X, Y, )

%% Grafici nox

madrid = shaperead('..\madrid-districtsgeojson.shp');
bo = madrid(8).BoundingBox

% concatenameto

function frequenza = contaValore(vettore, valore)
    % Inizializza la variabile di conteggio
    conteggio = 0;

    % Itera attraverso il vettore
    for i = 1:length(vettore)
        % Controlla se l'elemento corrente è uguale al valore desiderato
        if find(vettore(i) == valore)
            % Incrementa il conteggio
            conteggio = conteggio + 1;
        end
    end

    % Restituisci il numero di occorrenze
    frequenza = conteggio;
end


X = [];
Y = [];

for i = 1:size(madrid,1)
    X = [X madrid(i).X];
    Y = [Y madrid(i).Y];
    
end

coord = [X; Y]';

A = unique(coord, 'rows', 'stable');
B = setdiff(coord, A, 'stable'); 
C = setdiff(A, B, 'rows', 'stable');


%selezione dei non duplicati
X_cont = [];
Y_cont = [];

groupcounts()

%coord_cont = [];

for i=1:size(X,2)
    for j= 1:size(X,2)
        if coord(i,:) == coord(j,:) & (i ~= j)
            coord_cont = [coord_cont; coord(i,:)];
        end
    end
end

% in coord cont ho quelle che compaiono più di una volta
count = 0;
for i=1:size(X,2)
    for j= 1:size(coord_cont,2)
        if coord(i,1) == coord_cont(j,1) & coord(i,2) == coord_cont(j,2) 
            coord(i,1) = NaN;
            coord(i,2) = NaN;
            count = count +1;
        end
    end
end


% Rimuovere le righe con NaN
matrice_senza_nan = matrice(~indici_righe_nan, :);



%% codice (forse) funzionante
madrid = shaperead('../madrid-districtsgeojson.shp');
a = load('./../kriging_regulaGrid_elevations.csv');
Xg = reshape(a(1:end,2),56,56);
Yg = reshape(a(1:end,1),56,56);



X = [];
Y = [];

for i = 1:size(madrid,1)
    X = [X madrid(i).X];
    Y = [Y madrid(i).Y];
end

coord = [X; Y]';
isgood = isnan(coord);
coord = coord(not(isgood(:,1)),:);

A = unique(coord, 'rows', 'stable');

[ga,gr] = groupcounts(coord);

check = ga == 1 ;

X1 = gr{1,1}(check,:);
Y1 = gr{1,2}(check,:);

PUNTI = inpolygon(Xg, Yg, X1, Y1);
print = y_hat_univariato_NOX(:,:,1);

for k =1:56
    for j =1:56
       if PUNTI(k,j) == 0
            print(k,j) = NaN;
       end
    end
end

%%
figure
colorbar
caxis([0 100])
hold on

h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56),print,'DisplayType','texturemap');
set(h,'FaceColor','flat');
geoshow(madrid,'FaceColor','none');
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','y');
%%

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

%% Autoregressivo con LOGOCV (prendo la media delle stazioni di test)
% univariato ovviamente

% LOGOCV
% y è N*n
% X_reg è N*n*d
idx_totali = 1:24*365;

%prima generiamo tutte le matrici
y_reg = reshape(dati_NOX,1, 365*24)';
X_lat = reshape(ground.X_beta{1,1}(:,1,:),24, 365);
X_alt = reshape(ground.X_beta{1,1}(:,2,:),24, 365);
X_lat_def = reshape(X_lat,1,365*24)';
X_alt_def = reshape(X_alt,1,365*24)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    temp = [NaN(dim,1); y_reg(1:end-dim)];
    X_autoreg = [X_autoreg temp];
end

X_reg_NOx = [X_lat_def X_alt_def  X_autoreg];
rmse_cv_lm = [];
r2_cv_lm = [];

for t = 1:365:24*365-365
    idx_test = t:t+365;
    idx_train = setdiff(idx_totali, idx_test, 'stable');

    mdl_reg_nox = fitlm(X_reg_NOx(idx_train,:), y_reg(idx_train), 'Intercept', false);
    y_hat_nox_reg = predict(mdl_reg_nox, X_reg_NOx(idx_test,:));

    rmse_reg = nanstd(y_reg(idx_test)-y_hat_nox_reg);
    r2_reg = 1- (nanvar(y_reg(idx_test)-y_hat_nox_reg)/nanvar(y_reg(idx_train)));
    rmse_cv_lm = [rmse_cv_lm rmse_reg];
    

    r2_cv_lm = [r2_cv_lm r2_reg]
end

rmse_cv_nox_reg = mean(rmse_cv_lm)
rmse_r2_nox_reg = mean(r2_cv_lm)

%per pm10

% LOGOCV
% y è N*n
% X_reg è N*n*d
idx_totali = 1:12*365;

%prima generiamo tutte le matrici
y_reg = reshape(dati_PM25,1, 365*12)';
X_alt = reshape(ground.X_beta{1,2}(:,1,:),12, 365);
X_alt_def = reshape(X_alt,1,365*12)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    temp = [NaN(dim,1); y_reg(1:end-dim)];
    X_autoreg = [X_autoreg temp];
end

X_reg_pm25 = [X_alt_def  X_autoreg];
rmse_cv_lm = [];
r2_cv_lm = [];

for t = 1:365:12*365-365
    idx_test = t:t+365;
    idx_train = setdiff(idx_totali, idx_test, 'stable');

    mdl_reg_pm25 = fitlm(X_reg_pm25(idx_train,:), y_reg(idx_train), 'Intercept', false);
    y_hat_pm25_reg = predict(mdl_reg_pm25, X_reg_pm25(idx_test,:));

    rmse_reg = nanstd(y_reg(idx_test)-y_hat_pm25_reg);
    
    r2_reg = 1- (nanvar(y_reg(idx_test)-y_hat_pm25_reg)/nanvar(y_reg(idx_train)));

    r2_cv_lm = [r2_cv_lm r2_reg];
    rmse_cv_lm = [rmse_cv_lm rmse_reg];
    
end

rmse_cv_pm10_reg = mean(rmse_cv_lm);
r2_cv_pm10_reg = mean(r2_cv_lm);

%% rmse reg usando tutti i dati

%NOX
%valori precedenti
y_reg = reshape(dati_NOX,1, 365*24)';
X_lat = reshape(ground.X_beta{1,1}(:,1,:),24, 365);
X_alt = reshape(ground.X_beta{1,1}(:,2,:),24, 365);
X_lat_def = reshape(X_lat,1,365*24)';
X_alt_def = reshape(X_alt,1,365*24)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    X_autoreg = [X_autoreg y_reg(dim+1:end-(ord-dim)-1)];
end

X_reg_NOx = [X_lat_def(1:end-ord-1) X_alt_def(1:end-ord-1)  X_autoreg];

mean_nox = nanmean(X_reg_NOx,1);
std_nox = nanstd(X_reg_NOx,1);

mean_y_nox = nanmean(y_reg(1:end-ord-1));
std_y_nox = nanstd(y_reg(1:end-ord-1));

mdl_reg_nox = fitlm((X_reg_NOx-mean_nox)./std_nox, (y_reg(1:end-ord-1)-mean_y_nox)./std_y_nox, 'Intercept', false);
mdl_reg_nox.Rsquared.Adjusted



% non std
%NOX
%valori precedenti
y_reg = reshape(dati_NOX,1, 365*24)';
X_lat = reshape(ground.X_beta{1,1}(:,1,:),24, 365);
X_alt = reshape(ground.X_beta{1,1}(:,2,:),24, 365);
X_lat_def = reshape(X_lat,1,365*24)';
X_alt_def = reshape(X_alt,1,365*24)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    X_autoreg = [X_autoreg y_reg(dim+1:end-(ord-dim)-1)];
end

X_reg_NOx = [X_lat_def(1:end-ord-1) X_alt_def(1:end-ord-1)  X_autoreg];

mdl_reg_nox = fitlm(X_reg_NOx, y_reg(1:end-ord-1), 'Intercept', false);
mdl_reg_nox.Rsquared.Adjusted
mdl_reg_nox.RMSE

%PM10
%valori precedenti
y_reg = reshape(dati_PM25,1, 365*12)';
X_alt = reshape(ground.X_beta{1,2}(:,1,:),12, 365);
X_alt_def = reshape(X_alt,1,365*12)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    X_autoreg = [X_autoreg y_reg(dim+1:end-(ord-dim)-1)];
end

X_reg_PM10 = [X_alt_def(1:end-ord-1)  X_autoreg];
mdl_reg_pm10 = fitlm(normalize(X_reg_PM10), normalize(y_reg(1:end-ord-1)), 'Intercept', false);
mdl_reg_pm10.Rsquared.Adjusted
pm10_rmse = mdl_reg_pm10.RMSE;

%non std
y_reg = reshape(dati_PM25,1, 365*12)';
X_alt = reshape(ground.X_beta{1,2}(:,1,:),12, 365);
X_alt_def = reshape(X_alt,1,365*12)';
X_autoreg = []; 
ord = 5;
for dim=1:ord
    X_autoreg = [X_autoreg y_reg(dim+1:end-(ord-dim)-1)];
end

X_reg_PM10 = [X_alt_def(1:end-ord-1)  X_autoreg];
mdl_reg_pm10 = fitlm(X_reg_PM10, y_reg(1:end-ord-1), 'Intercept', false);
mdl_reg_pm10.Rsquared.Adjusted
pm10_rmse = mdl_reg_pm10.RMSE;
