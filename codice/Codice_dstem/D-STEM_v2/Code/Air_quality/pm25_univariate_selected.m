%% setup dei dati
%addpath('C:\Users\arici\Documents\GitHub\S4HDD_project\codice\Codice_dstem\D-STEM\Src'); %D-STEM
clc
clear all

addpath('../../Src'); %D-STEM
a = load("..\..\..\..\dati_pollutants_e_scraper\dati_pollutants_2019.mat");
dati_pollutants_2019 = a.dati_meteo_2019;
load ..\..\..\..\dati_traffico_2019.mat
load ..\..\..\..\dati_meteo_2019.mat

rng(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

PM25 = dati_pollutants_2019(5, 2:end);

%% RIMOZIONE STAZIONI CHE NON HANNO MISURAZIONI IN TUTTI I MESI DELLA ANNO
% PM25
min_col = size(PM25{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PM25{1,iii},1)
        min_col = size(PM25{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(PM25{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PM25{1,mese}(:,1));
    temp = [];
    for row = 1: size(PM25{1,mese}, 1)
        if ~ismember(PM25{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PM25{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PM25{1,mese}(:,1)) == table2array(temp(:,1)));
        PM25{1,mese} = PM25{1,mese}(rowremove, :);       
        disp("COCK")
    end
end
% PM25
min_col = size(PM25{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PM25{1,iii},1)
        min_col = size(PM25{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(PM25{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PM25{1,mese}(:,1));
    temp = [];
    for row = 1: size(PM25{1,mese}, 1)
        if ~ismember(PM25{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PM25{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PM25{1,mese}(:,1)) == table2array(temp(:,1)));
        PM25{1,mese} = PM25{1,mese}(rowremove, :);       
        disp("COCK")
    end
end

%% parameter setting
dati_PM25 = [];

for i = 1:size(PM25, 2)
    % pm2.5
    tabella_corrente = PM25{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PM25 = [dati_PM25 colonne_numeriche];
end

dati_PM25 = dati_PM25(:,14:24:end)  


%% script per divisione dati di training e dati di testing (per stazione)
%numero delle stazioni totali
n = size(dati_PM25, 1); 

% Specifica la percentuale desiderata di righe da estrarre
percentuale_righe = 0.75;

% Calcola il numero desiderato di righe
numero_righe = round(percentuale_righe * n);

% indici di train e test
indici_totali = 1:n;
indici_righe_train = randperm(n, numero_righe);
indici_righe_test = setdiff(indici_totali, indici_righe_train);

indici_righe_test = [1];
indici_righe_train = [5 3 6 4 2]';

% Estrazione dati train e test
dati_train_PM25 = dati_PM25(indici_righe_train, :);
dati_test_PM25 = dati_PM25(indici_righe_test, :);


%load PM25 obs
ground.Y{1} = dati_train_PM25;
ground.Y_name{1} = 'pm25';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

% da ripetere anche per il test questa procedura
% matrice [stazioni x numero_covariate x giorni]
PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};
PM25_alt = PM25{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n1, 1, T);
X_krig = zeros(size(dati_test_PM25, 1), 1, T);
for i=1:T     
    X(:,1,i) = PM25_lat(indici_righe_train);
    X(:,2,i) = PM25_long(indici_righe_train);  
    X(:,3,i) = PM25_alt(indici_righe_train);  
    %X_krig(:,1,i) = PM25_lat(indici_righe_test);
    X_krig(:,1,i) = ones(size(dati_test_PM25, 1),1);  
    %X_krig(:,2,i) = PM25_long(indici_righe_test);
    %X_krig(:,3,i) = ones(size(dati_test_PM25, 1),1);  
    %X_krig(:,4,i) = PM25_alt(indici_righe_test); 
end
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'lat', 'long','alt'};
ground.X_beta_name_krig{1} = {'constant'};
ground.X_beta_krig{1} = X_krig;


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};


obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                [], [], ... 
                                ground.X_z, ground.X_z_name);

PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [PM25_lat(indici_righe_train,:), PM25_long(indici_righe_train, :)];

obj_stem_grid = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid);



%% Prova digrafico
madrid = shaperead('madrid-districtsgeojson.shp');

figure
grid on
geoscatter(ground.coordinates{1}(:,1), ground.coordinates{1}(:,2), 'filled', 'b');
hold on
geoscatter(PM25_lat(indici_righe_test,:), PM25_long(indici_righe_test, :), 'filled', 'r');
for i=1:length(madrid)
    geoplot(madrid(i).Y, madrid(i).X, "k-");
end


%utilizzare geoplot per fare delle linee sulla mappa

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2019 00:00','31-12-2019 00:00', T);

shape = [];

obj_stem_validation = [];

obj_stem_modeltype = stem_modeltype('HDGM');
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, obj_stem_validation, obj_stem_modeltype, shape);

%stem_par object creation
obj_stem_par_constraints=stem_par_constraints();
obj_stem_par_constraints.time_diagonal=0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.theta_z = 0.1;
obj_stem_par.v_z = 0.2;
obj_stem_par.sigma_eta = 1;
obj_stem_par.G = 0.9;
obj_stem_par.sigma_eps = 0.1; 

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = 0.001;
max_iterations = 200;
obj_stem_EM_options = stem_EM_options();
obj_stem_EM_options.max_iterations = max_iterations;
obj_stem_EM_options.exit_tol_par = exit_toll;
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;
obj_stem_model.set_logL;


obj_stem_model.print;

%%
% kriging on validation stations
krig_coordinates = [PM25_lat(indici_righe_test, :), PM25_long(indici_righe_test, :)];

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'sparse','point');

obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid, ground.X_beta_krig{1,1}, ground.X_beta_name_krig{1,1}, []);
obj_stem_krig = stem_krig(obj_stem_model,obj_stem_krig_data);

obj_stem_krig_options = stem_krig_options();
obj_stem_krig_options.block_size = 1000;

obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);

%calcolo dell'RMSE e R2
y_hat = obj_stem_krig_result{1}.y_hat;

% prendiamo le y originali

r2 = [];
res = [];

rmse = nanstd(dati_test_PM25 - y_hat,1,2);


rmse
r2 = 1 - nanvar(dati_test_PM25 - y_hat,1,2)./nanvar(dati_test_PM25,1,2);
rmse_tot = mean(rmse);
mean(r2)


%%
% KRIGING GRID COMPLETA
a = load('./kriging_regulaGrid_elevations.csv');

X_krig = zeros(56*56, 1, T);
for i=1:T    
    X_krig(:,1,i) = a(1:end,1);
    X_krig(:,2,i) = a(1:end,2);
    X_krig(:,3,i) = ones(56*56, 1);  
    X_krig(:,4,i) = a(1:end,3);  
end
ground.X_beta_name_krig{1} = {'lat', 'long', 'constant', 'alt'};
ground.X_beta_krig{1} = X_krig;

krig_coordinates = [a(1:end,1), a(1:end,2)];

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'regular','pixel', [56, 56], 'square',0.75,0.75);
    
clear X_krig

obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid, ground.X_beta_krig{1,1}, ground.X_beta_name_krig{1,1}, []);
obj_stem_krig = stem_krig(obj_stem_model,obj_stem_krig_data);

obj_stem_krig_options = stem_krig_options();
obj_stem_krig_options.block_size = 1000;
clear ground

obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);


%%%%%%%%%% CODICE FUNZIONANTE INIZIO %%%%%%%%%%%%
%% ATTENZIOEN USARE DATI GIORNO 9 SONO MOLTO BELLI

figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56),obj_stem_krig_result{1,1}.y_hat(:,:,17),'DisplayType','texturemap')
set(h,'FaceColor','flat')
geoshow(madrid,'FaceColor','none')
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','b');

figure
colorbar
hold on
h = mapshow(reshape(a(1:end,2),56,56),reshape(a(1:end,1),56,56),obj_stem_krig_result{1,1}.diag_Var_y_hat(:,:,17),'DisplayType','texturemap')
set(h,'FaceColor','flat')
geoshow(madrid,'FaceColor','none')
geoshow(obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,1), obj_stem_krig_result{1,1}.stem_grid_sites.coordinate(:,2),...
                    'DisplayType','multipoint','Marker','*','MarkerEdgeColor','b');

%%%%%%%%%% CODICE FUNZIONANTE FINE %%%%%%%%%%%%


figure;
contourf(reshape(a(1:end,2), 56,56), reshape(a(1:end,1), 56,56), reshape(obj_stem_krig_result{1,1}.y_hat(:,16),56,56));
colorbar;
caxis([5 40])

figure
obj_stem_krig_result{1}.plot(200)

plot(obj_stem_krig_result, 0,'both')

figure
imagesc(reshape(obj_stem_krig_result{1,1}.y_hat(:,1),56,56))

figure
surf(reshape(a(1:9:end,2), 56,56)', reshape(a(1:9:end,1), 56,56)', ones(56), reshape(obj_stem_krig_result{1,1}.y_hat(:,1),56,56)')
hold on
for i=1:length(madrid)
    geoplot(madrid(i).Y, madrid(i).X, "k-");
end

figure
imagesc(reshape(obj_stem_krig_result{1, 1}.diag_Var_y_hat(:,1), 56,56))





