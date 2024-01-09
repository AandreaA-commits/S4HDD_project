clc
clearvars

%% setup dei dati
addpath('../../Src'); %D-STEM
a = load("..\..\..\..\dati_pollutants_e_scraper\dati_pollutants_2019.mat");
dati_pollutants_2019 = a.dati_meteo_2019;
load ..\..\..\..\dati_traffico_2019.mat
load ..\..\..\..\dati_meteo_2019.mat

rng(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NOX = dati_pollutants_2019(7, 2:end);
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
            temp = [temp; PM25{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PM25{1,mese}(:,1)) == table2array(temp(:,1)));
        PM25{1,mese} = PM25{1,mese}(rowremove, :);
    end
end
% NOX
min_col = size(NOX{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(NOX{1,iii},1)
        min_col = size(NOX{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(NOX{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(NOX{1,mese}(:,1));
    temp = [];
    for row = 1: size(NOX{1,mese}, 1)
        if ~ismember(NOX{1,mese}(row,1), stazioni_comuni)
           
            temp = [temp; NOX{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(NOX{1,mese}(:,1)) == table2array(temp(:,1)));
        NOX{1,mese} = NOX{1,mese}(rowremove, :);       
        
    end
end


dati_NOX = [];
dati_PM25 = [];

for i = 1:size(NOX, 2)
    % NOX
    tabella_corrente = NOX{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NOX = [dati_NOX colonne_numeriche];
    % pm2.5
    tabella_corrente = PM25{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PM25 = [dati_PM25 colonne_numeriche];    
end


% creazione vettore covariata is_weekend
% 1 gennaio 2019 era martedì
is_weekend = zeros(1,365*24);
flag_weekend = 48;
flag_week = 24*5;
for i = 4*24+1:size(is_weekend,2) % itero tutte le ore dell'anno 
    if flag_weekend ~= 0
        is_weekend(1,i) = 1;
        flag_weekend = flag_weekend -1; 
        flag_week = 24*5;
    else
        if flag_week == 0
            flag_weekend = 48;
        else
            flag_week = flag_week-1; 
        end
    end    
end

%% script per divisione dati di training e dati di testing (per stazione)
%numero delle stazioni totali
ns1 = size(dati_NOX, 1); 
ns2 = size(dati_PM25, 1); 

% Specifica la percentuale desiderata di righe da estrarre
percentuale_righe = 0.75;

% Calcola il numero desiderato di righe
numero_righe1 = round(percentuale_righe * ns1);
numero_righe2 = round(percentuale_righe * ns2);

% indici di train e test
indici_totali1 = 1:ns1;
indici_righe_train1 = randperm(ns1, numero_righe1);
indici_righe_test1 = setdiff(indici_totali1, indici_righe_train1);

indici_totali2 = 1:ns2;
indici_righe_train2 = randperm(ns2, numero_righe2);
indici_righe_test2 = setdiff(indici_totali2, indici_righe_train2);

% Estrazione dati train e test
dati_train_NOX = dati_NOX(indici_righe_train1, :);
dati_test_NOX = dati_NOX(indici_righe_test1, :);

dati_train_PM25 = dati_PM25(indici_righe_train2, :);
dati_test_PM25 = dati_PM25(indici_righe_test2, :);

%load no2 obs
ground.Y{1} = dati_train_NOX;
ground.Y_name{1} = 'nox';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load pm2.5 obs
ground.Y{2} = dati_train_PM25;
ground.Y_name{2} = 'pm25';
n2 = size(ground.Y{2}, 1);

NOx_lat = NOX{1,1}{:,3};
NOx_long = NOX{1,1}{:,4};
NOx_alt = NOX{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n1, 1, T);
X_krig = zeros(size(dati_test_NOX, 1), 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n1,1);
        X_krig(:,1,i) = zeros(size(dati_test_NOX, 1),1);
    else
        X(:,1,i) = ones(n1,1);
        X_krig(:,1,i) = ones(size(dati_test_NOX, 1),1);
    end 
    X(:,2,i) = NOx_lat(indici_righe_train1);
    X(:,3,i) = NOx_long(indici_righe_train1);  
    X(:,4,i) = NOx_alt(indici_righe_train1); 
    X_krig(:,2,i) = NOx_lat(indici_righe_test1);
    X_krig(:,3,i) = NOx_long(indici_righe_test1);
    X_krig(:,4,i) = ones(size(dati_test_NOX, 1),1);
    X_krig(:,5,i) = NOx_alt(indici_righe_test1);
end
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'weekend', 'lat', 'long','alt'};
ground.X_beta_name_krig{1} = {'weekend', 'lat', 'long', 'constant','alt'};
ground.X_beta_krig{1} = X_krig;


PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};
PM25_alt = PM25{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n2, 1, T);
X_krig = zeros(size(dati_test_PM25, 1), 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n2,1);
        X_krig(:,1,i) = zeros(size(dati_test_PM25, 1),1);      
    else
        X(:,1,i) = ones(n2,1);
        X_krig(:,1,i) = ones(size(dati_test_PM25, 1),1);        
    end
    X(:,2,i) = PM25_lat(indici_righe_train2);
    X(:,3,i) = PM25_long(indici_righe_train2);
    X(:,4,i) = PM25_alt(indici_righe_train2);
    X_krig(:,2,i) = PM25_lat(indici_righe_test2);
    X_krig(:,3,i) = PM25_long(indici_righe_test2); 
    X_krig(:,4,i) = ones(size(dati_test_PM25, 1),1); 
    X_krig(:,5,i) = PM25_alt(indici_righe_test2); 
end
ground.X_beta{2} = X;
ground.X_beta_name{2} = {'weekend', 'lat', 'long','alt'};
ground.X_beta_name_krig{2} = {'weekend', 'lat', 'long', 'constant','alt'};
ground.X_beta_krig{2} = X_krig;


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

ground.X_z{2} = ones(n2, 1);
ground.X_z_name{2} = {'constant'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);

%laod of the station coordinates
NOx_lat = NOX{1,1}{:,3};
NOx_long = NOX{1,1}{:,4};

PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [NOx_lat(indici_righe_train1), NOx_long(indici_righe_train1)];
ground.coordinates{2} = [PM25_lat(indici_righe_train2), PM25_long(indici_righe_train2)];


obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');


obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2019 00:00','31-12-2019 00:00',T);

%stem_data object creation
shape = [];
%obj_stem_validation=[];
S_val1=1:5:n1;
S_val2=1:5:n2;

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
obj_stem_par.v_z = eye(2)*0.1;
obj_stem_par.sigma_eta = diag([0.02 0.02]);
obj_stem_par.G = diag(0.9*ones(2,1));
obj_stem_par.sigma_eps = diag([0.01 0.3]); 

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

a = 0;

%% Kriging on validation stations

% KRINGING NOX
krig_coordinates = [NOx_lat(indici_righe_test1), NOx_long(indici_righe_test1)];

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'sparse','point');

obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid, ground.X_beta_krig{1,1}, ground.X_beta_name_krig{1,1});
obj_stem_krig = stem_krig(obj_stem_model, obj_stem_krig_data);

obj_stem_krig_options = stem_krig_options();
obj_stem_krig_options.block_size = 1000;

obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);

%calcolo dell'RMSE e R2
y_hat_nox = obj_stem_krig_result{1}.y_hat;

% prendiamo le y originali
rmse_nox = [];
r2_nox = [];

rmse_nox = nanstd(dati_test_NOX - y_hat_nox,1,2)

r2_nox = 1 - nanvar(dati_test_NOX - y_hat_nox,1,2)./nanvar(dati_test_NOX,1,2);
rmse_nox_mean = mean(rmse_nox)
r2_nox_mean = mean(r2_nox)


%%
% KRINGING PM25
krig_coordinates = [PM25_lat(indici_righe_test2), PM25_long(indici_righe_test2)];

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'sparse','point',[], 'square', 0.001, 0.001);

obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid, ground.X_beta_krig{1,2}, ground.X_beta_name_krig{1,2});
obj_stem_krig = stem_krig(obj_stem_model, obj_stem_krig_data);

obj_stem_krig_options = stem_krig_options();
obj_stem_krig_options.block_size = 1000;

obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);

%calcolo dell'RMSE e R2
y_hat_pm25 = obj_stem_krig_result{1}.y_hat;

% prendiamo le y originali
rmse_pm25 = nanstd(dati_test_PM25 - y_hat_pm25,1,2)
r2_pm25 = 1 - nanvar(dati_test_PM25 - y_hat_pm25,1,2)./nanvar(dati_test_PM25,1,2);

rmse_pm25_mean = mean(rmse_pm25);
r2_pm25_mean = mean(r2_pm25)
