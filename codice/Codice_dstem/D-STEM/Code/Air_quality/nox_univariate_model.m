%% setup dei dati
%addpath('C:\Users\arici\Documents\GitHub\S4HDD_project\codice\Codice_dstem\D-STEM\Src'); %D-STEM
clc
%clear all

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

BEN = dati_pollutants_2019(10, 2:end);

VEL_VENTO = dati_meteo_2019(2, 2:end)
PRECIPITAZIONI = dati_meteo_2019(8, 2:end)
UMIDITA = dati_meteo_2019(5, 2:end)


%% RIMOZIONE STAZIONI CHE NON HANNO MISURAZIONI IN TUTTI I MESI DELLA ANNO
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
            disp("COCCO")
            temp = [temp; NOX{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(NOX{1,mese}(:,1)) == table2array(temp(:,1)));
        NOX{1,mese} = NOX{1,mese}(rowremove, :);       
        disp("COCK")
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
            disp("COCCO")
            temp = [temp; NOX{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(NOX{1,mese}(:,1)) == table2array(temp(:,1)));
        NOX{1,mese} = NOX{1,mese}(rowremove, :);       
        disp("COCK")
    end
end
% BEN
min_col = size(BEN{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(BEN{1,iii},1)
        min_col = size(BEN{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(BEN{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(BEN{1,mese}(:,1));
    temp = [];
    for row = 1: size(BEN{1,mese}, 1)
        if ~ismember(BEN{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; BEN{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(BEN{1,mese}(:,1)) == table2array(temp(:,1)));
        BEN{1,mese} = BEN{1,mese}(rowremove, :);       
        disp("COCK")
    end
end

% VEL VENTO
min_col = size(VEL_VENTO{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(VEL_VENTO{1,iii},1)
        min_col = size(VEL_VENTO{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(VEL_VENTO{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(VEL_VENTO{1,mese}(:,1));
    temp = [];
    for row = 1: size(VEL_VENTO{1,mese}, 1)
        if ~ismember(VEL_VENTO{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; VEL_VENTO{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(VEL_VENTO{1,mese}(:,1)) == table2array(temp(:,1)));
        VEL_VENTO{1,mese} = VEL_VENTO{1,mese}(rowremove, :);       
        disp("COCK")
    end
end

% PRECIPITAZIONI
min_col = size(PRECIPITAZIONI{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PRECIPITAZIONI{1,iii},1)
        min_col = size(PRECIPITAZIONI{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(PRECIPITAZIONI{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PRECIPITAZIONI{1,mese}(:,1));
    temp = [];
    for row = 1: size(PRECIPITAZIONI{1,mese}, 1)
        if ~ismember(PRECIPITAZIONI{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PRECIPITAZIONI{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PRECIPITAZIONI{1,mese}(:,1)) == table2array(temp(1,1)));
        PRECIPITAZIONI{1,mese} = PRECIPITAZIONI{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end
min_col = size(PRECIPITAZIONI{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PRECIPITAZIONI{1,iii},1)
        min_col = size(PRECIPITAZIONI{1,iii},1);
        col = iii;
    end
end
stazioni_comuni = unique(PRECIPITAZIONI{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PRECIPITAZIONI{1,mese}(:,1));
    temp = [];
    for row = 1: size(PRECIPITAZIONI{1,mese}, 1)
        if ~ismember(PRECIPITAZIONI{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PRECIPITAZIONI{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PRECIPITAZIONI{1,mese}(:,1)) == table2array(temp(1,1)));
        PRECIPITAZIONI{1,mese} = PRECIPITAZIONI{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end




%% parameter setting

dati_NOX = [];

for i = 1:size(NOX, 2)
    % pm2.5
    tabella_corrente = NOX{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NOX = [dati_NOX colonne_numeriche];
end

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
n = size(dati_NOX, 1); 

% Specifica la percentuale desiderata di righe da estrarre
percentuale_righe = 0.75;

% Calcola il numero desiderato di righe
numero_righe = round(percentuale_righe * n);

% indici di train e test
indici_totali = 1:n;
indici_righe_train = randperm(n, numero_righe);
indici_righe_test = setdiff(indici_totali, indici_righe_train);

% Estrazione dati train e test
dati_train_NOX = dati_NOX(indici_righe_train, :);
dati_test_NOX = dati_NOX(indici_righe_test, :);


%load NOX obs
ground.Y{1} = dati_train_NOX;
ground.Y_name{1} = 'nox';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

% da ripetere anche per il test questa procedura
% matrice [stazioni x numero_covariate x giorni]
NOx_lat = NOX{1,1}{:,3};
NOx_long = NOX{1,1}{:,4};
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
    X(:,2,i) = NOx_lat(indici_righe_train);
    X(:,3,i) = NOx_long(indici_righe_train);  
    X_krig(:,2,i) = NOx_lat(indici_righe_test);
    X_krig(:,3,i) = NOx_long(indici_righe_test);
    X_krig(:,4,i) = ones(size(dati_test_NOX, 1),1);    
end
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'weekend', 'lat', 'long'};
ground.X_beta_name_krig{1} = {'weekend', 'lat', 'long', 'constant'};
ground.X_beta_krig{1} = X_krig;


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};


obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);



NOX_lat = NOX{1,1}{:,3};
NOX_long = NOX{1,1}{:,4};

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [NOX_lat(indici_righe_train,:), NOX_long(indici_righe_train, :)];

obj_stem_grid = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_gridlist_p.add(obj_stem_grid);



%% Prova digrafico
madrid = shaperead('madrid-districtsgeojson.shp');

figure
grid on
geoscatter(ground.coordinates{1}(:,1), ground.coordinates{1}(:,2), 'filled', 'b');
hold on
geoscatter(NOX_lat(indici_righe_test,:), NOX_long(indici_righe_test, :), 'filled', 'r');
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
krig_coordinates = [NOX_lat(indici_righe_test, :), NOX_long(indici_righe_test, :)];

obj_stem_krig_grid = stem_grid(krig_coordinates, 'deg', 'sparse','point');

X_krig = ones(size(krig_coordinates, 1), 1, T);
for i=1:size(is_weekend,2)
    %{
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X_krig(:,1,i) = zeros(size(krig_coordinates, 1),1);
    else
        X_krig(:,1,i) = ones(size(krig_coordinates, 1),1);
    end
    %}
    %X_krig(:,1,i) = ones(size(krig_coordinates, 1),1);
end


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

rmse = nanstd(dati_test_NOX - y_hat,1,2);


rmse
r2 = 1 - nanvar(dati_test_NOX - y_hat,1,2)./nanvar(dati_test_NOX,1,2);
rmse_tot = mean(rmse);
mean(r2)


%figure;
%plot(res(1,:));
%adftest(res(1,:)) % se 1 è stazionario

% Calcolo e plot della variabile latente z(s,t)

% Creazione processo gaussiano n(s,t)
%v = mvnrnd(zeros(1,size(dist,1)), obj_stem_model.stem_EM_result.stem_par.sigma_eta,1);

% in sospeso



