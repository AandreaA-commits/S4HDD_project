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

VEL_VENTO = dati_meteo_2019(2, 2:end)
UMIDITA = dati_meteo_2019(5, 2:end)
TEMPERATURA = dati_meteo_2019(4, 2:end)
PRESSIONE = dati_meteo_2019(6, 2:end)


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

% TEMPERATURA
min_col = size(TEMPERATURA{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(TEMPERATURA{1,iii},1)
        min_col = size(TEMPERATURA{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(TEMPERATURA{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(TEMPERATURA{1,mese}(:,1));
    temp = [];
    for row = 1: size(TEMPERATURA{1,mese}, 1)
        if ~ismember(TEMPERATURA{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; TEMPERATURA{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(TEMPERATURA{1,mese}(:,1)) == table2array(temp(1,1)));
        TEMPERATURA{1,mese} = TEMPERATURA{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end
min_col = size(TEMPERATURA{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(TEMPERATURA{1,iii},1)
        min_col = size(TEMPERATURA{1,iii},1);
        col = iii;
    end
end
stazioni_comuni = unique(TEMPERATURA{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(TEMPERATURA{1,mese}(:,1));
    temp = [];
    for row = 1: size(TEMPERATURA{1,mese}, 1)
        if ~ismember(TEMPERATURA{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; TEMPERATURA{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(TEMPERATURA{1,mese}(:,1)) == table2array(temp(1,1)));
        TEMPERATURA{1,mese} = TEMPERATURA{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end

% PRESSIONE
min_col = size(PRESSIONE{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PRESSIONE{1,iii},1)
        min_col = size(PRESSIONE{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(PRESSIONE{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PRESSIONE{1,mese}(:,1));
    temp = [];
    for row = 1: size(PRESSIONE{1,mese}, 1)
        if ~ismember(PRESSIONE{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PRESSIONE{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PRESSIONE{1,mese}(:,1)) == table2array(temp(1,1)));
        PRESSIONE{1,mese} = PRESSIONE{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end
min_col = size(PRESSIONE{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(PRESSIONE{1,iii},1)
        min_col = size(PRESSIONE{1,iii},1);
        col = iii;
    end
end
stazioni_comuni = unique(PRESSIONE{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(PRESSIONE{1,mese}(:,1));
    temp = [];
    for row = 1: size(PRESSIONE{1,mese}, 1)
        if ~ismember(PRESSIONE{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; PRESSIONE{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(PRESSIONE{1,mese}(:,1)) == table2array(temp(1,1)));
        PRESSIONE{1,mese} = PRESSIONE{1,mese}(rowremove, :); 
        rowremove
        disp("COCK")       
    end
end

% UMIDITA
min_col = size(UMIDITA{1,1},1);
col = 1;
for iii = 2:12
    if min_col > size(UMIDITA{1,iii},1)
        min_col = size(UMIDITA{1,iii},1);
        col = iii;
    end
end

stazioni_comuni = unique(UMIDITA{1,col}(:,1));
for mese = 1:12
    stazioni_mese = unique(UMIDITA{1,mese}(:,1));
    temp = [];
    for row = 1: size(UMIDITA{1,mese}, 1)
        if ~ismember(UMIDITA{1,mese}(row,1), stazioni_comuni)
            disp("COCCO")
            temp = [temp; UMIDITA{1,mese}(row, :)];
        end
    end
    if ~size(temp,1) == 0         
        rowremove = ~(table2array(UMIDITA{1,mese}(:,1)) == table2array(temp(:,1)));
        UMIDITA{1,mese} = UMIDITA{1,mese}(rowremove, :);       
        disp("COCK")
    end
end

dati_NOX = [];
dati_VELVENTO = [];
dati_UMIDITA = [];
dati_TEMPERATURA = [];
dati_PRESSIONE = [];


for i = 1:size(NOX, 2)
    % NOX
    tabella_corrente = NOX{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NOX = [dati_NOX colonne_numeriche];
    % VEL VENTO
    tabella_corrente = VEL_VENTO{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_VELVENTO = [dati_VELVENTO colonne_numeriche];   
    % UMIDITA
    tabella_corrente = UMIDITA{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_UMIDITA = [dati_UMIDITA colonne_numeriche];
    % TEMPERATURA
    tabella_corrente = TEMPERATURA{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_TEMPERATURA = [dati_TEMPERATURA colonne_numeriche];
    % PRESSIONE
    tabella_corrente = PRESSIONE{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PRESSIONE = [dati_PRESSIONE colonne_numeriche];
end

dati_NOX = dati_NOX(:,14:24:end)  
dati_VELVENTO = dati_VELVENTO(:,14:24:end)  
dati_UMIDITA= dati_UMIDITA(:,14:24:end)  
dati_TEMPERATURA = dati_TEMPERATURA(:,14:24:end)  
dati_PRESSIONE = dati_PRESSIONE(:,14:24:end)  

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

ns3= size(dati_TEMPERATURA, 1); 
ns4 = size(dati_UMIDITA, 1); 
ns5 = size(dati_VELVENTO, 1); 
ns6 = size(dati_PRESSIONE, 1); 


% Specifica la percentuale desiderata di righe da estrarre
percentuale_righe = 0.75;

% Calcola il numero desiderato di righe
numero_righe1 = round(percentuale_righe * ns1);

numero_righe3 = round(percentuale_righe * ns3);
numero_righe4 = round(percentuale_righe * ns4);
numero_righe5 = round(percentuale_righe * ns5);
numero_righe6 = round(percentuale_righe * ns6);

% indici di train e test
indici_totali1 = 1:ns1;
indici_righe_train1 = randperm(ns1, numero_righe1);
indici_righe_test1 = setdiff(indici_totali1, indici_righe_train1);

indici_righe_train1 = [18 22	19 14	13	6	20	3	7	23	17	16	9	10	12	21	11	15]';
indici_righe_test1 = [1 2 4 5 8 24]';


indici_totali3 = 1:ns3;
indici_righe_train3 = randperm(ns3, numero_righe3);
indici_righe_test3 = setdiff(indici_totali3, indici_righe_train3);

indici_totali4 = 1:ns4;
indici_righe_train4 = randperm(ns4, numero_righe4);
indici_righe_test4 = setdiff(indici_totali4, indici_righe_train4);

indici_totali5 = 1:ns5;
indici_righe_train5 = randperm(ns5, numero_righe5);
indici_righe_test5 = setdiff(indici_totali5, indici_righe_train5);

indici_totali6 = 1:ns6;
indici_righe_train6 = randperm(ns6, numero_righe6);
indici_righe_test6 = setdiff(indici_totali6, indici_righe_train6);

% Estrazione dati train e test
dati_train_NOX = dati_NOX(indici_righe_train1, :);
dati_test_NOX = dati_NOX(indici_righe_test1, :);

dati_train_TEMPERATURA = dati_TEMPERATURA(indici_righe_train3, :);
dati_test_TEMPERATURA = dati_TEMPERATURA(indici_righe_test3, :);

dati_train_UMIDITA = dati_UMIDITA(indici_righe_train4, :);
dati_test_UMIDITA = dati_UMIDITA(indici_righe_test4, :);

dati_train_VELVENTO = dati_VELVENTO(indici_righe_train5, :);
dati_test_VELVENTO = dati_VELVENTO(indici_righe_test5, :);

dati_train_PRESSIONE = dati_PRESSIONE(indici_righe_train6, :);
dati_test_PRESSIONE = dati_PRESSIONE(indici_righe_test6, :);

%load no2 obs
ground.Y{1} = dati_train_NOX;
ground.Y_name{1} = 'nox';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);


%load temp obs
ground.Y{2} = dati_TEMPERATURA;
ground.Y_name{2} = 'temp';
n3 = size(ground.Y{2}, 1);

%load umidity obs
ground.Y{3} = dati_UMIDITA;
ground.Y_name{3} = 'umidity';
n4 = size(ground.Y{3}, 1);

%load wind velocity obs
ground.Y{4} = dati_VELVENTO;
ground.Y_name{4} = 'wind_vel';
n5 = size(ground.Y{4}, 1);

%load precipitation obs
ground.Y{5} = dati_PRESSIONE;
ground.Y_name{5} = 'press';
n6 = size(ground.Y{5}, 1);

NOx_lat = NOX{1,1}{:,3};
NOx_long = NOX{1,1}{:,4};
NOx_alt = NOX{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n1, 1, T);
X_krig = zeros(size(dati_test_NOX, 1), 1, T);
for i=1:T
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

TEMPERATURA_lat = TEMPERATURA{1,1}{:,3};
TEMPERATURA_long = TEMPERATURA{1,1}{:,4};
TEMPERATURA_alt = TEMPERATURA{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n3, 1, T);
for i=1:T  
    X(:,1,i) = TEMPERATURA_lat;
    X(:,2,i) = TEMPERATURA_long;
    X(:,3,i) = TEMPERATURA_alt;
end
ground.X_beta{2} = X;
ground.X_beta_name{2} = {'lat', 'long','alt'};
ground.X_beta_krig{2} = X;

UMIDITA_lat = UMIDITA{1,1}{:,3};
UMIDITA_long = UMIDITA{1,1}{:,4};
UMIDITA_alt = UMIDITA{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n4, 1, T);
for i=1:T  
    X(:,1,i) = UMIDITA_lat;
    X(:,2,i) = UMIDITA_long;
    X(:,3,i) = UMIDITA_alt;
end
ground.X_beta{3} = X;
ground.X_beta_name{3} = {'lat', 'long','alt'};
ground.X_beta_krig{3} = X;


VEL_VENTO_lat = VEL_VENTO{1,1}{:,3};
VEL_VENTO_long = VEL_VENTO{1,1}{:,4};
VEL_VENTO_alt = VEL_VENTO{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n5, 1, T);
for i=1:T
    X(:,1,i) = VEL_VENTO_lat;
    X(:,2,i) = VEL_VENTO_long;
    X(:,3,i) = VEL_VENTO_alt;
end
ground.X_beta{4} = X;
ground.X_beta_name{4} = {'lat', 'long','alt'};
ground.X_beta_krig{4} = X;


PRESSIONE_lat = PRESSIONE{1,1}{:,3};
PRESSIONE_long = PRESSIONE{1,1}{:,4};
PRESSIONE_alt = PRESSIONE{1,1}{:,2};
%matrice [stazioni x numero_covariate x giorni]
X = zeros(n6, 1, T);
for i=1:T
    X(:,1,i) = PRESSIONE_lat;
    X(:,2,i) = PRESSIONE_long;
    X(:,3,i) = PRESSIONE_alt;
end
ground.X_beta{5} = X;
ground.X_beta_name{5} = {'lat', 'long','alt'};
ground.X_beta_krig{5} = X;


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

ground.X_z{2} = ones(n3, 1);
ground.X_z_name{2} = {'constant'};

ground.X_z{3} = ones(n4, 1);
ground.X_z_name{3} = {'constant'};

ground.X_z{4} = ones(n5, 1);
ground.X_z_name{4} = {'constant'};

ground.X_z{5} = ones(n6, 1);
ground.X_z_name{5} = {'constant'};


obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);

%laod of the station coordinates
NOx_lat = NOX{1,1}{:,3};
NOx_long = NOX{1,1}{:,4};

UMIDITA_lat = UMIDITA{1,1}{:,3};
UMIDITA_long = UMIDITA{1,1}{:,4};

PRESSIONE_lat = PRESSIONE{1,1}{:,3};
PRESSIONE_long = PRESSIONE{1,1}{:,4};

VEL_VENTO_lat = VEL_VENTO{1,1}{:,3};
VEL_VENTO_long = VEL_VENTO{1,1}{:,4};
                            
TEMPERATURA_lat = TEMPERATURA{1,1}{:,3};
TEMPERATURA_long = TEMPERATURA{1,1}{:,4};

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [NOx_lat(indici_righe_train1), NOx_long(indici_righe_train1)];
ground.coordinates{2} = [TEMPERATURA_lat, TEMPERATURA_long];
ground.coordinates{3} = [UMIDITA_lat, UMIDITA_long];
ground.coordinates{4} = [VEL_VENTO_lat, VEL_VENTO_long];
ground.coordinates{5} = [PRESSIONE_lat, PRESSIONE_long];


obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
obj_stem_grid3 = stem_grid(ground.coordinates{3}, 'deg', 'sparse', 'point');
obj_stem_grid4 = stem_grid(ground.coordinates{4}, 'deg', 'sparse', 'point');
obj_stem_grid5 = stem_grid(ground.coordinates{5}, 'deg', 'sparse', 'point');


obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);
obj_stem_gridlist_p.add(obj_stem_grid3);
obj_stem_gridlist_p.add(obj_stem_grid4);
obj_stem_gridlist_p.add(obj_stem_grid5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2019 00:00','31-12-2019 00:00',T);

%stem_data object creation
shape = [];
%obj_stem_validation=[];
S_val1=1:5:n1;
S_val2=1:5:n3;
S_val3=1:5:n4;
S_val4=1:5:n5;
S_val5=1:5:n6;


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
obj_stem_par.v_z = rand(5,5)*0.1 + eye(5);
obj_stem_par.sigma_eta = diag([0.02 0.1 0.1 0.1 0.1]);
obj_stem_par.G = diag(0.9*ones(5,1));
obj_stem_par.sigma_eps = diag([0.01 0.02 0.1 0.1 0.2]); 

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

d = sqrt(diag(obj_stem_model.stem_par.v_z).*eye(5))
R = inv(d)*obj_stem_model.stem_par.v_z*inv(d)




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
y_hat_nox = obj_stem_krig_result{1,1}.y_hat;

% prendiamo le y originali
rmse_nox = [];
r2_nox = [];

rmse_nox = nanstd(dati_test_NOX - y_hat_nox,1,2)

r2_nox = 1 - nanvar(dati_test_NOX - y_hat_nox,1,2)./nanvar(dati_test_NOX,1,2);
rmse_nox_mean = mean(rmse_nox)
r2_nox_mean = mean(r2_nox)

