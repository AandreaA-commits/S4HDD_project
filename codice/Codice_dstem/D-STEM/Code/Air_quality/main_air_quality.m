clc
clearvars

%{
repositoryPath = fullfile(pwd, 'GitHub\S4HDD_project\codice');
addpath(genpath(repositoryPath));

nomeFile = 'dati_meteo_2019.mat';

percorsoCompleto = fullfile(repositoryPath, nomeFile);

load(percorsoCompleto);
%}

%% setup dei dati
%addpath('C:\Users\arici\Documents\GitHub\S4HDD_project\codice\Codice_dstem\D-STEM\Src'); %D-STEM
addpath('../../Src'); %D-STEM
a = load("..\..\..\..\dati_pollutants_e_scraper\dati_pollutants_2019.mat");
dati_pollutants_2019 = a.dati_meteo_2019;
load ..\..\..\..\dati_traffico_2019.mat
load ..\..\..\..\dati_meteo_2019.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NOX = dati_pollutants_2019(7, 2:end);
PM25 = dati_pollutants_2019(5, 2:end);
BEN = dati_pollutants_2019(10, 2:end);

VEL_VENTO = dati_meteo_2019(2, 2:end)
PRECIPITAZIONI = dati_meteo_2019(8, 2:end)
UMIDITA = dati_meteo_2019(5, 2:end)


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
dati_PM25 = [];
dati_BEN = [];
dati_VELVENTO = [];
dati_PRECIPITAZIONI = [];
dati_UMIDITA = [];


for i = 1:size(NOX, 2)
    % NOX
    tabella_corrente = NOX{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NOX = [dati_NOX colonne_numeriche];
    % pm2.5
    tabella_corrente = PM25{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PM25 = [dati_PM25 colonne_numeriche];
    % BEN
    tabella_corrente = BEN{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_BEN = [dati_BEN colonne_numeriche];
    % VEL VENTO
    tabella_corrente = VEL_VENTO{1, i};
    colonne_numeriche = tabella_corrente{:, 4:end};
    dati_VELVENTO = [dati_VELVENTO colonne_numeriche];
    % PRECIPITAZIONI
    tabella_corrente = PRECIPITAZIONI{1, i};
    colonne_numeriche = tabella_corrente{:, 4:end};
    dati_PRECIPITAZIONI = [dati_PRECIPITAZIONI colonne_numeriche];
    % UMIDITA
    tabella_corrente = UMIDITA{1, i};
    colonne_numeriche = tabella_corrente{:, 4:end};
    dati_UMIDITA = [dati_UMIDITA colonne_numeriche];
end

traffico = str2double(totale(:, 4:end));

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

%load no2 obs
ground.Y{1} = dati_NOX;
ground.Y_name{1} = 'nox';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load traffic obs
ground.Y{2} = traffico;
ground.Y_name{2} = 'traffic';
n2 = size(ground.Y{2}, 1);

%load pm2.5 obs
ground.Y{3} = dati_PM25;
ground.Y_name{3} = 'pm25';
n3 = size(ground.Y{3}, 1);

%load benzene obs
ground.Y{4} = dati_BEN;
ground.Y_name{4} = 'ben';
n4 = size(ground.Y{4}, 1);

%load umidity obs
ground.Y{5} = dati_UMIDITA;
ground.Y_name{5} = 'umidity';
n5 = size(ground.Y{5}, 1);

%load precipitation obs
ground.Y{6} = dati_PRECIPITAZIONI;
ground.Y_name{6} = 'prec';
n6 = size(ground.Y{6}, 1);

%load wind velocity obs
ground.Y{7} = dati_VELVENTO;
ground.Y_name{7} = 'wind_vel';
n7 = size(ground.Y{7}, 1);

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n1, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n1,1);
    else
        X(:,1,i) = ones(n1,1);
    end
end


%load of sunday flags
ground.X_beta{1} = X;
ground.X_beta_name{1} = {'weekend'};


%matrice [stazioni x numero_covariate x giorni]
X = zeros(n2, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n2,1);
    else
        X(:,1,i) = ones(n2,1);
    end
end


ground.X_beta{2} = X;
ground.X_beta_name{2} = {'weekend'};

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n3, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n3,1);
    else
        X(:,1,i) = ones(n3,1);
    end
end


ground.X_beta{3} = X;
ground.X_beta_name{3} = {'weekend'};

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n4, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n4,1);
    else
        X(:,1,i) = ones(n4,1);
    end
end


ground.X_beta{4} = X;
ground.X_beta_name{4} = {'weekend'};

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n5, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n5,1);
    else
        X(:,1,i) = ones(n5,1);
    end
end

ground.X_beta{5} = X;
ground.X_beta_name{5} = {'weekend'};

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n6, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n6,1);
    else
        X(:,1,i) = ones(n6,1);
    end
end

ground.X_beta{6} = X;
ground.X_beta_name{6} = {'weekend'};

%matrice [stazioni x numero_covariate x giorni]
X = zeros(n7, 1, T);
for i=1:size(is_weekend,2)
    if is_weekend(i) == 0
        %creiamo una matrice n_stazioni x 1
        X(:,1,i) = zeros(n7,1);
    else
        X(:,1,i) = ones(n7,1);
    end
end

ground.X_beta{7} = X;
ground.X_beta_name{7} = {'weekend'};


%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

ground.X_z{2} = ones(n2, 1);
ground.X_z_name{2} = {'constant'};

ground.X_z{3} = ones(n3, 1);
ground.X_z_name{3} = {'constant'};

ground.X_z{4} = ones(n4, 1);
ground.X_z_name{4} = {'constant'};

ground.X_z{5} = ones(n5, 1);
ground.X_z_name{5} = {'constant'};

ground.X_z{6} = ones(n6, 1);
ground.X_z_name{6} = {'constant'};

ground.X_z{7} = ones(n7, 1);
ground.X_z_name{7} = {'constant'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);

%laod of the station coordinates
NOx_lat = NOX{1,1}{:,3};
NO2_long = NOX{1,1}{:,4};

traffic_lat = str2double(totale(:,2));
traffic_long = str2double(totale(:,3));

PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};

BEN_lat = BEN{1,1}{:,3};
BEN_long = BEN{1,1}{:,4};

UMIDITA_lat = UMIDITA{1,1}{:,2};
UMIDITA_long = UMIDITA{1,1}{:,3};

PRECIPITAZIONI_lat = PRECIPITAZIONI{1,1}{:,2};
PRECIPITAZIONI_long = PRECIPITAZIONI{1,1}{:,3};

VEL_VENTO_lat = VEL_VENTO{1,1}{:,2};
VEL_VENTO_long = VEL_VENTO{1,1}{:,3};
                            
obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [NOx_lat, NO2_long];
ground.coordinates{2} = [traffic_lat, traffic_long];
ground.coordinates{3} = [PM25_lat, PM25_long];
ground.coordinates{4} = [BEN_lat, BEN_long];
ground.coordinates{5} = [UMIDITA_lat, UMIDITA_long];
ground.coordinates{6} = [PRECIPITAZIONI_lat, PRECIPITAZIONI_long];
ground.coordinates{7} = [VEL_VENTO_lat, VEL_VENTO_long];

obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
obj_stem_grid3 = stem_grid(ground.coordinates{3}, 'deg', 'sparse', 'point');
obj_stem_grid4 = stem_grid(ground.coordinates{4}, 'deg', 'sparse', 'point');
obj_stem_grid5 = stem_grid(ground.coordinates{5}, 'deg', 'sparse', 'point');
obj_stem_grid6 = stem_grid(ground.coordinates{6}, 'deg', 'sparse', 'point');
obj_stem_grid7 = stem_grid(ground.coordinates{7}, 'deg', 'sparse', 'point');

obj_stem_gridlist_p.add(obj_stem_grid1);
obj_stem_gridlist_p.add(obj_stem_grid2);
obj_stem_gridlist_p.add(obj_stem_grid3);
obj_stem_gridlist_p.add(obj_stem_grid4);
obj_stem_gridlist_p.add(obj_stem_grid5);
obj_stem_gridlist_p.add(obj_stem_grid6);
obj_stem_gridlist_p.add(obj_stem_grid7);



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2019 00:00','31-12-2019 00:00',T);

%stem_data object creation
shape = [];
%obj_stem_validation=[];
S_val1=1:5:n1;
S_val2=1:5:n2;
S_val3=1:5:n3;
S_val4=1:5:n4;
S_val5=1:5:n5;
S_val6=1:5:n6;
S_val7=1:5:n7;

obj_stem_validation = stem_validation({'nox','traffic','pm25', 'ben', 'umidity', 'prec', 'wind_vel'}, ... 
    {S_val1,S_val2,S_val3,S_val4,S_val5,S_val6,S_val7},0,{'point','point','point','point','point','point','point'});

obj_stem_modeltype = stem_modeltype('HDGM');
obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
                          [], [], obj_stem_datestamp, obj_stem_validation, obj_stem_modeltype, shape);

%stem_par object creation
obj_stem_par_constraints=stem_par_constraints();
obj_stem_par_constraints.time_diagonal=0;
obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
%stem_model object creation
obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
clear ground

%Data transform
obj_stem_model.stem_data.log_transform;
obj_stem_model.stem_data.standardize;

%Starting values
obj_stem_par.beta = obj_stem_model.get_beta0();
obj_stem_par.theta_z = 1;
obj_stem_par.v_z = eye(7);
obj_stem_par.sigma_eta = diag([0.02 0.5 0.02 0.1 0.1 0.1 0.1]);
obj_stem_par.G = diag(0.9*ones(7,1));
obj_stem_par.sigma_eps = diag([0.02 0.5 0.02 0.02 0.1 0.1 0.1]); 

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
% validation result for no2
obj_stem_model.plot_val(0,'no2');
% validation result for pm2.5
obj_stem_model.plot_val(0,'traffic');
save;
