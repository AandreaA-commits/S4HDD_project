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

%% parameter setting

dati_PM25 = [];

for i = 1:size(PM25, 2)
    % pm2.5
    tabella_corrente = PM25{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PM25 = [dati_PM25 colonne_numeriche];
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

%load pm2.5 obs
ground.Y{1} = dati_PM25;
ground.Y_name{1} = 'pm25';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

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

%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};


obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                ground.X_beta, ground.X_beta_name, ... 
                                ground.X_z, ground.X_z_name);



PM25_lat = PM25{1,1}{:,3};
PM25_long = PM25{1,1}{:,4};

obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [PM25_lat, PM25_long];

obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');

obj_stem_gridlist_p.add(obj_stem_grid1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2019 00:00','31-12-2019 00:00', 8760);

S_val1=1:5:n1;

shape = [];

obj_stem_validation = stem_validation({'pm25'}, ... 
    {S_val1},0,{'point'});

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
obj_stem_par.theta_z = 0.1;
obj_stem_par.v_z = 1;
obj_stem_par.sigma_eta = 1;
obj_stem_par.G = 0.9;
obj_stem_par.sigma_eps = 1; 

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

