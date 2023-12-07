%clc
%clearvars

%{
repositoryPath = fullfile(pwd, 'GitHub\S4HDD_project\codice');
addpath(genpath(repositoryPath));

nomeFile = 'dati_meteo_2019.mat';

percorsoCompleto = fullfile(repositoryPath, nomeFile);

load(percorsoCompleto);
%}

%% setup dei dati
addpath('C:\Users\arici\Documents\GitHub\S4HDD_project\codice\Codice_dstem\D-STEM\Src'); %D-STEM

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

NO2 = dati_pollutants_2019(4, 2:end);
dati_NO2 = [];

for i = 1:size(NO2, 2)
    tabella_corrente = NO2{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NO2 = [dati_NO2 colonne_numeriche];
end

traffico = str2double(totale(:, 4:end));

%load no2 obs
ground.Y{1} = dati_NO2;
ground.Y_name{1} = 'no2';
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%load traffic obs
ground.Y{2} = traffico;
ground.Y_name{2} = 'traffic';
n2 = size(ground.Y{2}, 1);


%load covariates for the NO2 monitoring stations (ci mettiamo l'altezza delle stazioni)
% devo prendere solo la colonna della elevation
dati_elevation = tabella_corrente{:, 2};

ground.X_beta{1} = dati_elevation;
ground.X_beta_name{1} = {'elevation'};
ground.X_beta{1} = ones();

%X_z
ground.X_z{1} = ones(n1, 1);
ground.X_z_name{1} = {'constant'};

ground.X_z{2} = ones(n2, 1);
ground.X_z_name{2} = {'constant'};

obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                [], [], ... 
                                ground.X_z, ground.X_z_name);


NO2_lat = tabella_corrente{:,3};
NO2_long = tabella_corrente{:,4};

traffic_lat = str2double(totale(:,2));
traffic_long = str2double(totale(:,3));
                            
obj_stem_gridlist_p = stem_gridlist();

ground.coordinates{1} = [NO2_lat, NO2_long];
ground.coordinates{2} = [traffic_lat, traffic_long];

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
obj_stem_validation=stem_validation({'no2','traffic'},{S_val1,S_val2},0,{'point','point'});

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
obj_stem_par.v_z = [4 2.5; 2.5 2];
obj_stem_par.sigma_eta = diag([0.2 0.2]);
obj_stem_par.G = diag([0.9 -0.4]);
obj_stem_par.sigma_eps = diag([0.02 0.6]);

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
