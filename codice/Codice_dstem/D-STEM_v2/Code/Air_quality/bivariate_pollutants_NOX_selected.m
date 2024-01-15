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
PM25 = dati_pollutants_2019(6, 2:end);


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

dati_PM25 = dati_PM25(:,14:24:end);
dati_NOX = dati_NOX(:,14:24:end);


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

%% LOOGCV NOX
indici_totali = 1:size(dati_NOX, 1)+size(dati_PM25,1);

rmse_cv = [];
R2_cv = [];
beta_cv = [];
theta_z_cv = [];
v_z_cv = {};
sigma_eta_cv = [];
G_cv = {};
sigma_eps_cv = [];
diag_varcov_cv = {};
log_likelihood_cv = [];

%Setting parametri inziali
beta = [];
theta_z = 0.1;
v_z = 0.2*eye(2);
sigma_eta = diag([1 1]);
G = 0.9*eye(2);
sigma_eps = diag([0.1 0.1]); 

totali_lat = [NOX{1,1}{:,3}; PM25{1,1}{:,3}];
totali_long = [NOX{1,1}{:,4}; PM25{1,1}{:,4}];


for l = 1:size(dati_NOX, 1)
    indici_righe_test = l;
    indici_righe_train = setdiff(indici_totali, indici_righe_test);

    NOx_lat_train = NOX{1,1}{:,3}(indici_righe_train(1,1:size(dati_NOX, 1)-1));
    NOx_long_train = NOX{1,1}{:,4}(indici_righe_train(1,1:size(dati_NOX, 1)-1));
    NOx_alt_train = NOX{1,1}{:,2}(indici_righe_train(1,1:size(dati_NOX, 1)-1));
    
    NOx_lat_test =  NOX{1,1}{:,3}(indici_righe_test);
    NOx_long_test =  NOX{1,1}{:,4}(indici_righe_test);
    NOx_alt_test =  NOX{1,1}{:,2}(indici_righe_test);
    
    PM25_lat = PM25{1,1}{:,3};
    PM25_long = PM25{1,1}{:,4};
    PM25_alt = PM25{1,1}{:,2};    

    % controllo che la stazione di test non sia in nessun dataset di
    % training
    dati_train_PM25 = dati_PM25;  
    if sum(totali_lat == NOx_lat_test) > 1 & sum(totali_long == NOx_long_test) > 1
        PM25_lat = setdiff(PM25_lat, NOx_lat_test, 'stable');
        PM25_long = setdiff(PM25_long, NOx_long_test, 'stable');        
        righe_da_togliere = PM25{1,1}{:,3} == NOx_lat_test;
        PM25_alt = PM25_alt(not(righe_da_togliere));
        dati_train_PM25 = dati_train_PM25(not(righe_da_togliere),:);
    end

    % Estrazione dati train e test
    dati_train_NOX = dati_NOX(indici_righe_train(1,1:size(dati_NOX, 1)-1), :);
    dati_test_NOX = dati_NOX(indici_righe_test, :);    

    %load no2 obs
    ground.Y{1} = dati_train_NOX;
    ground.Y_name{1} = 'nox';
    n1 = size(ground.Y{1}, 1);
    T = size(ground.Y{1}, 2);
    
    %load pm2.5 obs
    ground.Y{2} = dati_train_PM25;
    ground.Y_name{2} = 'pm25';
    n2 = size(ground.Y{2}, 1);
    
    %matrice [stazioni x numero_covariate x giorni]
    X = zeros(n1, 1, T);
    X_krig = zeros(size(dati_test_NOX, 1), 1, T);
    for i=1:T
        X(:,1,i) = NOx_lat_train;
        X(:,2,i) = NOx_alt_train; 
        X_krig(:,1,i) = NOx_lat_test;
        X_krig(:,2,i) = ones(size(dati_test_NOX, 1),1);
        X_krig(:,3,i) = NOx_alt_test;
    end
    ground.X_beta{1} = X;
    ground.X_beta_name{1} = {'lat', 'alt'};
    ground.X_beta_name_krig{1} = {'lat', 'constant','alt'};
    ground.X_beta_krig{1} = X_krig;
    
    
   
    %matrice [stazioni x numero_covariate x giorni]
    X = zeros(n2, 1, T);
    for i=1:T
        X(:,1,i) = PM25_lat;
        X(:,2,i) = PM25_alt;
        %X_krig(:,2,i) = PM25_lat(indici_righe_test2);
        %X_krig(:,3,i) = PM25_long(indici_righe_test2); 
        %X_krig(:,4,i) = ones(size(dati_test_PM25, 1),1); 
        %X_krig(:,5,i) = PM25_alt(indici_righe_test2); 
    end
    ground.X_beta{2} = X;
    ground.X_beta_name{2} = {'lat', 'alt'};
    ground.X_beta_name_krig{2} = {'lat', 'constant','alt'};
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
    obj_stem_gridlist_p = stem_gridlist();
    
    ground.coordinates{1} = [NOx_lat_train, NOx_long_train];
    ground.coordinates{2} = [PM25_lat, PM25_long];
    
    
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
    beta = obj_stem_model.get_beta0();
    obj_stem_par.beta = beta;
    obj_stem_par.theta_z = theta_z;
    obj_stem_par.v_z = v_z;
    obj_stem_par.sigma_eta = sigma_eta;
    obj_stem_par.G = G;
    obj_stem_par.sigma_eps = sigma_eps; 
    
    obj_stem_model.set_initial_values(obj_stem_par);
    
    %Model estimation
    exit_toll = 0.001;
    max_iterations = 400;
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.max_iterations = max_iterations;
    obj_stem_EM_options.exit_tol_par = exit_toll;
    obj_stem_model.EM_estimate(obj_stem_EM_options);
    obj_stem_model.set_varcov;
    obj_stem_model.set_logL;    
    
    %obj_stem_model.print; 
    
    %d = sqrt(diag(obj_stem_model.stem_EM_result.stem_par.v_z).*eye(2));
    %R = inv(d)*obj_stem_model.stem_EM_result.stem_par.v_z*inv(d);
    
    %% Kriging on validation stations
    
    % KRINGING NOX
    krig_coordinates_NOx = [NOx_lat_test, NOx_long_test];
    
    obj_stem_krig_grid = stem_grid(krig_coordinates_NOx, 'deg', 'sparse','point');
    
    
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
    
    rmse_nox = nanstd(dati_test_NOX - y_hat_nox,1,2);
    
    r2_nox = 1 - nanvar(dati_test_NOX - y_hat_nox,1,2)./nanvar(dati_test_NOX,1,2);
    rmse_tot = mean(rmse_nox);
    r2_tot = mean(r2_nox);

    %concateniamo rmse_cv e R2
    rmse_cv = [rmse_cv rmse_tot];
    R2_cv = [R2_cv r2_tot];

    %aggiornamento parametri iterazione successiva
    beta = obj_stem_model.stem_EM_result.stem_par.beta;
    theta_z = obj_stem_model.stem_EM_result.stem_par.theta_z;
    v_z = obj_stem_model.stem_EM_result.stem_par.v_z;
    sigma_eta = obj_stem_model.stem_EM_result.stem_par.sigma_eta;
    G = obj_stem_model.stem_EM_result.stem_par.G;
    sigma_eps = obj_stem_model.stem_EM_result.stem_par.sigma_eps;

    %salvataggio delle distribuzioni
    
    beta_cv = [beta_cv beta];
    theta_z_cv = [theta_z_cv theta_z];
    v_z_cv{1,l} = v_z;
    sigma_eta_cv = [sigma_eta_cv sigma_eta];
    G_cv{1,l} = G;
    sigma_eps_cv = [sigma_eps_cv sigma_eps];
    diag_varcov_cv{1,l} = diag(obj_stem_model.stem_EM_result.stem_par.varcov);
    log_likelihood_cv = [log_likelihood_cv obj_stem_model.stem_EM_result.logL];
    disp(["CROSS-VALIDATION: Iterazione LOOGCV numero: ", num2str(l)]);

end


%Calcolo delle t_stat
t_stat = zeros(size(beta_cv,1), size(beta_cv, 2));
for c =1:size(beta_cv, 2)
    for r = 1:size(beta_cv,1)
        t_stat(r,c) = abs(beta_cv(r,c)/sqrt(diag_varcov_cv{1,c}(r,1)));
        disp(t_stat(r,c))
    end
end

mean(R2_cv)
mean(rmse_cv)

%media delle t_stat
mean(t_stat, 2)


%% salvataggio in .mat
result_data_biavariate_NOX_selected{1} = beta_cv;
result_data_biavariate_NOX_selected{2} = theta_z_cv;
result_data_biavariate_NOX_selected{3} = v_z_cv;
result_data_biavariate_NOX_selected{4} = sigma_eta_cv;
result_data_biavariate_NOX_selected{5} = G_cv;
result_data_biavariate_NOX_selected{6} = sigma_eps_cv;
result_data_biavariate_NOX_selected{7} = diag_varcov_cv;
result_data_biavariate_NOX_selected{8} = log_likelihood_cv;
result_data_biavariate_NOX_selected{9} = t_stat;

save("result_data_biavariate_NOX_selected.mat", 'result_data_biavariate_NOX_selected')

