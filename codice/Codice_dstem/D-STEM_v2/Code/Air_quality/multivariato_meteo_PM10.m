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

VEL_VENTO = dati_meteo_2019(2, 2:end);
UMIDITA = dati_meteo_2019(5, 2:end);
TEMPERATURA = dati_meteo_2019(4, 2:end);
PRESSIONE = dati_meteo_2019(6, 2:end);


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

%% Rimozione del doppio sensore su una sola stazione
for mese = 1:12
    TEMPERATURA{1,mese}(end, :) = [];
    TEMPERATURA{1,mese}(end-2, :) = [];
    
    UMIDITA{1,mese}(end, :) = [];
    UMIDITA{1,mese}(end-2, :) = [];
    
end

dati_NOX = [];
dati_PM25 = [];
dati_VELVENTO = [];
dati_UMIDITA = [];
dati_TEMPERATURA = [];
dati_PRESSIONE = [];


for i = 1:size(NOX, 2)
    % NOX
    tabella_corrente = NOX{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_NOX = [dati_NOX colonne_numeriche];
    % pm2.5
    tabella_corrente = PM25{1, i};
    colonne_numeriche = tabella_corrente{:, 5:end};
    dati_PM25 = [dati_PM25 colonne_numeriche];    
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

dati_PM25 = dati_PM25(:,14:24:end);
dati_NOX = dati_NOX(:,14:24:end);
dati_VELVENTO = dati_VELVENTO(:,14:24:end);  
dati_UMIDITA= dati_UMIDITA(:,14:24:end) ;
dati_TEMPERATURA = dati_TEMPERATURA(:,14:24:end);  
dati_PRESSIONE = dati_PRESSIONE(:,14:24:end);

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


% LOOGCV NOX
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
beta = [-0.1, 0.3, -0.6, 0.2, -0.2, -0.0, -0.0, -0.0, -0.0, -0.1, -0.1, -0.1, 0.0, 0.1, -0.4, -0.2, 0.3, -0.0, -0.1, -0.4];
theta_z = 0.65;
v_z = [
    1.4  0.2  0.4  0.1 -0.3  0.1;
    0.2  1.3  0.2  0.2 -0.4  0.1;
    0.4  0.2  0.5 -0.2 -0.1  0.2;
    0.1  0.2 -0.2  1.3  0.0 -0.0;
   -0.3 -0.4 -0.1  0.0  1.3  0.2;
    0.1  0.1  0.2 -0.0  0.2  1.1;
];
sigma_eta = 1;
G = eye(6,6).*[0.67, 0.47, 0.89, 0.76, 0.59, 0.70];
sigma_eps = eye(6,6).*[0.13, 0.43, 0.004, 0.04, 0.06, 0.08];

totali_lat = [NOX{1,1}{:,3}; PM25{1,1}{:,3}];
totali_long = [NOX{1,1}{:,4}; PM25{1,1}{:,4}];

for l = size(dati_NOX, 1)+1:size(totali_lat, 1)
    indici_righe_test = l;
    indici_righe_train = setdiff(indici_totali, indici_righe_test);

    NOx_lat = NOX{1,1}{:,3};
    NOx_long = NOX{1,1}{:,4};
    NOx_alt = NOX{1,1}{:,2};    
    
    PM25_lat_train = PM25{1,1}{:,3}(indici_righe_train(1,size(dati_NOX, 1)+1:end)-ones(size(dati_PM25,1)-1,1)'*size(dati_NOX,1));
    PM25_long_train = PM25{1,1}{:,4}(indici_righe_train(1,size(dati_NOX, 1)+1:end)-ones(size(dati_PM25,1)-1,1)'*size(dati_NOX,1));
    PM25_alt_train = PM25{1,1}{:,2}(indici_righe_train(1,size(dati_NOX, 1)+1:end)-ones(size(dati_PM25,1)-1,1)'*size(dati_NOX,1));
    
    PM25_lat_test =  PM25{1,1}{:,3}(indici_righe_test-size(dati_NOX,1));
    PM25_long_test =  PM25{1,1}{:,4}(indici_righe_test-size(dati_NOX,1));
    PM25_alt_test =  PM25{1,1}{:,2}(indici_righe_test-size(dati_NOX,1));
       
    % controllo che la stazione di test non sia in nessun dataset di
    % training
    dati_train_NOX = dati_NOX;  
    if sum(totali_lat == PM25_lat_test) > 1 & sum(totali_long == PM25_long_test) > 1        
        NOx_lat = setdiff(NOx_lat, PM25_lat_test);
        NOx_long = setdiff(NOx_long, PM25_long_test);        
        righe_da_togliere = NOX{1,1}{:,3} == PM25_lat_test;
        NOx_alt = NOx_alt(not(righe_da_togliere));
        dati_train_NOX = dati_train_NOX(not(righe_da_togliere),:);
    end

    % Estrazione dati train e test
    dati_train_PM25 = dati_PM25(indici_righe_train(1,size(dati_NOX, 1)+1:end)-ones(size(dati_PM25,1)-1,1)'*size(dati_NOX,1), :);
    dati_test_PM25 = dati_PM25(indici_righe_test-size(dati_NOX,1), :);    
     

    %load no2 obs
    ground.Y{1} = dati_train_NOX;
    ground.Y_name{1} = 'nox';
    n1 = size(ground.Y{1}, 1);
    T = size(ground.Y{1}, 2);
    
    %load pm2.5 obs
    ground.Y{2} = dati_train_PM25;
    ground.Y_name{2} = 'pm25';
    n2 = size(ground.Y{2}, 1);
    
    %load temp obs
    ground.Y{3} = dati_TEMPERATURA;
    ground.Y_name{3} = 'temp';
    n3 = size(ground.Y{3}, 1);
    
    %load umidity obs
    ground.Y{4} = dati_UMIDITA;
    ground.Y_name{4} = 'umidity';
    n4 = size(ground.Y{4}, 1);
    
    %load wind velocity obs
    ground.Y{5} = dati_VELVENTO;
    ground.Y_name{5} = 'wind_vel';
    n5 = size(ground.Y{5}, 1);
    
    %load precipitation obs
    ground.Y{6} = dati_PRESSIONE;
    ground.Y_name{6} = 'press';
    n6 = size(ground.Y{6}, 1);
    
    %matrice [stazioni x numero_covariate x giorni]
    X = zeros(n1, 1, T);
    %X_krig = zeros(size(dati_test_NOX, 1), 1, T);
    for i=1:T
        if is_weekend(i) == 0
            %creiamo una matrice n_stazioni x 1
            X(:,1,i) = zeros(n1,1);
            %X_krig(:,1,i) = zeros(size(dati_test_NOX, 1),1);
        else
            X(:,1,i) = ones(n1,1);
            %X_krig(:,1,i) = ones(size(dati_test_NOX, 1),1);
        end 
        X(:,2,i) = NOx_lat;
        X(:,3,i) = NOx_long; 
        X(:,4,i) = NOx_alt; 
        %X_krig(:,2,i) = NOx_lat_test;
        %X_krig(:,3,i) = NOx_long_test;
        %X_krig(:,4,i) = ones(size(dati_test_NOX, 1),1); 
        %X_krig(:,5,i) = NOx_alt_test;
    end
    ground.X_beta{1} = X;
    ground.X_beta_name{1} = {'weekend', 'lat', 'long','alt'};
    ground.X_beta_name_krig{1} = {'weekend', 'lat', 'long', 'constant','alt'};
    ground.X_beta_krig{1} = X;
    
    
    %matrice [stazioni x numero_covariate x giorni]
    X = zeros(n2, 1, T);
    X_krig = zeros(size(dati_test_PM25, 1), 1, T);
    for i=1:T
        if is_weekend(i) == 0
            %creiamo una matrice n_stazioni x 1
            X(:,1,i) = zeros(n2,1);
            X_krig(:,1,i) = zeros(size(dati_test_PM25, 1),1);      
        else
            X(:,1,i) = ones(n2,1);
            X_krig(:,1,i) = ones(size(dati_test_PM25, 1),1);        
        end
        X(:,2,i) = PM25_lat_train;
        X(:,3,i) = PM25_long_train;
        X(:,4,i) = PM25_alt_train;
        X_krig(:,2,i) = PM25_lat_test;
        X_krig(:,3,i) = PM25_long_test; 
        X_krig(:,4,i) = ones(size(dati_test_PM25, 1),1); 
        X_krig(:,5,i) = PM25_alt_test; 
    end
    ground.X_beta{2} = X;
    ground.X_beta_name{2} = {'weekend', 'lat', 'long','alt'};
    ground.X_beta_name_krig{2} = {'weekend', 'lat', 'long', 'constant','alt'};
    ground.X_beta_krig{2} = X_krig;
    
    
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
    ground.X_beta{3} = X;
    ground.X_beta_name{3} = {'lat', 'long','alt'};
    ground.X_beta_krig{3} = X;
    
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
    ground.X_beta{4} = X;
    ground.X_beta_name{4} = {'lat', 'long','alt'};
    ground.X_beta_krig{4} = X;
    
    
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
    ground.X_beta{5} = X;
    ground.X_beta_name{5} = {'lat', 'long','alt'};
    ground.X_beta_krig{5} = X;
    
    
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
    ground.X_beta{6} = X;
    ground.X_beta_name{6} = {'lat', 'long','alt'};
    ground.X_beta_krig{6} = X;
    
    
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
    
    
    obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
                                    ground.X_beta, ground.X_beta_name, ... 
                                    ground.X_z, ground.X_z_name);
    
    %laod of the station coordinate   
    UMIDITA_lat = UMIDITA{1,1}{:,3};
    UMIDITA_long = UMIDITA{1,1}{:,4};
    
    PRESSIONE_lat = PRESSIONE{1,1}{:,3};
    PRESSIONE_long = PRESSIONE{1,1}{:,4};
    
    VEL_VENTO_lat = VEL_VENTO{1,1}{:,3};
    VEL_VENTO_long = VEL_VENTO{1,1}{:,4};
                                
    TEMPERATURA_lat = TEMPERATURA{1,1}{:,3};
    TEMPERATURA_long = TEMPERATURA{1,1}{:,4};
    
    obj_stem_gridlist_p = stem_gridlist();
    
    ground.coordinates{1} = [NOx_lat, NOx_long];
    ground.coordinates{2} = [PM25_lat_train, PM25_long_train];
    ground.coordinates{3} = [TEMPERATURA_lat, TEMPERATURA_long];
    ground.coordinates{4} = [UMIDITA_lat, UMIDITA_long];
    ground.coordinates{5} = [VEL_VENTO_lat, VEL_VENTO_long];
    ground.coordinates{6} = [PRESSIONE_lat, PRESSIONE_long];
    
    
    obj_stem_grid1 = stem_grid(ground.coordinates{1}, 'deg', 'sparse', 'point');
    obj_stem_grid2 = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
    obj_stem_grid3 = stem_grid(ground.coordinates{3}, 'deg', 'sparse', 'point');
    obj_stem_grid4 = stem_grid(ground.coordinates{4}, 'deg', 'sparse', 'point');
    obj_stem_grid5 = stem_grid(ground.coordinates{5}, 'deg', 'sparse', 'point');
    obj_stem_grid6 = stem_grid(ground.coordinates{6}, 'deg', 'sparse', 'point');
    
    
    obj_stem_gridlist_p.add(obj_stem_grid1);
    obj_stem_gridlist_p.add(obj_stem_grid2);
    obj_stem_gridlist_p.add(obj_stem_grid3);
    obj_stem_gridlist_p.add(obj_stem_grid4);
    obj_stem_gridlist_p.add(obj_stem_grid5);
    obj_stem_gridlist_p.add(obj_stem_grid6);
    
    
    
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
    max_iterations = 200;
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.max_iterations = max_iterations;
    obj_stem_EM_options.exit_tol_par = exit_toll;
    obj_stem_model.EM_estimate(obj_stem_EM_options);
    obj_stem_model.set_varcov;
    obj_stem_model.set_logL;
    
    
    %obj_stem_model.print; 
    
    %d = sqrt(diag(obj_stem_model.stem_par.v_z).*eye(6));
    %R = inv(d)*obj_stem_model.stem_par.v_z*inv(d);  
    
    
    %% Kriging on validation stations
    
    % KRINGING PM25
    krig_coordinates_PM25 = [PM25_lat_test, PM25_long_test];
    
    obj_stem_krig_grid = stem_grid(krig_coordinates_PM25, 'deg', 'sparse','point');
    
    
    obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid, ground.X_beta_krig{1,2}, ground.X_beta_name_krig{1,2});
    obj_stem_krig = stem_krig(obj_stem_model, obj_stem_krig_data);
    
    obj_stem_krig_options = stem_krig_options();
    obj_stem_krig_options.block_size = 1000;
    
    obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);
    
    %calcolo dell'RMSE e R2
    y_hat_pm25 = obj_stem_krig_result{2,1}.y_hat;
    
    % prendiamo le y originali
    rmse_pm25 = [];
    r2_pm25 = [];
    
    rmse_pm25 = nanstd(dati_test_PM25 - y_hat_pm25,1,2);
    
    r2_pm25 = 1 - nanvar(dati_test_PM25 - y_hat_pm25,1,2)./nanvar(dati_test_PM25,1,2);
    rmse_tot = mean(rmse_pm25);
    r2_tot = mean(r2_pm25);

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
    disp("CROSS-VALIDATION: Iterazione LOOGCV numero: ", num2str(l));
end


%Calcolo delle t_stat
t_stat = zeros(size(beta_cv,1), size(beta_cv, 2));
t = size(dati_NOX,1);
for c =1:size(beta_cv, 2)
    for r = 1:size(beta_cv,1)
        t_stat(r,c) = abs(beta_cv(r,c)/sqrt(diag_varcov_cv{1,c+t}(r,1)));
        disp(t_stat(r,c))
    end
end

mean(R2_cv)
mean(rmse_cv)

%media delle t_stat
mean(t_stat, 2)

%% run 6 OUTLIER
mean([R2_cv(1:5) R2_cv(7:end)])
mean([rmse_cv(1:5) rmse_cv(7:end)])
mean([t_stat(:, 1:5) t_stat(:, 7:end)], 2)


%% salvataggio in .mat
result_data_multivariato_meteo_PM10{1} = beta_cv;
result_data_multivariato_meteo_PM10{2} = theta_z_cv;
result_data_multivariato_meteo_PM10{3} = v_z_cv;
result_data_multivariato_meteo_PM10{4} = sigma_eta_cv;
result_data_multivariato_meteo_PM10{5} = G_cv;
result_data_multivariato_meteo_PM10{6} = sigma_eps_cv;
result_data_multivariato_meteo_PM10{7} = diag_varcov_cv;
result_data_multivariato_meteo_PM10{8} = log_likelihood_cv;
result_data_multivariato_meteo_PM10{9} = t_stat;

save("result_data_multivariato_meteo_PM10.mat", 'result_data_multivariato_meteo_PM10')
 
