load("Data\EQN_dataset\ground.mat");
load("Data\EQN_dataset\station.mat");

addpath('C:/Users/arici/Desktop/SHDD/Finazzi/D-STEM/Src/');

%implementiamo modello multivariato DCM
T = 1; % perchè non abbiamo dati temporali

%'INGV PGA'	'EQN PGA'	'EQN felt intensity'	'EMSC felt intensity'

MSE = [];

y_vera = ground.Y{1};
coord_vera= ground.coordinates{1};


ground.X_p{1} = ground.X_p{1}(1:end-1);

%inserire il ciclo for di indice i

i =1;
y_estratta = y_vera(i);
coord_estratta = coord_vera(i, :);

y_vera_meno1 = [y_vera(1:i-1); y_vera(i+1:end)];
coord_vera_meno1 = [coord_vera(1:i-1,:); coord_vera(i+1:end, :)];
 

ground.Y{1} = y_vera_meno1;

ground.coordinates{1} = coord_vera_meno1;
obj_stem_varset_p = stem_varset(ground.Y, ground.Y_name, [], [], ...
    [], [], ...
    [], [], ...
    ground.X_p,ground.X_p_name);

obj_stem_gridlist_p = stem_gridlist();
    obj_stem_grid = stem_grid(coord_vera_meno1, 'deg', 'sparse', 'point');
    obj_stem_gridlist_p.add(obj_stem_grid);
    obj_stem_grid = stem_grid(ground.coordinates{2}, 'deg', 'sparse', 'point');
    obj_stem_gridlist_p.add(obj_stem_grid);
    obj_stem_grid = stem_grid(ground.coordinates{3}, 'deg', 'sparse', 'point');
    obj_stem_gridlist_p.add(obj_stem_grid);
    obj_stem_grid = stem_grid(ground.coordinates{4}, 'deg', 'sparse', 'point');
    obj_stem_gridlist_p.add(obj_stem_grid);

    %Creare lo shape
    shape = []; %prenderlo dal file del prof

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Model building     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Mettiamo la parte iniziale e finale allo stesso valore perchè tutte le registrazioni che abbiamo non dipendono dal tempo
    obj_stem_datestamp = stem_datestamp('01-01-2009 00:00', '01-01-2009 00:00', T);

    obj_stem_modeltype = stem_modeltype('DCM');
    obj_stem_data = stem_data(obj_stem_varset_p, obj_stem_gridlist_p, ...
        [], [], obj_stem_datestamp, [], obj_stem_modeltype,shape);

    %stem_par object creation
    obj_stem_par_constraints=stem_par_constraints();
    obj_stem_par = stem_par(obj_stem_data, 'exponential',obj_stem_par_constraints);
    %stem_model object creation
    obj_stem_model = stem_model(obj_stem_data, obj_stem_par);
    
    %Data transform
    obj_stem_model.stem_data.log_transform;
    obj_stem_model.stem_data.standardize;

    %Starting values
    %obj_stem_par.beta = b_acc;% [b_acc; b_int];
    obj_stem_par.theta_p = 5;
    obj_stem_par.v_p = eye(4)*0.7;
    obj_stem_par.sigma_eps = eye(4)*0.7;

    obj_stem_model.set_initial_values(obj_stem_par);

    %Model estimation
    obj_stem_EM_options = stem_EM_options();
    obj_stem_EM_options.exit_tol_par=0.005;
    obj_stem_EM_options.max_iterations=50;
    obj_stem_model.EM_estimate(obj_stem_EM_options);
    

    %Kriging
    obj_stem_krig_grid = stem_grid(coord_estratta, 'deg', 'regular','pixel', 1 ,'square',0.001,0.001);
    
    obj_stem_krig_data = stem_krig_data(obj_stem_krig_grid,ones(length(coord_estratta),1),{'constant'});
    obj_stem_krig = stem_krig(obj_stem_model,obj_stem_krig_data);
    
    obj_stem_krig_options = stem_krig_options();
    obj_stem_krig_options.block_size = 1000;
    
    obj_stem_krig_result = obj_stem_krig.kriging(obj_stem_krig_options);
