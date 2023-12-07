%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D-STEM - Distributed Space Time Expecation Maximization              %
%%%                                                                      %
%%% Author: Francesco Finazzi                                            %
%%% E-mail: francesco.finazzi@unibg.it                                   %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?francesco.finazzi          %
%%%                                                                      %
%%% Author: Yaqiong Wang                                                 %
%%% E-mail: yaqiongwang@pku.edu.cn                                       %
%%% Affiliation: Peking University,                                      %
%%%              Guanghua school of management,                          %
%%%              Business Statistics and Econometrics                    %
%%%                                                                      %
%%% Author: Alessandro Fassò                                             %
%%% E-mail: alessandro.fasso@unibg.it                                    %
%%% Affiliation: University of Bergamo                                   %
%%%              Dept. of Management, Information and                    %
%%%              Production Engineering                                  %
%%% Author website: http://www.unibg.it/pers/?alessandro.fasso           %
%%%                                                                      %
%%% Code website: https://github.com/graspa-group/d-stem                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This file is part of D-STEM.
% 
% D-STEM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 2 of the License, or
% (at your option) any later version.
% 
% D-STEM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with D-STEM. If not, see <http://www.gnu.org/licenses/>.

classdef stem_EM_result < handle
    
    %PROPERTIES
    %Each class property or method property is defined as follows
    %
    %"Name"="Default value";    %["type"]    "dimension"     "description" 
    %
    %DIMENSION NOTATION
    %(1 x 1) is a scalar
    %(N x 1) is a Nx1 vector
    %(N x T) is a NxT matrix
    %(N x B x T) is a NxBxT array
    %{q} is a cell array of length q
    %{q}{p} is a cell array of length q, each cell is a cell array of length p
    %{q}(NxT) is a cell array of length q, each cell is a NxT matrix
    %
    %CONSTANTS
    %N_p = n1_p+...+nq_p - total number of point sites
    %N_b = n1_b+...+nq_b - total number of pixel sites
    %N   = N_p+N_b - total number of observation sites
    %T   - number of temporal steps
    %H   - total number of scalar model parameters
    %I   - number of EM iterations at convergence 
    
    properties
        stem_par=[];                    %[stem_par object]  (1x1) estimated model parameters after the last EM itaration
        stem_kalmansmoother_result=[]   %[stem_kalmansmoother_result object]    (1x1) stem_kalmansmoother_result object
   
        y_hat=[];                       %[double]           {q}(n_i x T) Estimated variable data with missing data filled
        res=[];                         %[double]           {q}(n_i x T) Model residuals
        diag_Var_y_hat_back=[]          %[double]           {q}(n_i x T) diagonals of Var[y|Y(1)] backtrasformed
        E_wp_y1=[];                     %[double]           (N_pxTxK) E[wp|Y(1)]
        E_wb_y1=[];                     %[double]           (N_pxTxK) E[wp_k|Y(1)]
        diag_Var_wp_y1=[];              %[double]           (N_pxTxK) diagonals of Var[wp|Y(1)]
        diag_Var_wb_y1=[];              %[double]           (N_bxT)   diagonals of Var[wb|Y(1)]
        diag_Var_e_y1=[];               %[double]           {q}(n_i x T)  diagonals of Var[y|Y(1)]
        
        y_hat_back=[];                  %[double]           {q}(n_i x T) Backtransformed estimated variable data with missing data filled
        y_back=[];                      %[double]           {q}(n_i x T) Backtransformed original variable data
        res_back=[];                    %[double]           {q}(n_i x T) Backtransformed model residuals
        AIC = [];                       %[double]           (1x1) AIC at convergence 
        logL=[];                        %[double]           (1x1) observed-data log-likelihood at convergence
        logL_all=[];                    %[double]           (Ix1) observed-data log-likelihood at all EM iterations
        
        R2=[];                          %[double >=0]       (dNx1) R-square for in-sample prediction
        
        iterations=[];                  %[integer >0]       (1x1) number of EM iterations at convergence
        max_iterations=[];              %[integer >0]       (1x1) maximum number of EM iterations
        exit_tol_par=[];                %[double] >0]       (1x1) exit tolerance on parameters
        exit_tol_loglike=[];            %[double] >0]       (1x1) exit tolerance on log-likelihood
        partitions=0;                   %[integer>0]        (1x1) partitioning used in model estimation
        
        computation_time=[];            %[double]           (1x1) EM estimation time in seconds
        date_start=[];                  %[double]           (1x1) date of the start of the EM estimation
        machine=[];                     %[string]           (1x1) operative system of the machine
        version=[];                     %[string]           (1x1) the version of D-STEM
    end
    
    methods
        function obj = stem_EM_result()
            if nargin>0
                error('stem_EM_result constructor does not need input arguments');
            end
        end
    end
end