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

classdef stem_kalmanfilter_result < handle
    
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
    %p   - dimension of the latent temporal variable z
    %T   - number of time steps
    %N   - the total number of observation sites (points and pixels)
    
    properties
        zk_f   = [];    %[double]     (pxT+1)    the filtered state
        zk_u   = [];    %[double]     (pxT+1)    the updated state
        Pk_f   = [];    %[double]     {T+1}(pxp) variance-covariance matrix of the filtered state
        Pk_u   = [];    %[double]     {T+1}(pxp) variance-covariance matrix of the updated state
        J_last = [];    %[double]     (pxN)      innovation vector at time t=T
        J      = [];    %[double]     {T+1}(pxN) innovation vector from time t=0 to time t=T
        logL   = [];    %[double]     (1x1)      observed-data log-likelihood
    end
    
    methods
        function obj = stem_kalmanfilter_result(zk_f,zk_u,Pk_f,Pk_u,J_last,J,logL)
            %DESCRIPTION: constructor of the class stem_kalmanfilter_result
            %
            %INPUT 
            %zk_f   = [];    -[double]     (pxT+1)    the filtered state
            %zk_u   = [];    -[double]     (pxT+1)    the updated state
            %Pk_f   = [];    -[double]     {T+1}(pxp) variance-covariance matrix of the filtered state
            %Pk_u   = [];    -[double]     {T+1}(pxp) variance-covariance matrix of the updated state
            %J_last = [];    -[double]     (pxN)      innovation vector at time t=T
            %J      = [];    -[double]     {T+1}(pxN) innovation vector from time t=0 to time t=T
            %logL   = [];    -[double]     (1x1)      observed-data log-likelihood
            %
            %OUTPUT
            %obj             - [stem_kalmanfilter_result object]   (1x1) stem_kalmanfilter_result object  
            
            obj.zk_f=zk_f;
            obj.zk_u=zk_u;
            obj.Pk_f=Pk_f;
            obj.Pk_u=Pk_u;
            obj.J_last=J_last;
            obj.J=J;
            obj.logL=logL;
        end
    end
end