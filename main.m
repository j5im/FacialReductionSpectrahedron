clear;

addpath(genpath('.\data'))
addpath(genpath('.\src'))

% load the data; Gamma, gamma for Gamma(X) = gamma

% example 1
load data1.mat;   % real case 

% example 2
%load data2.mat;  % complex case




opt = [];
opt.verbose = 1;
 % verbose =0 for no display; =1 for minimal display; =2 full display

%%% Facial reduction for the set 
%%% {rho : Gamma(rho)=gamma, rho psd}, the usual spectrahedron
%%%% MGamma is a matrix representation of the linear map Gamma

% Call the FR routine  (facial reduction Gauss Newton Predictor Corrector)
% Note : the data matrices in Gamma must be linearly independent
[V,Out] = FRGNPC(Gamma,gamma,opt);
% V: Facial range vector 
% Out : flag = -2 : inconclusive; -1 : infeasible; 0 : Slater holds; 
%               1 : reducible (obtain a nontrivial V)
%       (.W, .y, .lambda, .S) = final iterates
       
