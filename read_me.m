%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% READ ME FILE for epidemic network game m-files
% m-files in the folder
% epidemic_game_main.m : main file that is used to generate Figure 6 data
%                        You have to run it with different \beta values and 
%                        initial conditions to get the data for figure 6 
%                        This file uses two helper m files: 
%                        1) iterated_elimination.m: solves Nash equilibrium
%                        2) preferential_attachment_network.m: generates
%                        scale free network.
% epidemic_runs_star_network_demo.m : Runs the demo for the star network in
%                                    Figure 3 of the paper. 
% R_0_computation_new.m : Gets R_0 plot data for Figure 4 of the manuscript. 
% R_x_computation.m : Gets R_* plot data for Figure 5 of the manuscript. 
% 
% manuscript_figures.m : Allows to get figures from data. Either load the
%                        associated .mat files or use the data from your
%                        latest runs to get the figures. 
%
% iterated_elimination.m: Solves for the stage game Nash equilibrium given
%                         state. 
% preferential_attachment_network.m: generates scale free network.

% Summary for generating manuscript figures from scratch.
% Figure 3: 1) Run epidemic_runs_star_network_demo.m
%           2) Open manuscript_figures.m. Run Figure 3 plots code section

% Figure 4: 1) Run R_0_computation_new.m
%           2) Open manuscript_figures.m. Go to Figure 4 plots. Comment out
%             load .mat file. Run the code for Figure 4 plots

% Figure 5: 1) Run R_x_computation.m
%           2) Open manuscript_figures.m. Go to Figure 5 plots. Comment out
%             load .mat file. Run the code for Figure 5 plots

% Figure 6: 1) Run epidemic_game_main.m for (rate = {0.1,0.2,0.3}), all
%              infected or single infected initially
%           2) Open manuscript_figures.m. Go to Figure 6 plots. Comment out
%             load .mat file. Run the corresponding code in Figure 6 plots

% Summary for generating supplementary figures from scratch.
% Figure 3 (supplementary): 1) Run epidemic_runs_star_network_demo.m
%                           2) Open manuscript_figures.m. Comment out
%                              load .mat file. Run Figure 3 (Supplement)
%                              code section.

% Figure 4 (supplementary): 1) Run epidemic_runs_star_network_demo.m
%                           2) Open manuscript_figures.m. Comment out
%                              load .mat file. Run Figure 4 (Supplement)
%                              code section.

% 