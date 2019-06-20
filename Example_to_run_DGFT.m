clc; clear; close all
A = [0,1,0,0,0;1,0,1,0,0;0,1,0,1,0;0,0,1,0,1;0,0,0,1,0]; % The given adjacency matrix
num_init = 10; % number of initializations for finding the maximum frequency
[DGFT_frequencies,DGFT_basis] = DGFT(A,10) % outputs the frequencies along with the basis