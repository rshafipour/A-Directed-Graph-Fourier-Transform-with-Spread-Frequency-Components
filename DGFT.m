% Authors: Rasoul Shafipour, Ali Khodabakhsh
% Reference: 
% R. Shafipour, A. Khodabakhsh, G. Mateos, and E. Nikolova. A directed graph Fourier transform with spread frequency components. IEEE Trans. Signal Process., 67(4):946?960, Feb 2019
% This function outputs the DGFT frequencies and basis vectorts given an N by N adjacency matrix
% DGFT_frequencies: a 1 by N vector collecting the frequencies (directed variation of basis signals)
% DGFT_basis: an N by N matrix storing the basis vectors as columns
% A: (directed) adjacency matrix
% itr: number of initializations that is used to pick the maximum frequency
function [DGFT_frequencies,DGFT_basis] = DGFT(A,itr)
N = size(A,1);
temp=zeros(N,itr);
temp2=zeros(1,itr);
for i=1:itr
    temp(:,itr) = f_max_fn(A); % finds a stationary solution for the problem maximize_{||x|| = 1} DV(x)
    temp2(itr) = sum(sum(A .* max(repmat(temp(:,itr)',N,1) - repmat(temp(:,itr),1,N),0).^2));
end
[~,index]=max(temp2); % pick a stationary solution that has the maximum variation
last_col=temp(:,index); % last basis vector
first_col = ones(N,1)/sqrt(N); % first basis vector
[DGFT_basis,DGFT_frequencies] = dispersion(A,first_col,last_col); % fixing the first and last basis vectors,
% find the rest of the basis vectors such that the spectral dispersion is minimized
end