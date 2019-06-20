% Authors: Rasoul Shafipour, Ali Khodabakhsh
% Reference: 
% R. Shafipour, A. Khodabakhsh, G. Mateos, and E. Nikolova. A directed graph Fourier transform with spread frequency components. IEEE Trans. Signal Process., 67(4):946?960, Feb 2019

% For a given adjacency matrix, Adj, this function outputs a stationary solution X
% for the problem maximize_{||x|| = 1} DV(x)
function [X]=f_max_fn(Adj)
    N = size(Adj,1);
    k = 1;
    X0 = randn(N,k);
    X0 = orth(X0);
    opts.record = 0;
    opts.mxitr  = 1000;
    opts.xtol = 1e-7;
    opts.gtol = 1e-7;
    opts.ftol = 1e-10;
    out.tau = 1e-7;
    [X, out]= OptStiefelGBB(X0, @funDV, opts, Adj);
    function [F, G] = funDV(x,  Adj)
            NN = length(x);
            grad = zeros(NN,1);
            for j = 1:NN
                grad(j)= 2 * (  Adj(:,j)' * max(x(j) - x,0) -  Adj(j,:) * max(x - x(j),0)  );
            end
            G = - grad;
            F = -sum(sum(Adj .* max(repmat(x',NN,1) - repmat(x,1,NN),0).^2));
    end
end
