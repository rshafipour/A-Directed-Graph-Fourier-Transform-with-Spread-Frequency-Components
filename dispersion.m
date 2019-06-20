% Authors: Rasoul Shafipour, Ali Khodabakhsh
% Reference: 
% R. Shafipour, A. Khodabakhsh, G. Mateos, and E. Nikolova. A directed graph Fourier transform with spread frequency components. IEEE Trans. Signal Process., 67(4):946?960, Feb 2019


% fixing the first and last basis vectors, this function
% finds the rest of the basis vectors such that the spectral dispersion is minimized
function [X,DV_final]=dispersion(Adj,first_col,last_col)
    N = size(Adj,1);
    X0 = randn(N,N);
    X0 = orth(X0);
    opts.record = 0;
    opts.mxitr  = 1000;
    opts.xtol = 1e-7;
    opts.gtol = 1e-7;
    opts.ftol = 1e-10;
    out.tau = 1e-7;
    [X, out]= OptStiefelGBB(X0, @funDV, opts, Adj);
    DV_final = zeros(1,N);
    for ii = 1:N
        vecc = X(:,ii);
        DV_final(ii) = sum(sum(Adj .* max(repmat(vecc',N,1) - repmat(vecc,1,N),0).^2));
    end
    function [F, G] = funDV(x,  Adj)
            NN = size(x,1);
            lambda = 10000;
            grad = zeros(NN,NN);
            for i=2:NN-1
                DV_prev = sum(sum(Adj .* max(repmat(x(:,i-1)',NN,1) - repmat(x(:,i-1),1,NN),0).^2));
                DV_curr = sum(sum(Adj .* max(repmat(x(:,i)',NN,1) - repmat(x(:,i),1,NN),0).^2));
                DV_next = sum(sum(Adj .* max(repmat(x(:,i+1)',NN,1) - repmat(x(:,i+1),1,NN),0).^2));
                for j = 1:NN
                    grad(j,i)= 4 * (  Adj(:,j)' * max(x(j,i) - x(:,i),0) -  Adj(j,:) * max(x(:,i) - x(j,i),0)  ) * (2*DV_curr-DV_prev-DV_next);
                end
            end
            % first col of gradient
            DV_curr = sum(sum(Adj .* max(repmat(x(:,1)',NN,1) - repmat(x(:,1),1,NN),0).^2));
            DV_next = sum(sum(Adj .* max(repmat(x(:,2)',NN,1) - repmat(x(:,2),1,NN),0).^2));
            for j = 1:NN
                grad(j,1)= 4 * (  Adj(:,j)' * max(x(j,1) - x(:,1),0) -  Adj(j,:) * max(x(:,1) - x(j,1),0)  ) * (DV_curr-DV_next);
            end
            grad(:,1)=grad(:,1)+lambda*(x(:,1)-first_col);
            % last col of gradient
            DV_curr = sum(sum(Adj .* max(repmat(x(:,NN)',NN,1) - repmat(x(:,NN),1,NN),0).^2));
            DV_prev = sum(sum(Adj .* max(repmat(x(:,NN-1)',NN,1) - repmat(x(:,NN-1),1,NN),0).^2));
            for j = 1:NN
                grad(j,NN)= 4 * (  Adj(:,j)' * max(x(j,NN) - x(:,NN),0) -  Adj(j,:) * max(x(:,NN) - x(j,NN),0)  ) * (DV_curr-DV_prev);
            end
            grad(:,NN)=grad(:,NN)+lambda*(x(:,NN)-last_col);
            G = grad;

            DVs = zeros(1,NN);
            for i = 1:NN
                vec = x(:,i);
                DVs(i) = sum(sum(Adj .* max(repmat(vec',NN,1) - repmat(vec,1,NN),0).^2));
            end
            F = sum(diff(DVs).^2) + lambda*(norm(x(:,1)-first_col)^2 + norm(x(:,NN)-last_col)^2);
    end
end
