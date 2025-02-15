function [L, S] = tpcps_p_order1(M, W, varargin)
% p>3

Ms = size(M);
[m,n,h]=size(M);
len=length(Ms);
M_fnorm = norm(M(:),'fro');

if nargin > 3
    lambda = varargin{1};
    kappa = varargin{2};
else
    ivar = sqrt(Ms(3)*max(Ms(1),Ms(2)));
    for i=4:len
    % for i=4:len
        ivar = (Ms(i))^(i-3)*ivar;
    end
    lambda = 3/ivar;
    if nargin > 2
        kappa = varargin{1};
    else
        kappa = 0.5; 
    end
end
    Z = M;

theta = .3/ sqrt(max(m,n));
mu = 1e-4;
mu_bar = mu*1e7;
rho = 1.05;          % this one can be tuned
par.alpha = 1;
par.theta = theta;
% mu = 1e-4; 
alpha = 1.5;
mu_max = 1e18;

L = zeros(Ms);
E = L;
Z = L;
N = L;

tol = 1e-10;   

iter = 0;
% iter_max = 1000;
iter_max = 10;
hasConverged = false;

while ~hasConverged && iter < iter_max
    iter = iter + 1;
    
    S = get_S_update;
    L= get_L_update;
    E = get_E_update;
    
    Temp1 = M - S - L;
    Temp2 = L - E - W;
    
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    
    mu = min(mu*alpha, mu_max);  
    
    hasConverged = max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'))/M_fnorm < tol;
    if 1
        if iter == 1 || mod(iter, 10) == 0
            obj = lambda*norm(S(:),12);
            err = norm(Temp1(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                ', obj=' num2str(obj) ', err=' num2str(err)]);
        end
    end
end

function [new_L] = get_L_update
 
% new_L = fun_fval (M,L,S,Z,t_now);
   par.mu = mu;
    new_L =Ten_gamma((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);
   % new_L =fun_fval(M, L,Z,par,lambda);
    L = L+mu*(M-new_L-L);
    mu = min(mu*rho, mu_bar);
% new_L = truncatedL1L2((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);
end

function new_E = get_E_update

% new_E = Ten_gamma(L - W + N/mu, kappa/mu); 
new_E = Ten_gamma(L - W + N/mu, kappa/mu); 

end

function new_S = get_S_update

% new_S = prox_l1(M - L + Z/mu, lambda/mu);
new_S = solve_l1l2(M - L + Z/mu, lambda/mu);
end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu*Dir;

end

function new_N = get_N_update(Dir)

new_N = N + mu*Dir;

end

end