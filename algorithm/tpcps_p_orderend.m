function [L, S] = tpcps_p_orderend(M, W, varargin)
% p>3
Ms = size(M);
% len=length(Ms);
M_fnorm = norm(M(:),'fro');
if nargin > 3
    lambda = varargin{1};
    kappa = varargin{2};
else
    ivar = sqrt(Ms(3)*max(Ms(1),Ms(2)));
    % for i=3:len
    %     ivar = (Ms(i))^(i-3)*ivar;
    % end
    lambda = 5/ivar;
    if nargin > 2
        kappa = varargin{1};
    else
        kappa =0.5; 
    end
end


% opt.lambda = 40;
% opt.tau = 30;%%%
% opt.p = 1/2;
lambda1=5/sqrt(Ms(3)*max(Ms(1),Ms(2)));
gamma = 1.15;% gamma = 1.15;
mu2 = 0.03/mean(abs(M(:)));%mu2 = 0.03/mean(abs(M(:)));
mu =1e-4; %试一下mu=1e-4
alpha = 1.5;%1.6
mu_max = 1e10;

%初始化
L = zeros(Ms);
E = L;
Z = L;
N = L;

V=zeros(3*prod(Ms), 1);
p1= zeros(3*prod(Ms), 1);
% V2 = zeros(3*prod(Ms), 1);
tol = 1e-7;   

iter = 0;
% iter_max = 1000;
iter_max = 10;
hasConverged = false;

%迭代
while ~hasConverged && iter < iter_max
    iter = iter + 1;

    S = get_S_update;
    L= get_L_update;
    E = get_E_update;
    F = get_F_update;
    Temp1 = M - S - L;
    Temp2 = L - E - W; 
    Temp3=S-F;
    Z = get_Z_update(Temp1);
    N = get_N_update(Temp2);
    mu = min(mu*alpha, mu_max);  
    % V = get_V_update(Temp3); 
     % G= norm(M(:),'fro');
    A= max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'));
    B=max(norm(Temp2(:), 'fro'),norm(Temp3(:),'fro'));
    hasConverged = max([A,B])/M_fnorm < tol;
      % hasConverged = max(norm(Temp1(:), 'fro'), norm(Temp2(:), 'fro'))/M_fnorm < tol;
    if 1
        if iter == 1 || mod(iter, 10) == 0
            % obj = tnnL+lambda*norm(S(:),1);
             obj = lambda*norm(S(:),12); 

             % obj = l112L+lambda*norm(S(:),1);
            err = norm(Temp1(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                ', obj=' num2str(obj) ', err=' num2str(err)]);
        end
    end
end

    function new_L = get_L_update

new_L= Ten_gamma((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);
% [new_L,tnnL] =WNNM((M - S + W + E + (Z - N)/mu)/2, 0.5/mu);
% new_L = Mytest_prox_tnn((M - S + W + E + (Z - N)/mu)/2,0.5/mu);
end


function new_E = get_E_update

%[new_E,~] = My_prox_tnn1(L - W + N/mu, kappa/mu);
 % [new_E,~] = Mytest_prox_tnn(L - W + N/mu, kappa/mu);
 new_E = Ten_gamma(L - W + N/mu, kappa/mu);
end

function new_S = get_S_update
  
% new_S = opt_L(M - L + Z/mu, lambda/mu);
 new_S= solve_l1l2(M - L + Z/mu, lambda/mu);

% new_S = RPCA_HQF(M - L + Z/mu, lambda/mu);

end
    function new_F=get_F_update   %%mu、mu2为惩罚参数
        % new_F=TTV(M,Z,mu,S+V/mu,lambda/mu);
        % function [C,p1] = TTV(M,L,Z,V,mu,p1)
[h,w,d] = size(M);
sizeD = [h,w,d];
Eny_x  = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y  = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z  = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z  =  permute(Eny_z, [3, 1 2]);
denom1 =sqrt( (Eny_x).^2  +(Eny_y).^2 + (Eny_z).^2);
% denom1 =( (Eny_x)  +(Eny_y) + (Eny_z));
 temp_S  = M(:) - L(:) ;
diffT_p1 = diffT3(mu2*p1+V, sizeD);
Q = reshape(diffT_p1 + mu*temp_S+Z(:),sizeD);
% Q = diffT_p1 + mu*(M-L)+Z(:);
new_F= real( ifftn( fftn(Q) ./ (mu2*denom1 + mu) ) );
diff_x2 = diff3(new_F(:), sizeD);
p1 = softThres( diff_x2 - V/mu2, lambda1/mu2 );
V=V- gamma*mu2*(p1- diff_x2)  ;
    end

function new_Z = get_Z_update(Dir)

new_Z = Z + mu*Dir;

end

function new_N = get_N_update(Dir)

new_N = N + mu*Dir;

end
    % function new_V=get_V_update(Temp3)
    %     new_V=V+mu*(-Temp3);
    % end
  % V= gamma*mu2*(p1- diff_x2)  ;
end