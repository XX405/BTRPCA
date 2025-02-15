function L = logL(M)
   tau = 30;%%%

mu = 1e-4; 


    lambda = 1/mu;
    inputImgpath = 'C:\Users\37276\Desktop\input\';

%% 读取图片
img_path_list = dir(strcat(inputImgpath, '*.jpg'));
num = length(img_path_list);
    %lambda = opt.lamda;
    M=(M(:,:,:,num));
    Ms = size(M);
    S=zeros(Ms);
    Z = zeros(Ms);
    Z=(Z(:,:,:,num));
    S=(S(:,:,:,num));
Z=rgb2gray(Z);
S=rgb2gray(S);



    G = M - S + Z/(mu);
    %G = X + lambda*Z;
    [Q,Sigma,R] = svd(G);
    p = 1/2;
    [m,n] = size(G);
    sigma = diag(Sigma);
    delta = zeros(size(sigma));
%      
%% self    
    
    v = (2*lambda*(1-p))^(1.0/(2-p));
    v1 = v + lambda*p*v^(p-1);
    for i = 1:m
        s = sigma(i);
        if s >= v1
            x_ = GST(lambda,p,s);
        else
            x_ = 0;
        end
        tau_ = ((1.0/(2*lambda))*(x_-s)^2 + x_^p)^(1.0/p);
        if tau <= tau_
            delta(i) = s;
        else
            delta(i) = x_;
        end
    end
    
    Delta = zeros(size(Sigma));
    Delta(:, 1:m) = diag(delta);
    L = Q*Delta*R';
end