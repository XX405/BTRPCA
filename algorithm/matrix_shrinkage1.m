function [ Z, trnorm, U ] = matrix_shrinkage1( X, lambda)
% X is a matrix
% Shrinkage operation on a matrix X = U*diag(s)*V'
% Z = U*diag(max(0,s-lambda))*V'
% opt = 'matrix' | 'tensor'

% uncomment the code to use PROPACK to compute partial SVD

global sv;
global tmode;
global use_propack;

n = min(size(X,1),size(X,2));
i = tmode;
if i==1
    sv_local=n;
else
sv_local = sv(i);
end

if n==max(size(n), size(sv_local)) % choosvd( n, sv_local) 
    opt.delta = eps;
    opt.eta = eps;
    [ U, S, V ] = lansvd( X, sv_local, 'L', opt);
else
    [ U, S, V ] = pagesvd( X, 'econ' );     %fprintf( 'sv(%d) = %d \n', i, sv_local );
end
inputImgpath = 'C:\Users\37276\Desktop\input\';

%% 读取图片
img_path_list = dir(strcat(inputImgpath, '*.jpg'));
num = length(img_path_list);
indimg=num;
U=rgb2gray(U(:,:,:,indimg));
V=rgb2gray(V(:,:,:,indimg));
S=rgb2gray(S(:,:,:,indimg));
norm_two = lansvd(X, 1, 'L');
mu = 0.1/sqrt(norm_two); % this one can be tuned
myeps = 1e-6;
%myeps = lambda;
%C = sqrt(size(X,1));
C = lambda;

s = diag(S);

% [tempDiagS,svp]=ClosedWNNM(s,C/mu,myeps);
% Z =  U(:,1:svp)*diag(tempDiagS)*V(:,1:svp)';  

 temp   = (s-myeps).^2-4*(C/mu-myeps*s);
 ind    = find (temp>0);
 svp    = length(ind);
 SigmaX = max(s(ind)-myeps+sqrt(temp(ind)),0)/2;
 Z =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';  


% s = diag(S);
% s = s - lambda;    
% sMax = max(SigmaX);
 posI = SigmaX > 0;   %sum(posI)
% Z = scale_matrix( U(:,posI), s(posI), 1 ) * V(:,posI)';
 trnorm = sum(SigmaX(posI));
  
% % adjust sv
% % svp should be no. of sv's larger than a threshold!!
svp = sum(posI);
if svp < sv_local
    sv_local = min(svp + 1, n);
else
    sv_local = min(svp + round(0.05*n), n);
end
sv(i+1) = sv_local;
end