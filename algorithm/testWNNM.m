function  [X] = testWNNM( Y,NSig,C )
% [X,trank] = testWNNM( SigmaY,NSig,Y,C, p )
sa = size(Y);
la = length(sa);
% trank = 0;
[U,SigmaY,V] =   svd(full(Y),'econ');   
for i = 3:2
    Y = fft(Y,[],i);
end
Y = fft(Y,[],3);
[~,~,n3] = size(Y);

for i=1:n3
    % sU = sa; 
    % sV= sa; 
% U = zeros(sU); 
% S = zeros(sa);
% V = zeros(sV);
[U(:,:,i),SigmaY(:,:,i),V(:,:,i)] =   svd(Y(:,:,i),'econ');
 s(:,i) = diag(SigmaY(:,:,i));
end
PatNum       = size(Y,2);
    Temp=sqrt(max(s(:,i).^2 - PatNum*NSig^2, 0 ));
     
s1 = zeros(size(s));
    for j=1:4
        for i=1:n3
             % W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( s(:,i)+eps );  
             W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp+eps );
        % W_Vec    =   (C*sqrt(h*d))./( s(:,i) + eps );
       	% s1(:,i)       =   solve_Lp_w(s(:,i), W_Vec, p);
        s1(:,i)       =   solve_Lp_w(s(:,i), W_Vec);
       	s(:,i)      =   s1(:,i);
        Temp     = diag(SigmaX);
        end
    end
   for i=1:n3
   X(:,:,i) =  U(:,:,i)*diag(s1(:,i))*V(:,:,i)' ;
   end
   X = ifft(X,[],3);