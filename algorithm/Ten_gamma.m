function X = Ten_gamma(Y,tau)

[n1,n2,n3] = size(Y);
I=eye([n1,n2]);
Y = fft(Y,[],3);
if mod(n3,2)==1
    for i=1:floor((n3+1)/2)
        [U,S,V] = pagesvd(Y(:,:,i),'econ');
        gamma=10*max(S(:));
        I=eye(size(S));
        delta=max(I-S/gamma,0);
        S1=max(S-tau*delta,0);      
        W(:,:,i)=U*S1*V';
    end
    for i=floor((n3+1)/2)+1:n3
        W(:,:,i)=conj(W(:,:,n3-i+2));
    end
else
    for i=1:floor((n3+1)/2)+1
        [U,S,V] = pagesvd(Y(:,:,i),'econ');
        gamma=100*max(S(:));
        % t=1:3;
        % S=S(:,:,3);
        I=eye(size(S));
        delta=max(I-S/gamma,0);
        S1=max(S-tau*delta,0);  
       
        W(:,:,i)=U*S1*V';                                  
    end
    for i=floor((n3+1)/2)+2:n3
        W(:,:,i)=conj(W(:,:,n3-i+2));
    end
end
X = ifft(W,[],3);
