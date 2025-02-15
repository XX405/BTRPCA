function  [X,SigmaY] =  WNNM( Y, C, NSig, m, Iter )
    % [U,SigmaY,V] =   pagesvd(full(Y),'econ');  
    sa = size(Y);
la = length(sa);
faces = prod(sa(3:la));
sU = sa;         
sU(2) = sU(1);
sV = sa;        
sV(1) = sV(2);
U = zeros(sU); 
S = zeros(sa);
V = zeros(sV);
 [U,SigmaY,V] =   pagesvd(full(Y),'econ'); 
    PatNum       = size(Y,2);
    Temp         =   sqrt(max( diag(SigmaY).^2 - PatNum*NSig^2, 0 ));
    for i=1:faces
        W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp + eps );               % Weight vector
       	SigmaX   =  soft(SigmaY, diag(W_Vec));
       	Temp     = diag(SigmaX);
    end
               X =  U*SigmaX*V' + m;     
return;


