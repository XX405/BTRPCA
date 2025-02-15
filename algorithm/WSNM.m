function  [X,U,V,S] =  WSNM( Y, C, p,h,d )
   [U,S,V] =   takeSVDs(full(Y),'econ');    
    % [U,SigmaY,V] =   pagesvd(full(Y),'econ');    
%     PatNum       = size(Y,2);
%     Temp         =   sqrt(max( diag(SigmaY).^2 - PatNum*NSig^2, 0 ));
    s = diag(S);
    s1 = zeros(size(s));
 
%     for i=1:Iter
%         W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp + eps );               % Weight vector
%        	SigmaX   =  soft(SigmaY, diag(W_Vec));
%        	Temp     = diag(SigmaX);
%     end
%     X =  U*SigmaX*V' + m;     


    for i=1:4
        W_Vec    =   (C*sqrt(h*d))./( s + eps );               % Weight vector
        %W_Vec    =   (C*sqrt(PatNum)*NSig^2)./( Temp.^(1/p) + eps );
       	s1       =   solve_Lp_w(s, W_Vec, p);
       	s     =   s1;
    end
    SigmaX = diag(s1);
    X =  U*SigmaX*V' ;
%      X =  U*SigmaX*V' + m; 
return;