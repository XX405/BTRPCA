% function [C,p1] = TTV(T,Z,A2,A3,lamda3,mu,p1)
function [C,p1,h,w,d] = TTV(T,B,A3,lamda3,mu,p1)
% [h,w,d] = size(M);
[h,w,d] = size(T);
sizeD = [h,w,d];
Eny_x  = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y  = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z  = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;

Eny_z  =  permute(Eny_z, [3, 1 2]);
denom1 =  Eny_x + Eny_y + Eny_z;
[diffT_p1,A3] = diffT3(mu*p1+A3, sizeD,denom1,Eny_x,Eny_y,Eny_z,T,p1,A3,B,lamda3);
diffT_p1 = reshape(diffT_p1,sizeD);
% Q = diffT_p1 + mu*(T-Z)+A2;
Q = diffT_p1 + mu*(T-B)+T;
C = real( ifftn( fftn(Q) ./ (mu*denom1 + mu) ) );
diff_x2 = diff3(C(:), sizeD);
p1 = softThres( diff_x2 - A3/mu, lamda3/mu );
end

