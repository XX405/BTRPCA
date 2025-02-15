function v=fun_fval(M, L, E, Z, lambda)
alpha=1e-4;mu=1e-4; t=20;
X_temp=pagesvd(L,'econ');
% v3=TL1L2norm(X_temp,t,alpha);
v1 = lambda*(TL1L2norm(X_temp,t,alpha));
v2 = (mu/2)*norm(M-L-E+Z/mu,'fro')^2;
v = v1+v2;
end