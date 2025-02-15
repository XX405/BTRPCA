function x = softThres(a, tau)

x = sign(a).* max( abs(a) - tau, 0);
end