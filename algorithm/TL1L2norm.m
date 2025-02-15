function val=TL1L2norm(Z,t,alpha)
%compute truncated l_{1-2} metric
Z=sort(abs(Z),'descend');
val=sum(Z(t+1:end))-alpha*sqrt(sum(Z(t+1:end).^2));
end