function K=rbf(x1,x2,tau)
sigma=tau(1);
l=tau(2);

K=sigma*ones(length(x1),length(x2));
for ic=1:size(x1,2)
K=K.*exp(-1/2*((x1(:,ic)-x2(:,ic)')/l).^2);
end