clear all;close all;clc

load test_data

sigma_pos=0.1;
l_pos=1;
tau_pos=[sigma_pos;l_pos];
kernelFcn_pos=@(x1,x2)rbf(x1,x2,tau_pos);
muFcn_pos=@(x)zeros(size(x,1),size(x,2));

sigma_vel=1;
l_vel=100;
tau_vel=[sigma_vel;l_vel];
kernelFcn_vel=@(x1,x2)rbf(x1,x2,tau_vel);
muFcn_vel=@(x)zeros(size(x,1),size(x,2));


kernelFcn=@(x1,x2)kernelFcn_pos(x1(:,1),x2(:,1)).*kernelFcn_vel(x1(:,2),x2(:,2));
muFcn=@(x)muFcn_pos(x(:,1)).*muFcn_vel(x(:,2));

p=LUT_position;
v=linspace(0,max(abs(velocity)),10)';


c={p,v};
[c{:}]=ndgrid(c{:});
n=length(c);

LUT_Coloumb_friction_nominal=0.1.*ones(length(LUT_position),1);
x2 = reshape(cat(n+1,c{:}),[],n);
for idx=1:length(x2)
    y2(idx,1)=LUTfriction(x2(idx,1),x2(idx,2),LUT_position,LUT_Coloumb_friction_nominal,LUT_viscous_friction);
end
direction=tanh(velocity*100);

nw=5;
mu2_1_approx=zeros(length(x2),1);

for iF=nw+1:size(position,1)
    idxs=iF+((1-nw):0);

    x1=[position(idxs),velocity(idxs).*direction(idxs)];
    y1=friction(idxs).*direction(idxs)   ;

    x1mean=mean(x1);
    [~,i1]=min(abs(p-x1mean(1)));
    [~,i2]=min(abs(v-x1mean(2)));

    idx=(i2-1)*10+i1;
    K11=kernelFcn(x1,x1);
    K12=kernelFcn(x1,x2(idx,:));
    mu1=muFcn(x1);
    mu2=muFcn(x2(idx,:));

    mu2_1_approx(idx,1)=K12'*(pinv(K11)*y1)+mu2;
end
X2=reshape(x2,[10,10,2]);
P=X2(:,:,1);
V=X2(:,:,2);
Y2=reshape(y2,10,10);
Y2s=reshape(mu2_1_approx,10,10);
%%
figure(1)

nfig=10
for ifig=0:nfig
    Y2t=(ifig*Y2s+(nfig-ifig)*Y2)/nfig;
%surf 1
surf(P,V,Y2,'FaceLighting','gouraud',...
    'MeshStyle','column',...
    'SpecularColorReflectance',0,...
    'SpecularExponent',5,...
    'SpecularStrength',0.2,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'FaceAlpha',0.2,...
    'FaceColor',[0.07 0.6 1],...
    'EdgeAlpha',0.2);
hold on
%surf 2
surf(P,V,Y2t,'SpecularExponent',1,...
    'SpecularStrength',1,...
    'DiffuseStrength',1,...
    'AmbientStrength',0.4,...
    'FaceColor',[0.5 0.5 .5],...
    'AlignVertexCenters','on',...
    'LineWidth',0.2,...
    'EdgeAlpha',1);
hold off


xlabel('Cabin Position')
ylabel('Cabin Speed')
zlabel('Friction')
legend('Nominal','Estimated')

export_fig(sprintf('friction%d.png',ifig),'-r 600','-transparent')

end

%%
figure(2)
plot([y2 mu2_1_approx])
%%
K22=kernelFcn(x2,x2);
testing=[p 4*ones(length(p),1)];
K2t=kernelFcn(x2,testing);
mu_testing=muFcn(testing);
mu_testing_posteriori=K2t'*(pinv(K22)*mu2_1_approx)+mu_testing;
plot(p,mu_testing_posteriori)