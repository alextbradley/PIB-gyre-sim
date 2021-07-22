clc; close all; clear all;

%Plot the barotropic stream function. Includes the redblue colormap.
%brings X, Y, Z, U, V, bsf, topo, and bathy into global scope
% load('two_gyre_sim.mat');
% U2=U; V2=V;
% load('two_gyre_sim_2.mat');
load('gyre_train_data_shiC1e-3_re.mat');
Y = Y - 40*1e3; %centre Y
% clear U V
melt=double(melt); U=double(U); V=double(V);
aa=find(melt==0); melt(aa)=NaN;
bb=find(U==0); U(bb)=NaN;
cc=find(V==0); V(cc)=NaN;
bsf(bb)=NaN;

%%
figure('position',[200 200 700 500]); hold on
% figure('position',[200 200 900 900]); hold on
%contour plot the bsf
contourf(X/1e3, Y/1e3, bsf', 20, 'linestyle', 'none');
colormap(redblue);
c = colorbar;
caxis([-0.5 0.5])
c.Label.String = 'Barotropic Stream Function [sv]';

%add the grounding line and ice front
contour(X/1e3,Y/1e3, topo', [0,0], 'k', 'linewidth', 1.5)

% plot([100 100],[-40 40],'--m','linewidth',2)
% plot([280 280],[-40 40],'--m','linewidth',2)
% plot([100 280],[-40 -40],'--m','linewidth',2)
% plot([100 280],[40 40],'--m','linewidth',2)

%tidy the plot
box on
xlabel('X (km)');
ylabel('Y (km)');
set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[100:25:300],'xticklabel',{'100','','150','','200','','250','','300'},'fontsize',15)
ylim([-40 40])
xlim([100, 280]); %remove left of the grounding line, which is the ice shelf
% set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[0:100:1000],'xticklabel',{'0','','200','','400','','600','','800','','1000'},'fontsize',13)
% ylim([-40 40])
% xlim([0 1000]); %remove left of the grounding line, which is the ice shelf

% [Xi Yi]=meshgrid(X,Y);
% quiver(Xi',Yi',U,V,'k');
axis equal
set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',15)
% set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',13)
set(gcf,'color','w')

print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_barotropic_stream_function.pdf
% print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_barotropic_stream_function_full.pdf

%%
figure('position',[200 200 700 500]); hold on
% figure('position',[200 200 900 900]); hold on
contourf(X/1e3, Y/1e3, bsf', 20, 'linestyle', 'none');
hold on
%tidy the plot
colormap(redblue);
c = colorbar;
caxis([-0.5 0.5])
c.Label.String = 'Barotropic Stream Function [sv]';

[Xi Yi]=meshgrid(X/1e3,Y/1e3);
U(86,39)=-0.05;  V(86,39)=0; % reference = -0.05 
quiver(Xi(1:2:end,1:5:end)',Yi(1:2:end,1:5:end)',U(1:5:end,1:2:end),V(1:5:end,1:2:end),'k','linewidth',0.8,'autoscalefactor',0.6);
box on

% plot([100 100],[-40 40],'--m','linewidth',2)
% plot([280 280],[-40 40],'--m','linewidth',2)
% plot([100 280],[-40 -40],'--m','linewidth',2)
% plot([100 280],[40 40],'--m','linewidth',2)

xlabel('X (km)');
ylabel('Y (km)');

set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[100:25:300],'xticklabel',{'100','','150','','200','','250','','300'},'fontsize',15)
ylim([-40 40])
xlim([100, 280]); %remove left of the grounding line, which is the ice shelf
text(155,35,'5 cm/s','fontsize',12)
% set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[0:100:1000],'xticklabel',{'0','','200','','400','','600','','800','','1000'},'fontsize',13)
% ylim([-40 40])
% xlim([0 1000]); %remove left of the grounding line, which is the ice shelf
% text(155,35,'5 cm/s','fontsize',12)

axis equal
set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',15)
% set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',13)
set(gcf,'color','w')

print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_UV.pdf
% print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_UV_full.pdf

%%
figure('position',[200 200 700 500]); hold on
% figure('position',[200 200 900 900]); hold on
%contour plot the bsf
contourf(X/1e3, Y/1e3, melt', 20, 'linestyle', 'none');
% colormap(redblue);
colormap(jet(55))
caxis([0 55])
c = colorbar;
c.Label.String = 'Melt Rate [m/yr]';

%add the grounding line and ice front
contour(X/1e3,Y/1e3, topo', [0,0], 'k', 'linewidth', 1.5)

%tidy the plot
box on
xlabel('X (km)');
ylabel('Y (km)');
set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[100:25:300],'xticklabel',{'100','','150','','200','','250','','300'},'fontsize',15)
ylim([-40 40])
xlim([100, 280]); %remove left of the grounding line, which is the ice shelf
% set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[0:100:1000],'xticklabel',{'0','','200','','400','','600','','800','','1000'},'fontsize',13)
% ylim([-40 40])
% xlim([0 1000]); %remove left of the grounding line, which is the ice shelf

% [Xi Yi]=meshgrid(X,Y);
% quiver(Xi',Yi',U,V,'k');
axis equal
set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',15)
% set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',13)
set(gcf,'color','w')

print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_melt_rate.pdf
% print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_melt_rate_full.pdf

%%
figure('position',[200 200 700 500]); hold on
% figure('position',[200 200 900 900]); hold on
%contour plot the bsf
ag=find(bathy<0); az=find(bathy==0);
bathy(ag)=NaN;
bathy(az)=45;
contourf(X/1e3,Y/1e3, bathy',20);
% colormap(redblue);
colormap(gray(55))
caxis([0 55])
c = colorbar;
c.Label.String = 'Melt Rate [m/yr]';
contour(X/1e3,Y/1e3, topo', [0,0], 'k', 'linewidth', 1.5)

%tidy the plot
box on
xlabel('X (km)');
ylabel('Y (km)');
set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[100:25:300],'xticklabel',{'100','','150','','200','','250','','300'},'fontsize',15)
ylim([-40 40])
xlim([100, 280]); %remove left of the grounding line, which is the ice shelf
% set(gca,'ytick',[-40:10:40],'yticklabel',{'40','','20','','0','','-20','','-40'},'xtick',[0:100:1000],'xticklabel',{'0','','200','','400','','600','','800','','1000'},'fontsize',13)
% ylim([-40 40])
% xlim([0 1000]); %remove left of the grounding line, which is the ice shelf

% [Xi Yi]=meshgrid(X,Y);
% quiver(Xi',Yi',U,V,'k');
axis equal
set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',15)
% set(gca,'tickdir','out','ydir','reverse','xdir','reverse','fontsize',13)
set(gcf,'color','w')

print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_bathy.pdf
% print -painters -dpdf -r300 D:\Dropbox\STYOON_KOPRI\ANA09B\backup\Paper13_KOPRI\figure\New\LM_Train_Alex_bathy_full.pdf

%% mean meltrate
mm=nanmean(nanmean(melt))

%% Redblue function
function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 
end