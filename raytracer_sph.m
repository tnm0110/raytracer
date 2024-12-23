function [X_t,Z_t,T]= raytracer_sph(dist,depth,wave)

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% this the main function to trace the ray path with given distance 
% and depth of and event. It will automattically plot the raypth in both
% spherical and flat earth system.
%  input::
%  dist  ==> distance in degree
%  depth ==> source depth in km 
%  output ::
%  X_t = increment and total X
%  Z_t = depth interval
%  T   = increment and total X travel time
%        ex. [ raytracer_sph(50,100, 'Vp')
% % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % % % %  

%dist=90;            % input distance in degree
%depth=0;          

data =load('ak135.mantle.vmod5');
if strcmp(wave,'Vp')
    Vp=data(:,2);
elseif strcmp(wave, 'Vs')
    Vp=data(:,3);
else 
    disp('Given wave is not coreect, allowed: Vp/Vs')
    disp('Plotting P wave ...')
    Vp=data(:,2);
    wave='Vp';
end

Z=data(:,1);%Vp=data(:,2);
dr=5;               % depth increment

% calculate where limits of rayp to look for 

XX=[];
st=2;
for k=st:length(Z)
    z_t=Z(k);
    p=1/Vp(k);
    x=p;
    [~,~,x(2)]=get_dist_sph(z_t,wave);
    XX=[XX;x];
end

dist_range=XX(:,2);
[~,idx] = (min(abs(dist_range(:) - dist)));

p_max=1/Vp(idx-2);
p_min=1/Vp(idx+2);


% loop to find the correct rayparameter for the given distance
d=0;
p=p_min;
iter=1;
i_max=5000; % usually converges within 1000   
while abs(dist-d) > 0.01 && iter < i_max  
    [~,~,d] = get_dist_sph_2(p,wave);    
    %t_d=t_d+(t_d_max-t_d_min)/10;
    p=p+(p_max-p_min)/100;
iter=iter + 1;
end

%% calculte the correct raypath with the corresponding distance

[X_t,Z_t,~,T] =get_dist_sph_2(p,wave);

% souce depth subtraction 
 if exist("depth",'var')
     depth_idx=round(depth/dr);
     X_t(1:depth_idx,:) = [];
     Z_t(1:depth_idx,:) = [];
     T(1:depth_idx,:) = [];
 end

%% plot the raypath

% Plot the ray on the flat earth 
figure(1)
plot(X_t(:,2),Z_t,'LineWidth',1.5);
set(gca,'Ydir','reverse');
set(gca,'XAxisLocation','top')
hold on
scatter(X_t(1,2),Z_t(1),30,'*','r')
hold on
scatter(X_t(end,2),Z_t(end),30,'v','MarkerEdgeColor','k',...
    'MarkerFaceColor','red')    
xlabel('Distance (X_{1})','FontSize',14,'FontWeight','bold')
ylabel('Depth (X_{3})','FontSize',14,'FontWeight','bold')
hold on
%plotting the boundaries

yline(410,'k--','LineWidth',1.0);
hold on 
yline(660,'k--','LineWidth',1.0);
yline(2900,'k--','LineWidth',1.0);
yline(5150,'k--','LineWidth',1.0);
text(82,320,'410');
text(82,760,'660');
text(45,3000,'CMB');
text(41,5250,'IC-OC Boundary')
title('Ray paths with different take-of angles',...
    'FontSize',16,'FontWeight','bold')
% 

% % plot the ray for the spherical earth
% please ucomment all the lines below to plot the raypath

figure(2)
polarplot(deg2rad(X_t(:,2)),Z_t);
hold on 
polarscatter(deg2rad(X_t(1,2)),Z_t(1),30,'*','r')
hold on
polarscatter(deg2rad(X_t(end,2)),Z_t(end),30,'v','MarkerEdgeColor','k',...
    'MarkerFaceColor','red')
ax=gca;
ax.RDir='r';
ax.ThetaGrid='off';
ax.RTick=[410 660 2900 5150 6371];
ax.ThetaZeroLocation='top';
ax.ThetaDir = 'clockwise';
ax.RColor='k';
ax.RAxisLocation = 120;
ax.RTickLabelRotation = 60;
ax.RLim=[0 6371];
ax.LineWidth=1.5;
tit=['Ray paths with different take-of angles for souce depth ', ...
    num2str(depth), ' km'];
title(tit,'FontSize',14,'FontWeight','bold')
hold on

end