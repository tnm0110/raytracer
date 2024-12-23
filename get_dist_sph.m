function [X_t,Z_t,X_max] = get_dist_sph(z_t,wave)

% %  % %%  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%
%    This code gives raypath and total X taking input as a turning points.
% 
%    input  : z_t == > turing point depth (km)
%           : wave ==> P wave or S wave ('Vp' or 'Vs')
%    output : X_t ==>  increment and total raypath segments (degree)
%             Z_t == > depth intervals (km)
%             X_max == > total distance (degree)   
%
%
% % % % % % % % % % %  % %% % %  %% % % % % % %%  %%% % % % % % % %  % %% 

data=load ('ak135.mantle.vmod5');
Z=data(:,1);
%Vp=data(:,2);

if strcmp(wave,'Vp')
    Vp=data(:,2);
else 
    Vp=data(:,3);
end


dr=5; %cell size
r_o=6371; % radious of the eart in km

r_t=r_o-z_t;  % depth of the turing point in km

idx_tp = ceil(z_t/dr);

p=r_t/Vp(idx_tp);         % ray parameter


    if  all(Vp(1:idx_tp) == Vp(1))
            %d='cannot trace the ray; all layers velocity are the same';
            %disp(d)
            X_t=0;
            Z_t=0;
            X_max=0;
            %break;

        elseif ismember(Vp(1:idx_tp-1),1/p) == 1
            %k=1;
            %while Vp(i) ~= Vp(i-k)
            d='identical velocity exists in the upper layer; try different depth';
            disp(d)
            X_t=0;
            Z_t=0;
            X_max=0;
            %break;
    else
        R1=[];
        for i=1:idx_tp-1
            r=r_o-(i)*dr;
            zeta=r/Vp(i);
            X1(i,1)=rad2deg((p*dr)/(r*(sqrt(zeta^2 - p^2))));
            %X1
            R1=[R1;r];

            X1(i,2)=sum(X1(:,1));
            
        end 
        %X1=[s;S];
        % constructing the mirror part
        t=length(X1);
        X1_s=sum(X1(:,1));
        c=0;
        for j=1:length(X1)
            X2(j,1)=X1(t-c,1);
            X2(j,2)=X1_s+sum(X2(:,1));
            c=c+1;
        end
        X_t=[X1;X2];
        Z_t=[Z(1:length(X1));flip(Z(1:length(X1)))];
        X_max=X_t(end);
    end
    %R=[R1;R1];
    %X_t=rad2deg(X_t);
end