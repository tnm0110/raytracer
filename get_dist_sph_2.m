function [X_t,Z_t,X_max,T] =get_dist_sph_2(p,wave)

% % % % % %  %% % % % % % % % % % %% %% % % % % % % % % % %% % %  %%% % % 
%  
%  This fuction calculate the raypath and travel time with a given 
%  ray parameter that have been used for the search program. It finds 
%  the velocity interval of that correspondins rayp and calculates X and T.  
%  
%  input  : p   == > rayparameter (s/km)
%           : wave ==> P wave or S wave ('Vp' or 'Vs')
%  Output : X_t == > increment and total raypath segments (degree)
%           Z_t == > depth intervals (km)
%           X_max == > total distance (degree)
%           T   == > increment and total travel time
%
% % % % % %  %% % % % % % % % % % %% %% % % % % % % % % % %% % %  %%% % %

data=load ('ak135.mantle.vmod5');

if strcmp(wave,'Vp')
    Vp=data(:,2);
else 
    Vp=data(:,3);
end


Z=data(:,1);
%Vp=data(:,2); 





%p=1/(Vp(562)-0.01);

[~,idx_tp]=min(abs(Vp(:) - 1/p));
dr=5 ;            % depth increment
r_o=6371;         % radious of the eart in km
z_t= Z(idx_tp);   % turing point depth
r_t=r_o-z_t;      % raidus of the tp
p=r_t*p;

% dealing with fist few layers of identical velocities

if  all(Vp(1:idx_tp) == Vp(1))
            d='cannot trace the ray; all layers velocity are the same';
            disp(d);
            X_t=0;
            Z_t=0;
            X_max=0;
            %break;

     elseif ismember(Vp(1:idx_tp-1),1/p) == 1
            %k=1;
            %while Vp(i) ~= Vp(i-k)
            d='identical velocity exists in the upper layer; try different depth';
            disp(d);
            X_t=0;
            Z_t=0;
            X_max=0;
            %break;
    else
        R1=[];
        for i=1:idx_tp-1
           % delta calculation loop
            r=r_o-(i)*dr;
            zeta=r/Vp(i);
            X1(i,1)=rad2deg((p*dr)/(r*(sqrt(zeta^2 - p^2))));
            %X1
            R1=[R1;r];

            X1(i,2)=sum(X1(:,1));
            
        end
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

% calculate travel time
           for j=1:idx_tp-1

            r=r_o-(j)*dr;
            zeta=r/Vp(j);
            T1(j,1)= deg2rad(X1(j,1))*p + (sqrt(zeta^2 - p^2)*dr)/r;
            T1(j,2)=sum(T1(:,1));
           end

        % the mirror part
        t=length(T1);
        T1_s=sum(T1(:,1));
        c=0;
        for j=1:t
            T2(j,1)=T1(t-c,1);
            T2(j,2)=T1_s+sum(T2(:,1));
            c=c+1;
        end
         T=[T1;T2];

end