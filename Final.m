
clc;
clear;
close all

%Define input files for grid, depth, porosity, permeability and fluid
%properties
geomfile = 'Grid.xlsx';
por_file = 'Porosity.xlsx';
perm_file = 'Permeability.xlsx';
props_file = 'Fluid properties.xlsx';
thickness_file = 'Thickness.xlsx';
data_file = 'Data.xlsx';
well_props_file = 'Wells.xlsx';

D = xlsread(geomfile);
D=-1*0.3048*D;
[Ny Nx] = size(D);
phi = xlsread(por_file); 

time_step = input('Enter time step in days: ');
del_t = 24*60*60*time_step;
t_max=10000*del_t;
steps = t_max/del_t;

disp('The fluid properties, and initial pressure at water-oil contact are taken from the input file, "Fluid properties.xlsx"')
disp('The file can be modified for different fluid properties as needed');
%Read in fluid properties from input file
props = xlsread(props_file);

mu1 = props(1); mu2 = props(2); B1 = props(3)*ones(Ny,Nx);
B2 = props(4)*ones(Ny,Nx); c1 = props(5); c2 = props(6);
cf = props(7); P1_init = props(8)*ones(Ny,Nx)*6894; rho1 = props(9);
rho2 = props(10); g = 9.8; 

gam1=rho1*g./B1; gam2=rho2*g./B2; 
k_x=xlsread(perm_file)*1e-15;
k_y=xlsread(perm_file)*1e-15;

del_x=125*.3048*ones(Ny,Nx);
del_y=125*.3048*ones(Ny,Nx);

h=0.3048*xlsread(thickness_file);
wat_oil_contact = 14000*0.3048; 

Vp = phi.*del_x.*del_y.*h;


%Getting well properties
well_props = xlsread(well_props_file);
well_xcoord = well_props(:,1);
well_ycoord = well_props(:,2);
water_rates = well_props(:,3);
oil_rates = well_props(:,4);
well_pwf = well_props(:,5);
well_skin = well_props(:,6);
well_radii = well_props(:,7);

q1 = zeros(Ny,Nx);q2 = zeros(Ny,Nx);
q1(sub2ind(size(q1),well_ycoord,well_xcoord)) = water_rates*0.3048^3*5.615/(24*60*60);
q2(sub2ind(size(q1),well_ycoord,well_xcoord)) = oil_rates*0.3048^3*5.615/(24*60*60);
Pwf=zeros(Ny,Nx);
Pwf(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_pwf;
s = zeros(Ny,Nx);
s(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_skin;
rw = zeros(Ny,Nx);
rw(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_radii;
Gama = 1.73; CA = 31;

num_wells = size(well_xcoord,1);

%some additional variables needed for generating final results
S_avg = zeros(steps,1);
water_prod = zeros(num_wells,steps); oil_prod=zeros(num_wells,steps); Water_cut=zeros(num_wells,steps);

ind = 0;

%read permeability and capillary pressure data, fit functions and plot the results
data = xlsread(data_file);

S_data = data(:,1);
kr_o = data(:,2);
kr_w = data(:,3);
Pc_data = data(:,4);

S1r = S_data(1);
S2r = 1-S_data(end);

kr10 = kr_w(end);
kr20 = kr_o(1);

g1 = @(n1, x) kr10*((x-S1r)./(1-S2r-S1r)).^n1;
g2 = @(n2, x) kr20*((1-x-S2r)./(1-S2r-S1r)).^n2;

kr1_f = fit(S_data,kr_w, g1,'StartPoint',1.5,'Robust','LAR');
kr2_f = fit(S_data,kr_o, g2,'StartPoint',1.5,'Robust','LAR');
Pc_f = fit(S_data, Pc_data,'poly6');
S_Pc_f = fit(Pc_data,S_data,'exp1');

n1 = coeffvalues(kr1_f);
n2 = coeffvalues(kr2_f);

P1_init= P1_init-gam1.*(wat_oil_contact-D);
Pc_init = (gam1-gam2).*(wat_oil_contact-D); %calculate initial distribution of Pc
Pc_init = Pc_init/6894;
S1_init = S_Pc_f(Pc_init); %get saturations for given Pc
%{
for i = 1:Nx
    for j = 1:Ny
        objective = @(satn)Pc_f(satn)*6894-Pc_init(j,i);
        S1_init(j,i) = fzero(objective, 0.2);               %calculate saturation for given capillary pressure from fitted Pc_f function
    end
end
%}

S1_init(S1_init<0.2) = 0.2; %if calculated saturation is less than S1r, set it to S1r
S1_init = reshape(S1_init,[Nx Ny])';
S2_init = 1-S1_init;
S1_init(h==0) = 0;
S2_init(h==0) = 0;

P = P1_init;
S = S1_init;
Init_oil = sum(sum(S2_init.*Vp./B2))/5.615/0.3046^3;
fprintf('\nAmount of oil initially in place is %4.2f million barrels\n',Init_oil/1e6);


% Transmissibility matrix terms
%T_x(i+1/2,j)
term1 = (del_x)./(k_x.*h.*del_y);
term2 = circshift(term1, [0 -1]);

T_x_right = 2.*(term1+term2).^(-1);
T_x_right(:,Nx) = 0;

%T_x(i-1/2,j)
term2 = circshift(term1, [0 1]);
T_x_left = 2.*(term1+term2).^(-1);
T_x_left(:,1) = 0;

%T_y(i,j+1/2)
term1 = (del_y)./(k_y.*h.*del_x);
term2 = circshift(term1, [-1 0]);
T_y_top = 2.*(term1+term2).^(-1);
T_y_top(Ny,:) = 0;

%T_y(i,j-1/2)
term2 = circshift(term1, [1 0]);
T_y_bot = 2.*(term1+term2).^(-1);
T_y_bot(1,:) = 0;

for t = del_t:del_t:t_max        %Main time loop
    ind = ind+1;
    ct = c1*S + c2*(1-S) + cf;
 
    kr1 = g1(n1,S);
    kr1(S<=0.2) = 0;
    kr1(S>0.8) = 0.2;

    kr2 = g2(n2,S);
    kr2(S<=0.2) = 0.8;
    kr2(S>=0.8) = 0;
    
    Pc = reshape(feval(Pc_f,S),size(P));
    Pc(S<=S1r) = 10;
    Pc(S>=(1-S1r)) = 0;
    Pc = Pc*6894;
    
    q1(sub2ind(size(q1),well_ycoord,well_xcoord)) = water_rates*0.3048^3*5.615/(24*60*60);
    q2(sub2ind(size(q1),well_ycoord,well_xcoord)) = oil_rates*0.3048^3*5.615/(24*60*60);
    Pwf(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_pwf;
    Jl=(2*pi*sqrt(k_x.*k_y).*h)./(0.5*log((4*del_x.*del_y)./(Gama*CA.*(rw.^2)))+1/4+s);
    Jl(Pwf == 0) = 0;

    Jl((P+Pc<=Pwf) & Pwf ~= 0) = 0;

    w=abs(Jl.*(kr1./mu1).*(Pwf-P)./B1);
    o=abs(Jl.*(kr2./mu2).*(Pwf-P-Pc)./B2);
    wor = w./(w+o);

    q1(wor>=0.95) = 2000*5.615*.3048^3/(24*60*60);
    Jl(wor>=0.95) = 0;
    Pwf(wor>=0.95) = 0;
    
    if(wor(well_ycoord(end),well_xcoord(end))> 0.95)
        break;
    end
    
    water_prod_cur = abs(Jl.*(kr1./mu1).*(Pwf-P)./B1);
    oil_prod_cur = abs(Jl.*(kr2./mu2).*(Pwf-P-Pc)./B2);
    water_cut_cur = water_prod_cur./(water_prod_cur+oil_prod_cur);
    water_cut_cur(Jl == 0) = 0;
    
    water_prod(:,ind)= water_prod_cur(sub2ind(size(q1),well_ycoord,well_xcoord));
    oil_prod(:,ind)= oil_prod_cur(sub2ind(size(q1),well_ycoord,well_xcoord));
    Water_cut(:,ind) = water_cut_cur(sub2ind(size(q1),well_ycoord,well_xcoord));
    
    
    lamb_r1 = kr1/mu1;
    lamb_r2 = kr2/mu2;
    
    %Upstreaming of mobility
    %Potential calculation at each cell - Phase 1
    phi1 = (P./gam1)-D;
    
    
    phi1_x_r = diff(phi1,1,2);      %Ny by Nx-1 matrix, get phi(i+1)-phi(i)
    phi1_x_r(:,Nx) = phi1_x_r(:,Nx-1);
    phi1_y_t = diff(phi1,1,1);      %Ny - 1 by Nx matrix get phi(j+1)-phi(j)
    if(Ny~=1)
        phi1_y_t(Ny,:) = phi1_y_t(Ny-1,:);
    end


    w_x_r = zeros(Ny,Nx);
    w_x_l = zeros(Ny,Nx);

    w_x_r(phi1_x_r>0) = 0;
    w_x_r(phi1_x_r<=0) = 1;

    w_x_r(:,Nx) = 1;      %No right cell, so use value of current cell

    lamb1_right = w_x_r.*lamb_r1 + (1-w_x_r).*circshift(lamb_r1, [0 -1]);
    
    phi1_x_l = circshift(phi1_x_r, [0 1]);
    w_x_l(phi1_x_l>=0) = 1;       %Since it is phi_(i+1)-phi_i, the values reverse from the lamb_right
    w_x_l(phi1_x_l<0) = 0;

    w_x_l(:,1) = 1;       %No left cell, so use value of current cell
    lamb1_left = w_x_l.*lamb_r1 + (1-w_x_l).*circshift(lamb_r1, [0 1]);

    w_y_t = zeros(Ny,Nx);
    w_y_b = zeros(Ny,Nx);

    w_y_t(phi1_y_t>0) = 0;
    w_y_t(phi1_y_t<=0) = 1;
    w_y_t(Ny,:) = 1;      %No cell above this, so use current cell

    lamb1_top = w_y_t.*lamb_r1 + (1-w_y_t).*circshift(lamb_r1, [-1 0]);
  
    phi1_y_b = circshift(phi1_y_t, [1 0]);
    w_y_b(phi1_y_b>=0) = 1;
    w_y_b(phi1_y_b<0) = 0;

    w_y_b(1,:) = 1;       %No cell below this, so use current cell
    lamb1_bot = w_y_b.*lamb_r1 + (1-w_y_b).*circshift(lamb_r1, [1 0]);
    
    %Potential calculation at each cell - Phase 2
    phi2 = (P./gam2)-D;

    phi2_x_r = diff(phi2,1,2);      %Ny by Nx-1 matrix, get phi(i+1)-phi(i)
    phi2_x_r(:,Nx) = phi2_x_r(:,Nx-1);
    phi2_y_t = diff(phi2,1,1);      %Ny - 1 by Nx matrix get phi(j+1)-phi(j)
    if(Ny~=1)
        phi2_y_t(Ny,:) = phi2_y_t(Ny-1,:);
    end
    w_x_r = zeros(Ny,Nx);
    w_x_r(phi2_x_r>0) = 0;
    w_x_r(phi2_x_r<=0) = 1;

    w_x_r(:,Nx) = 1;      %No right cell, so use value of current cell

    lamb2_right = w_x_r.*lamb_r2 + (1-w_x_r).*circshift(lamb_r2, [0 -1]);
    
    phi2_x_l = circshift(phi2_x_r, [0 1]);
    w_x_l(phi2_x_l>=0) = 1;       %Since it is phi_(i+1)-phi_i, the values reverse from the lamb_right
    w_x_l(phi2_x_l<0) = 0;

    w_x_l(:,1) = 1;       %No left cell, so use value of current cell
    lamb2_left = w_x_l.*lamb_r2 + (1-w_x_l).*circshift(lamb_r2, [0 1]);
    
    w_y_t = zeros(Ny,Nx);
    w_y_t(phi2_y_t>0) = 0;
    w_y_t(phi2_y_t<=0) = 1;
    w_y_t(Ny,:) = 1;      %No cell above this, so use current cell

    lamb2_top = w_y_t.*lamb_r2 + (1-w_y_t).*circshift(lamb_r2, [-1 0]);
    
    phi2_y_b = circshift(phi2_y_t, [1 0]);
    w_y_b(phi2_y_b>=0) = 1;
    w_y_b(phi2_y_b<0) = 0;

    w_y_b(1,:) = 1;       %No cell below this, so use current cell
    lamb2_bot = w_y_b.*lamb_r2 + (1-w_y_b).*circshift(lamb_r2, [1 0]);
    
    gradD_x_r = diff(D,1,2);
    gradD_x_r(:,Nx) = gradD_x_r(:,Nx-1);
    gradPc_x_r = diff(Pc,1,2);
    gradPc_x_r(:,Nx) = gradPc_x_r(:,Nx-1);

    gradD_x_l = circshift(gradD_x_r, [0 1]);
    gradD_x_l(:,1) = gradD_x_l(:,2);
    gradPc_x_l = circshift(gradPc_x_r, [0 1]);
    gradPc_x_l(:,1) = gradPc_x_l(:,2);

    if(Ny~=1)
        gradD_y_t = diff(D,1,1);
        gradD_y_t(Ny,:) = gradD_y_t(Ny-1,:);
        gradPc_y_t = diff(Pc,1,1);
        gradPc_y_t(Ny,:) = gradPc_y_t(Ny-1,:);

        gradD_y_b = circshift(gradD_y_t,[1 0]);
        gradD_y_b(1,:) = gradD_y_b(2,:);
        gradPc_y_b = circshift(gradPc_y_t,[1 0]);
        gradPc_y_b(1,:) = gradPc_y_b(2,:);
    else
        gradD_y_t = zeros(Ny,Nx);
        gradPc_y_t = zeros(Ny,Nx);

        gradD_y_b = zeros(Ny,Nx);
        gradPc_y_b = zeros(Ny,Nx);
    end
    
    T = zeros(Nx*Ny,Nx*Ny);
    
    for i2 = 1:Nx
        for j2 = 1:Ny
            k2 = (j2-1)*Nx +i2;
            if(h(j2,i2) == 0)
                T(k2,k2) = 1;
            else

                T(k2,k2) = ct(j2,i2)*Vp(j2,i2)+del_t*T_x_right(j2,i2)*(lamb1_right(j2,i2)+lamb2_right(j2,i2))+del_t*T_x_left(j2,i2)*(lamb1_left(j2,i2)+lamb2_left(j2,i2))+del_t*T_y_top(j2,i2)*(lamb1_top(j2,i2)+lamb2_top(j2,i2))+del_t*T_y_bot(j2,i2)*(lamb1_bot(j2,i2)+lamb2_bot(j2,i2)) ...
                    +del_t*Jl(j2,i2)*lamb_r1(j2,i2)+del_t*Jl(j2,i2)*lamb_r2(j2,i2);
                if (k2 > 1) 
                    T(k2, k2-1) = -del_t*T_x_left(j2,i2)*(lamb1_left(j2,i2)+lamb2_left(j2,i2));
                end

                if (k2 < Nx*Ny)
                    T(k2, k2+1) = -del_t*T_x_right(j2,i2)*(lamb1_right(j2,i2)+lamb2_right(j2,i2));
                end

                if (k2-Nx > 0)
                    T(k2, k2 - Nx) = -del_t*T_y_bot(j2,i2)*(lamb1_bot(j2,i2)+lamb2_bot(j2,i2));
                end

                if (k2+Nx <= Nx*Ny)
                    T(k2, k2 + Nx) = -del_t*T_y_top(j2,i2)*(lamb1_top(j2,i2)+lamb2_top(j2,i2));
                end
            end
        end
    end
    
    B = ct.*Vp.*P - del_t.*T_x_right.*gradD_x_r.*(lamb1_right.*gam1+lamb2_right.*gam2) + del_t.*T_x_left.*gradD_x_l.*(lamb1_left.*gam1+lamb2_left.*gam2)+del_t.*T_y_bot.*gradD_y_b.*(lamb1_bot.*gam1+lamb2_bot.*gam2) ...    
        -del_t.*T_y_top.*gradD_y_t.*(lamb1_top.*gam1+lamb2_top.*gam2)+del_t.*(T_x_right.*lamb2_right.*gradPc_x_r -T_x_left.*lamb2_left.*gradPc_x_l+T_y_top.*lamb2_top.*gradPc_y_t ...
        -T_y_bot.*lamb2_bot.*gradPc_y_b)+ del_t*B1.*q1 + del_t*B2.*q2 +del_t*Jl.*lamb_r1.*Pwf+del_t*Jl.*lamb_r2.*(Pwf-Pc);
    
    B(h==0) = 1e17;
    B = reshape(B', [Nx*Ny 1]);
    X = T\B;
    P = (reshape(X, [Nx Ny]))';
    
    
    
    gradP_x_r = diff(P,1,2);
    gradP_x_r(:,Nx) = gradP_x_r(:,Nx-1);
    gradP_x_l = circshift(gradP_x_r, [0 1]);
    gradP_x_l(:,1) = gradP_x_l(:,2);
    
    if(Ny~=1)
        gradP_y_t = diff(P,1,1);
        gradP_y_t(Ny,:) = gradP_y_t(Ny-1,:);

        gradP_y_b = circshift(gradP_y_t,[1 0]);
        gradP_y_b(1,:) = gradP_y_b(2,:);
    else
        gradP_y_t = zeros(Ny,Nx);
        
        gradP_y_b = zeros(Ny,Nx);
    end

    S_avg(ind) = sum(sum(S.*Vp))/sum(sum(Vp));
    
    kr_d = (n1-1)*(kr10/(1-S2r-S1r))*((S-S1r)./(1-S2r-S1r)).^(n1-1);
    kr_d(S <0.2) =0 ;
                  
    S = S +((del_t./Vp).*(T_x_right.*lamb1_right.*gradP_x_r-T_x_left.*lamb1_left.*gradP_x_l) +(del_t./Vp).*(T_y_top.*lamb1_top.*gradP_y_t-T_y_bot.*lamb1_bot.*gradP_y_b) ...
    -(del_t./Vp).*(T_x_right.*lamb1_right.*gam1.*gradD_x_r-T_x_left.*lamb1_left.*gam1.*gradD_x_l) -(del_t./Vp).*(T_y_top.*lamb1_top.*gam1.*gradD_y_t-T_y_bot.*lamb1_bot.*gam1.*gradD_y_b))./(1-del_t.*Jl.*(Pwf-P).*kr_d./mu1./Vp)+kr1.*del_t.*Jl.*(Pwf-P)./(Vp.*mu1-del_t.*Jl.*kr_d.*(Pwf-P))+(del_t./Vp).*q1;

    S(isnan(S)) = 0.2;

end

P_end = P;
S1_final = S;

cumul_oil_prod=sum(sum(Vp./B2))*(.8-(1-S_avg))/5.615/.3046^3;
fprintf('\nUltimate recovery efficiency = %4.2f %%\n',max(max(cumul_oil_prod))/(Init_oil)*100);
fprintf('\nProducing life of the field, in days = %5.0f\n',ind*time_step);
fprintf('\nCumulative oil produced over entire life of reservoir is %4.2f million barrels\n',max(max(cumul_oil_prod))/1e6);


P_end(P_end ==-10^17)=NaN;
S1_final(h == 0)=NaN;
P_end(h == 0)=NaN;

x=1:del_x(1,1):Nx*del_x(1,1);
y=1:del_y(1,1):Ny*del_y(1,1);
[X Y]=meshgrid(x,y); 

fig1 = figure;
plot(kr1_f,'r',S_data,kr_w,'r*');
hold on;
plot(kr2_f,'g',S_data,kr_o,'g*');
legend('Water rel perm data','Fitted curve','Oil rel perm data','Fitted curve');
title('Relative permeability curves');
xlim([0 1]);
ylim([0 1]);
xlabel('Saturation');
ylabel('Rel. perm');

fig2 = figure;
plot(Pc_f,S_data,Pc_data);
title('Capillary pressure vs saturation curve');
legend('Pc-S data','Fitted function');
xlabel('Saturation');
ylabel('Capillary pressure');

fig3 = figure;
surf(X,Y,P_end);view(2);
colorbar
title('Water pressure map (Pa)')
xlabel('X (ft)')
ylabel('Y (ft)')
zlabel('Pressure (Pa)')

fig4 = figure;
surf(X,Y,S1_final);view(2);
colorbar
title('Water saturation map')
xlabel('X(ft)')
ylabel('Y(ft)')
zlabel('Saturation')

time=1:time_step:(ind-1)*time_step;

fig5 = figure('Name','Water production rates for each well','NumberTitle','off');
for well_num=2:num_wells
    subplot(2,2,well_num-1)
    name = ['Well number: ' num2str(well_num)];
    plot(time,water_prod(well_num,time)/.3048^3/5.615*24*3600);
    title(name)
    xlabel('Time(days)')    
    ylabel('Prod. rate(STB/day)')   
end

fig6 = figure('Name','Oil Production rates for each well','NumberTitle','off');
for well_num=2:num_wells
    subplot(2,2,well_num-1)
    name = ['Well number: ' num2str(well_num)];
    plot(time,oil_prod(well_num,time)/.3048^3/5.615*24*3600);
    title(name)
    xlabel('Time(days)')    
    ylabel('Prod. rate(STB/day)')   
end

fig7 = figure('Name','Water cuts for each well','NumberTitle','off');
for well_num=2:num_wells
    subplot(2,2,well_num-1)
    name = ['Well number: ' num2str(well_num)];
    plot(time,Water_cut(well_num,time));
    title(name)
    xlabel('Time(days)')    
    ylabel('Water cut')   
end

fig8 = figure;
plot(time,S_avg(time));
title('Average water saturation vs Time')    
xlabel('time(days)')    
ylabel('Average water saturation')   

fig9 = figure;
plot(time,cumul_oil_prod(time));
title('Cumulative oil recovered vs. Time')    
xlabel('Time(days)')    
ylabel('Oil (STB)')

