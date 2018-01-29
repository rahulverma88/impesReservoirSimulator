clc;
clear;

disp('The fluid properties are taken from the input file, "Fluid properties.xlsx"')
disp('The file can be modified for different fluid properties as needed');
%Read in fluid properties from input file
props = xlsread('Fluid properties.xlsx');

mu1 = props(1); mu2 = props(2); B1 = props(3);
B2 = props(4); c1 = props(5); c2 = props(6);
cf = props(7); P1_init = props(8); rho1 = props(9);
rho2 = props(10); g = 9.8; 
%%
%Read geometry data from file, or generate geometry data from input
ans1 = 1;%input('Enter 1 for reading from a geometry file, or 2 for entering dimensions manually: ');
if (ans1 == 1)
    geomfile = 'Grid.xlsx';%input('Enter file name which has depths at each cell, with file extension: ','s');
    D = xlsread(geomfile);
    D=-1*0.3048*D;
    [Ny Nx] = size(D);
else
    Nx = input('Enter number of cells in x-direction: ');
    Ny = input('Enter number of cells in y-direction: ');
    D = input('Enter thickness (assumed uniform: ')*ones(Ny,Nx);
end
%%
%Getting well properties
well_props = xlsread('Wells.xlsx');
well_xcoord = well_props(:,1);
well_ycoord = well_props(:,2);
water_rates = well_props(:,3);
oil_rates = well_props(:,4);
well_pwf = well_props(:,5);
well_skin = well_props(:,6);
well_radii = well_props(:,7);

%%
%read permeability and capillary pressure data, and plot the results

data = xlsread('Data.xlsx');

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
n1 = coeffvalues(kr1_f);
n2 = coeffvalues(kr2_f);

%{
figure;
plot(kr1_f,S_data,kr_w);
hold on;
plot(kr2_f,S_data,kr_o);

figure;
plot(Pc_f,S_data,Pc_data);
%}
%%
del_t=24*60*60;
t_max=del_t*2;
P1_init=6000*6894*ones(Ny,Nx); 
gam1=rho1*g*ones(Ny,Nx);gam2=rho2*g*ones(Ny,Nx);
ct=10^-20*ones(Ny,Nx);
k_x=10^-15*xlsread('Permeability.xlsx'); k_y=10^-15*xlsread('Permeability.xlsx');

del_x=125*.3048*ones(Ny,Nx);del_y=125*.3048*ones(Ny,Nx);
%del_x=1*ones(Ny,Nx);del_y=1*ones(Ny,Nx);
h=0.3048*xlsread('Thickness.xlsx');
por = 0.26;
Vp = por.*del_x.*del_y.*h;
Pc=zeros(Ny,Nx);
dWOC=14000*.3048;

%Note that the the J values will determine if each well is open or not
q1 = zeros(Ny,Nx);q2 = zeros(Ny,Nx);
q1(sub2ind(size(q1),well_ycoord,well_xcoord)) = water_rates*0.3048^3*5.615/(24*60*60);
q2(sub2ind(size(q1),well_ycoord,well_xcoord)) = oil_rates*0.3048^3*5.615/(24*60*60);
%q1(1,1) = 100/(24*60*60);
Pwf=zeros(Ny,Nx);
%Pwf(1,Nx)=10^5;
Pwf(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_pwf;
s = zeros(Ny,Nx);
s(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_skin;
rw = zeros(Ny,Nx);
rw(sub2ind(size(q1),well_ycoord,well_xcoord)) = well_radii;
Gama = 1.73; CA = 31;
Jl=(2*pi*sqrt(k_x.*k_y).*h)./(0.5*log((4*del_x.*del_y)./(Gama*CA.*(rw.^2)))+1/4+s);
Jl(Pwf == 0) = 0;

S = 0.2*ones(Ny,Nx);
P1_init= P1_init-gam1.*(dWOC-D);
P = P1_init;

w_x = zeros(Ny,Nx);
w_y = zeros(Ny,Nx);
%%
tic
for t = del_t:del_t:t_max

ct = c1*S + c2*(1-S) + cf;    
Pc = reshape(feval(Pc_f,S),size(P));
Pc(S<=0.2) = 10;
Pc(S>=0.2) = 0;
Pc = Pc*6894;
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
%%
%Upstreaming of mobility
%Potential calculation at each cell - Phase 1
phi1 = (P./gam1)-D;

phi1_x_r = diff(phi1,1,2);      %Ny by Nx-1 matrix, get phi(i+1)-phi(i)
phi1_x_r(:,Nx) = phi1_x_r(:,Nx-1);
phi1_y_t = diff(phi1,1,1);      %Ny - 1 by Nx matrix get phi(j+1)-phi(j)
if(Ny~=1)
    phi1_y_t(Ny,:) = phi1_y_t(Ny-1,:);
end

kr1 = g1(n1,S);
kr1(S<=0.2) = 0;
kr1(S>0.8) = 0.2;

kr2 = g2(n2,S);
kr2(S<=0.2) = 0.8;
kr2(S>=0.8) = 0;

lamb_r1 = kr1./mu1;
lamb_r2 = kr2./mu2;

w_x_r = zeros(Ny,Nx);
w_x_l = zeros(Ny,Nx);

w_x_r(phi1_x_r>0) = 0;
w_x_r(phi1_x_r<=0) = 1;

w_x_r(:,Nx) = 1;      %No right cell, so use value of current cell

lamb1_right = w_x_r.*lamb_r1 + (1-w_x_r).*circshift(lamb_r1, [0 -1]);

phi1_x_l = circshift(phi1_x_r, [0 1]);
w_x_l(phi1_x_l>=0) = 1;       %Since it is phi_(i)-phi_(i-1), the values reverse from the lamb_right
w_x_l(phi1_x_l<0) = 0;

w_x_l(:,1) = 1;       %No left cell, so use value of current cell
lamb1_left = w_x_l.*lamb_r1 + (1-w_x_l).*circshift(lamb_r1, [0 1]);
%%
w_y_t = zeros(Ny,Nx);
w_y_b = zeros(Ny,Nx);

w_y_t(phi1_y_t>0) = 0;
w_y_t(phi1_y_t<=0) = 1;
w_y_t(Ny,:) = 1;      %No cell above this, so use current cell

lamb1_top = w_y_t.*lamb_r1 + (1-w_y_t).*circshift(lamb_r1, [-1 0]);
%%
phi1_y_b = circshift(phi1_y_t, [1 0]);
w_y_b(phi1_y_b>=0) = 1;
w_y_b(phi1_y_b<0) = 0;

w_y_b(1,:) = 1;       %No cell below this, so use current cell
lamb1_bot = w_y_b.*lamb_r1 + (1-w_y_b).*circshift(lamb_r1, [1 0]);

%Potential calculation at each cell - Phase 2
phi2 = (P./gam2)-D;

%%

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
%%
phi2_x_l = circshift(phi2_x_r, [0 1]);
w_x_l(phi2_x_l>=0) = 1;       %Since it is phi_(i+1)-phi_i, the values reverse from the lamb_right
w_x_l(phi2_x_l<0) = 0;

w_x_l(:,1) = 1;       %No left cell, so use value of current cell
lamb2_left = w_x_l.*lamb_r2 + (1-w_x_l).*circshift(lamb_r2, [0 1]);
%%
w_y_t = zeros(Ny,Nx);
w_y_t(phi2_y_t>0) = 0;
w_y_t(phi2_y_t<=0) = 1;
w_y_t(Ny,:) = 1;      %No cell above this, so use current cell

lamb2_top = w_y_t.*lamb_r2 + (1-w_y_t).*circshift(lamb_r2, [-1 0]);
%%
phi2_y_b = circshift(phi2_y_t, [1 0]);
w_y_b(phi2_y_b>=0) = 1;
w_y_b(phi2_y_b<0) = 0;

w_y_b(1,:) = 1;       %No cell below this, so use current cell
lamb2_bot = w_y_b.*lamb_r2 + (1-w_y_b).*circshift(lamb_r2, [1 0]);

%%
%Building T matrix
 T = zeros(Nx*Ny,Nx*Ny);
       for i = 1:Nx
            for j = 1:Ny
            k = (i-1)*Ny +j;
                    if(h(j,i) == 0)
                        T(k,k) = 1;
                    else
                  
                        T(k,k) = ct(j,i)*Vp(j,i)+del_t*T_x_right(j,i)*(lamb1_right(j,i)+lamb2_right(j,i))+del_t*T_x_left(j,i)*(lamb1_left(j,i)+lamb2_left(j,i))+del_t*T_y_top(j,i)*(lamb1_top(j,i)+lamb2_top(j,i))+del_t*T_y_bot(j,i)*(lamb1_bot(j,i)+lamb2_bot(j,i)) ...
                            +del_t*Jl(j,i)*lamb_r1(j,i)+del_t*Jl(j,i)*lamb_r2(j,i);
                        if (k > 1) 
                            T(k, k-1) = -del_t*T_x_left(j,i)*(lamb1_left(j,i)+lamb2_left(j,i));
                        end

                        if (k < Nx*Ny)
                            T(k, k+1) = -del_t*T_x_right(j,i)*(lamb1_right(j,i)+lamb2_right(j,i));
                        end

                        if (k-Ny > 0)
                            T(k, k - Ny) = -del_t*T_y_bot(j,i)*(lamb1_bot(j,i)+lamb2_bot(j,i));
                        end

                        if (k+Ny <= Nx*Ny)
                            T(k, k + Ny) = -del_t*T_y_top(j,i)*(lamb1_top(j,i)+lamb2_top(j,i));
                        end
                    end
            end
       end
%%
%Build B
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
%%

B = ct.*Vp.*P  ...
           -del_t*T_x_right.*gradD_x_r.*(lamb1_right.*gam1+lamb2_right.*gam2) ...
            +del_t.*T_x_left.*gradD_x_l.*(lamb1_left.*gam1+lamb2_left.*gam2) ...
            +del_t*T_y_bot.*gradD_y_b.*(lamb1_bot.*gam1+lamb2_bot.*gam2) ...    
            -del_t*T_y_top.*gradD_y_t.*(lamb1_top.*gam1+lamb2_top.*gam2) ...
            +del_t*(T_x_right.*lamb2_right.*gradPc_x_r ...
            -T_x_left.*lamb2_left.*gradPc_x_l...
            +T_y_top.*lamb2_top.*gradPc_y_t ...
            -T_y_bot.*lamb2_bot.*gradPc_y_b)...
            + del_t*B1.*q1 + del_t*B2.*q2 ...
            +del_t*Jl.*lamb_r1.*Pwf+del_t*Jl.*lamb_r2.*(Pwf-Pc);

B(h==0) = 9999;

%B is still a Ny x Nx matrix. Needs to be a column matrix for inversion.

B = reshape(B,[Nx*Ny 1]);
%%
%Inverting T
X = T\B;
P = reshape(X,[Ny Nx]);
%%
%getting new saturations
%kr_d = reshape(differentiate(kr1_f,S),size(S));       %getting derivative of relative permeability from the fitted function
kr_d=.2*2.639/.6*(((S-0.2)/.6).^1.639);
kr_d(S<=0.2) = 0;
gradP_x_r = diff(P,1,2);
gradP_x_r(:,Nx) = gradP_x_r(:,Nx-1);
gradP_x_l = circshift(gradP_x_r, [0 1]);
gradP_x_l(:,1) = gradP_x_l(:,2);

if(Ny ~= 1)
    gradP_y_t = diff(P,1,1);
    gradP_y_t(Ny,:) = gradP_y_t(Ny-1,:);
    gradP_y_b = circshift(gradP_y_t, [1 0]);
    gradP_y_b(1,:) = gradP_y_b(2,:);
else
    gradP_y_t = zeros(Ny,Nx);
    gradP_y_b = zeros(Ny,Nx);
end
%%

S =S+(del_t./Vp).*(T_x_right.*lamb1_right.*gradP_x_r-T_x_left.*lamb1_left.*gradP_x_l ...
                        +T_y_top.*lamb1_top.*gradP_y_t-T_y_bot.*lamb1_bot.*gradP_y_b ...
                    -(T_x_right.*lamb1_right.*gam1.*gradD_x_r-T_x_left.*lamb1_left.*gam1.*gradD_x_l) ...
                    -(T_y_top.*lamb1_top.*gam1.*gradD_y_t-T_y_bot.*lamb1_bot.*gam1.*gradD_y_b) ...
                    )./(1-del_t*Jl.*(Pwf-P).*kr_d./mu1./Vp)+kr1.*del_t.*Jl.*(Pwf-P)./(Vp.*mu1-del_t*Jl.*kr_d.*(Pwf-P)) ...
                    +(del_t./Vp).*q1;
S(isnan(S)) =  0.2;
S(S<0.2) = 0.2;
if(~isreal(S))
    break;
end
end
toc