
clc;
clear;
close all

%Define input file fluid properties
props_file = 'Fluid properties.xlsx';

Ny = 1;
Nx = 200;
D = ones(Ny,Nx);
phi = 0.26*ones(Ny,Nx);

time_step = 0.01;
del_t = 24*60*60*time_step;
t_max = 1.73281e6; %value of tmax is obtained from fw graphs, and knowing t_breakthrough
steps = t_max/del_t;

disp('This is the validation case. Well inputs and geometry are fixed, and can be modified within the code')
disp('The fluid properties are taken from the input file, "Fluid properties.xlsx"')
disp('The file can be modified for different fluid properties as needed');

%Read in fluid properties from input file
props = xlsread(props_file);

mu1 = props(1); mu2 = props(2); 
B1 = props(3)*ones(Ny,Nx);
B2 = props(4)*ones(Ny,Nx); c1 = props(5); c2 = props(6);
cf = props(7); P1_init = 1e7*ones(Ny,Nx); rho1 = props(9);
rho2 = props(10); g = 10; 

gam1=rho1*g./B1; gam2=rho2*g./B2; 
k_x=1e-12*ones(Ny,Nx);
k_y=ones(Ny,Nx);

del_x=ones(Ny,Nx);
del_y=ones(Ny,Nx);
h = 100*ones(Ny,Nx);
Vp = phi.*del_x.*del_y.*h;

%Getting well properties
q1 = zeros(Ny,Nx);q2 = zeros(Ny,Nx);
q1(1,1) = 100/(24*60*60);

Pwf=zeros(Ny,Nx);
Pwf(1,Nx) = 1e5;
Jl = zeros(Ny,Nx);
Jl(1,Nx)=1.5e-10;

ind = 0;

%read permeability and capillary pressure data, fit functions and plot the results
data = xlsread('Data.xlsx');

S_data = data(:,1);
kr_o = data(:,2);
kr_w = data(:,3);
Pc_data = data(:,4);

S1r = 0.2;
S2r = 0.3;

kr10 = 0.2;
kr20 = 0.8;
n1 = 2; n2 = 2;
g1 = @(x) kr10*((x-S1r)./(1-S1r-S2r)).^n1;
g2 = @(x) kr20*((1-x-S2r)./(1-S1r-S2r)).^n2;
S1_init = 0.2*ones(Ny,Nx);

S2_init = 1-S1_init;
S1_init(h==0) = 0;
S2_init(h==0) = 0;

P = P1_init;
S = S1_init;


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
    ct = 1e-20*ones(Ny,Nx);
 
    kr1 = g1(S);
    kr2 = g2(S);
    
    Pc = zeros(Ny,Nx);
    
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
    
    kr_d = (n1)*(kr10/(1-S2r-S1r))*((S-S1r)./(1-S2r-S1r)).^(n1-1);
    kr_d(S <0.2) =0 ;
                  
    S = S +((del_t./Vp).*(T_x_right.*lamb1_right.*gradP_x_r-T_x_left.*lamb1_left.*gradP_x_l) +(del_t./Vp).*(T_y_top.*lamb1_top.*gradP_y_t-T_y_bot.*lamb1_bot.*gradP_y_b) ...
    -(del_t./Vp).*(T_x_right.*lamb1_right.*gam1.*gradD_x_r-T_x_left.*lamb1_left.*gam1.*gradD_x_l) -(del_t./Vp).*(T_y_top.*lamb1_top.*gam1.*gradD_y_t-T_y_bot.*lamb1_bot.*gam1.*gradD_y_b))./(1-del_t.*Jl.*(Pwf-P).*kr_d./mu1./Vp)+kr1.*del_t.*Jl.*(Pwf-P)./(Vp.*mu1-del_t.*Jl.*kr_d.*(Pwf-P))+(del_t./Vp).*q1;

    S(isnan(S)) = 0.2;

end

P_end = P;
S1_final = S;

plot(S1_final);
hold on;

% Buckley Leverett solution
por=0.26;delx = 1;dely = 1;thickness = 100;
q = 100/(24*60*60);
s = 0.7:-0.01:0.2;

kr1_d = (n1)*(kr10/(1-S2r-S1r))*((s-S1r)./(1-S2r-S1r)).^(n1-1);
kr2_d = -1*(n2)*(kr20/(1-S2r-S1r))*((1-s-S2r)./(1-S2r-S1r)).^(n2-1);
fw = g1(s)/mu1./(g1(s)/mu1+g2(s)/mu2);
fw_d = (kr1_d/mu1)./(g1(s)/mu1+g2(s)/mu2)-((g1(s)/mu1+g2(s)/mu2).^(-2)).*(kr1_d/mu1+kr2_d/mu2).*(g1(s)/mu1);
x_a = zeros(size(s));
bt_index = find(abs(s-0.53)<1e-5);
x_a(:) = q*t_max.*fw_d/(por*dely*thickness);
x_a(s<=0.53) = q*t_max*fw_d(bt_index)/(por*dely*thickness);
plot(x_a,s,'r');
