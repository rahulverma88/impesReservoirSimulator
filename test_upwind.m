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
%%

Nx = 4; Ny = 3;
P(Ny,Nx) = 0;
P(1,:) = [2;3;3;1];
P(2,:) = [3;5;4;2];
P(3,:) = [3;3;2;2];

P = P*1e6;
%%
D(Ny,Nx) = 0;
D(1,:) = [1;2;2;1];
D(2,:) = [1.5;2.5;3;2];
D(3,:) = [2;1.5;1;1];
D = D*1e2;
%%
gam1 = ones(Ny,Nx);
gam1 = gam1*1000*9.8;
gam2 = ones(Ny,Nx);
gam2 = gam2*800*9.8;
%%
%Upstreaming of mobility
%Potential calculation at each cell - Phase 1
phi1 = (P./gam1)-D;

%phi = phi*(-1);
%%

phi1_x_r = diff(phi1,1,2);      %Ny by Nx-1 matrix, get phi(i+1)-phi(i)
phi1_x_r(:,Nx) = phi1_x_r(:,Nx-1);
phi1_y_t = diff(phi1,1,1);      %Ny - 1 by Nx matrix get phi(j+1)-phi(j)
if(Ny~=1)
    phi1_y_t(Ny,:) = phi1_y_t(Ny-1,:);
end
%%
S = zeros(Ny,Nx);
S(1,:) = [0.6;0.8;0.8;0.5];
S(2,:) = [0.7;0.9;0.8;0.4];
S(3,:) = [0.8;0.7;0.6;0.5];

S = S - 0.1;
Pc = reshape(feval(Pc_f,S),size(S));
%%
mu1 = 1e-3;
mu2 = 5e-3;

kr1 = g1(n1,S);
kr2 = g2(n2,S);
lamb_r1 = kr1./mu1;
lamb_r2 = kr2./mu2;
%%
w_x_r = zeros(Ny,Nx);
w_x_l = zeros(Ny,Nx);

w_x_r(phi1_x_r>0) = 0;
w_x_r(phi1_x_r<=0) = 1;

w_x_r(:,Nx) = 1;      %No right cell, so use value of current cell

lamb1_right = w_x_r.*lamb_r1 + (1-w_x_r).*circshift(lamb_r1, [0 -1]);
%%
phi1_x_l = circshift(phi1_x_r, [0 1]);
w_x_l(phi1_x_l>=0) = 1;       %Since it is phi_(i+1)-phi_i, the values reverse from the lamb_right
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
%%
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
kr_d = reshape(differentiate(kr1_f,S),size(S));       %getting derivative of relative permeability from the fitted function
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

