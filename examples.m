clear all
clf

% This file has four examples that uses the algorithm newtonR. 

s = 1; % 1, 2, 3 or 4

switch s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
% Example 1: Magnus' linear case
L = 10;
n = 25;
d = L/n;
xc = 0:d:L;
yc = 0:d:L;
zc = 0:d:L;
[xx,yy,zz] = meshgrid(xc,yc,zc);
rho   =   xx - 4;
K     =   xx + yy - 8;
mu    =   xx - yy + 4*zz - 20;
rho_g = 0;
K_g   = 0;
mu_g  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
% Example 2: Magnus' ellipsoid case
L = 10;
n = 25;
d = L/n;
xc = 0:d:L;
yc = 0:d:L;
zc = 0:d:L;
[xx,yy,zz] = meshgrid(xc,yc,zc);
rho   =   1/15*(   (xx-4.01).^2 +   (yy-5.04).^2 + 6*(zz-6.07).^2) -1;
K     =   1/15*(   (xx-4.02).^2 + 6*(yy-5.05).^2 +   (zz-6.08).^2) -1;
mu    =   1/15*( 6*(xx-4.03).^2 +   (yy-5.06).^2 +   (zz-6.09).^2) -1;
rho_g = 0;
K_g   = 0;
mu_g  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
% Example 3: Erling's linear test case 
filename1 = fullfile(pwd, 'dat', 'LinX_TestModel.dat');
filename2 = fullfile(pwd, 'dat', 'LinY_TestModel.dat');
filename3 = fullfile(pwd, 'dat', 'LinZ_TestModel.dat');
[Vol1 Specs1] = readVelrockCube(filename1, 0);
[Vol2 Specs2] = readVelrockCube(filename2, 0);
[Vol3 Specs3] = readVelrockCube(filename3, 0);
xc = Specs1.AxisInfo{1}.Values;
yc = Specs1.AxisInfo{2}.Values;
zc = Specs1.AxisInfo{3}.Values;
rho = Vol1;
K   = Vol2;
mu  = Vol3;
rho_g = 42.2;
K_g   = 45.2;
mu_g  = 35.2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
% Example 4: Erling's real case 
syms x real
syms y real
syms z real
filename1 = fullfile(pwd, 'dat', 'Bulk_QGCH_LGB2050R_A_DEM_26.dat');
filename2 = fullfile(pwd, 'dat', 'Shear_QGCH_LGB2050R_A_DEM_26.dat');
filename3 = fullfile(pwd, 'dat', 'Dsty_QGCH_LGB2050R_A_DEM_26.dat');
[Vol1 Specs1] = readVelrockCube(filename1, 0);
[Vol2 Specs2] = readVelrockCube(filename2, 0);
[Vol3 Specs3] = readVelrockCube(filename3, 0);
xc = Specs1.AxisInfo{2}.Values; %%%%%%%%%%
yc = Specs1.AxisInfo{1}.Values; %% note %%
zc = Specs1.AxisInfo{3}.Values; %%%%%%%%%%
rho = Vol1;
K   = Vol2;
mu  = Vol3;
rho_g = 20;
K_g   = 15;
mu_g  = 2.45;
end % switch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the three surfaces
trans = 0.35;
[xx,yy,zz] = meshgrid(xc,yc,zc);
patch(isosurface(xx,yy,zz,rho,rho_g),'FaceColor','magenta','EdgeColor','none','FaceAlpha',trans);
patch(isosurface(xx,yy,zz,K  ,K_g  ),'FaceColor','yellow' ,'EdgeColor','none','FaceAlpha',trans);
patch(isosurface(xx,yy,zz,mu ,mu_g ),'FaceColor','cyan'   ,'EdgeColor','none','FaceAlpha',trans);
figure(1)
AZ = -40; EL = 17; 
view(AZ, EL); axis vis3d equal on square
camlight; lighting phong
hold on
drawnow

% finding solutions
points = newtonR(xc,yc,zc,rho,K,mu,rho_g,K_g,mu_g);

% plotting solutions
for q = 1:size(points,1)
x = points(q,1);
y = points(q,2);
z = points(q,3);
plot3(x,y,z,'s');
end




