dx = 50;
nx = 203;
ny = 803;
x=(0:nx-1)*dx;
y=(0:ny-1)*dx;
sim2D=karman2d(x,y,0.1,3000);

vp = 2500;
sigma=200;

vel = vp + sigma*sim2D;

figure(1)
pcolor(vel');
shading flat;
colorbar;


fileID = fopen('randomvel.bin','w');
fwrite(fileID,vel','double');
fclose(fileID);
