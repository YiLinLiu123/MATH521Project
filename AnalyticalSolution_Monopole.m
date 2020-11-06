%% The purpose is to try and plot the exact solutions to the monopole, dipole and quadrapole distributions
x = [-23:-1,1:23];
y = [-23:-1,1:23];
x = x;
y = y;

c0 = 343;
dx = 1; 
rho0 = 1.225; 
w = 2*pi/30;
%from https://www.engineeringtoolbox.com/sound-speed-water-d_598.html at 20 degrees celsius

[X,Y] = meshgrid(x,y);

% define several functions for the distributions, lets just start with
% monopole, so source terms looks about right
S = @(x1,x2) 0.5.*exp(-log(2)/2.*(x1.^2+x2.^2));

source = S(X,Y);
surf(X,Y,source)
xlabel("x coordinate");
ylabel("y coordinate");
zlabel(" S value");
title("Surface Plot of Function S");
%%
% now lets try evualting the dG/dt for monpole source 
%G = @(x1,x2,t,w,c) besselh(0,w/c.*sqrt(x1.^2+x2.^2)) %
G = @(x1,x2,t,w,c) (w/(4*c^2)).*exp(-1i*w.*t).*besselh(0,w/c.*sqrt(x1.^2+x2.^2)).*(-1i*w);
greenFunc = G(X,Y,40*0.5,w,c0);
surf(X,Y,abs(greenFunc))
xlabel("x coordinate");
ylabel("y coordinate");
zlabel(" G value");
title("Surface Plot of Function G");

%% 

% now lets do the convolution for the exact solution

monopole = conv2(source, greenFunc,"same");
surf(X,Y,abs(monopole))
xlabel("x coordinate");
ylabel("y coordinate");
zlabel("monpole value");
axis square
title("Surface Plot of Function monpole");


%% 
% lets try to draw a polar plot
% so we need to get magnitude and angles
%}