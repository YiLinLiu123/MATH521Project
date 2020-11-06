%% grid definition
x = [-200:-1,1:200];
y = [-200:-1,1:200];
x = 1*x;
y = 1*y;

c0 = 1;%1481;%damping from C
w = 2*pi/50;
timeStep = 0.001;

[X,Y] = meshgrid(x,y);
F1 = @(x1,x2,t,w) 0.01.*exp(-log(2)/5.*(x2.^2)).*sin(w*t).*cos(pi/10.*x1)*1;
G_Dipole = @(x1,x2,t,w,c) (-(1i*w)/(4*c^3)).*exp(-1i*w*t).*besselh(1,(w/c).*sqrt(x1.^2+x2.^2)).*(x1./sqrt(x1.^2+x2.^2));

%% Try convolutions
%need to pad zeros outside 5,5
F1_Data = F1(X,Y,560*timeStep,w);
F1_Data(:,1:195) =0; 
F1_Data(:,205:end)= 0; 



G_Dipole_Data = G_Dipole(X,Y,560*timeStep,w,c0);
G_Dipole_Data_Real = real(G_Dipole_Data);


product = conv2(-1.*F1_Data, G_Dipole_Data,"same");
contour(X,Y,real(product))
xlabel("x coordinate");
ylabel("y coordinate");
zlabel("dipole Data value");
%set(gca,'DataAspectRatio',[1 1 1E-7])
title("Countor Plot of Function Dipole");

%% Try poltting at y at zero as much as possible. 
Y_0_Data = real(product(200,:));
fig = figure()
plot(x,Y_0_Data)
title("Density along y= 1")
xlabel("x axis")
ylabel("Dipole density perturbation")

%{
%% checking the source function
x_1 = -50:50;
y_1 = -50:50;

x_1 = 0.1*x_1;
y_1 = 1*y_1;

[X1,Y1] = meshgrid(x_1,y_1);
F1_Data2 = F1(X1,Y1,560,w);
f3 = figure;
surf(X1,Y1,F1_Data2)

%% check green's function
%first check hankels's first order 
f4 = figure;
x_2 = -1000:1000;
y_2 = -50:50;

x_2 = 0.2*x_2;
y_2 = 1*y_2;
[X2,Y2] = meshgrid(x_2,y_2);

oneDimensionalHankel = besselh(1,x_2);
plot(x_2,real(oneDimensionalHankel))
title("real hankel 1st order")
xlabel("x axis")
ylabel("value")

f5 = figure 
plot(x_2,imag(oneDimensionalHankel))
title("imag hankel 1st order")
xlabel("x axis")
ylabel("value")
%}