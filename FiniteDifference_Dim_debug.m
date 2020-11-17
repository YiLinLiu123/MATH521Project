i%% set up the domain including 3 ghost cells
clear
close all
dx = 1; 
dy = 1; 
c = 1; %speed of sound in air
gamma = 1.4; %specific heat ratio
dt = 1/(c);

xRange = -23:dx:23; 
yRange = -23:dy:23;

% constants
u_Mean = 0;
v_Mean = 0; 
rho_Mean = 1.225; 
p_Mean = 1; 
mean_Values = [rho_Mean, u_Mean, v_Mean, p_Mean];

% initialize the matrix
[X,Y]= meshgrid(xRange,yRange);
Un(:,:,1) = zeros(size(X,1),size(X,2));
Un(:,:,2) = zeros(size(X,1),size(X,2));
Un(:,:,3) = zeros(size(X,1),size(X,2));
Un(:,:,4) = zeros(size(X,1),size(X,2));


%source function
epi= 0.5; 
alpha = log(2)/5;
w = 2*pi/30;
S = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t).*[1;0;0;1]; %monopole only  density perturbation
STest = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t); %monopole only for density perturbation
%% coefficients for spatial scheme
%matrix of coefficients
optimized7Centered = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    %website: https://bellaard.com/tools/Finite%20difference%20coefficient%20calculator/
back1 = (1/dx).*[1/60,-2/15,1/2,-4/3,7/12,2/5,-1/30]; % -4,-3,-2,-1,0,1,2
back2 = (1/dx).*[-1/30,1/4,-5/6,5/3,-5/2,77/60,1/6]; % -5,-4,-3,-2,-1,0,1
back3 = (1/dx).*[1/6,-6/5,15/4,-20/3,15/2,-6,49/20];% -6,-5,-4,-3,-2,-1,0
forward1 = (1/dx).*[1/30,-2/5,-7/12,4/3,-1/2,2/15,-1/60]; %-2,-1,0,1,2,3,4
forward2 = (1/dx).*[-1/6,-77/60,5/2,-5/3,5/6,-1/4,1/30]; %-1,0,1,2,3,4,5
forward3 = (1/dx).*[-49/20,6,-15/2,20/3,-15/4,6/5,-1/6];%0,1,2,3,4,5,6




%a_nm_j: n is number of points to left, m = number of poitns to right, j is
%the node location
    %optimized coefficients from paper, thank god I didnt have to solve it
a_06_0 = -2.192280339;
a_60_0 = -1*a_06_0;
a_06_1 = 4.748611401; 
a_60_n1 = -1*a_06_1;
a_06_2 = -5.108851915; 
a_60_n2 = -1*a_06_2;
a_06_3 = 4.461567104; 
a_60_n3 = -1*a_06_3; 
a_06_4 = -2.833498741; 
a_60_n4 = -1*a_06_4; 
a_06_5 = 1.128328861; 
a_60_n5 = -1*a_06_5; 
a_06_6 = -0.203876371; 
a_60_n6 = -1*a_06_6; 

a_15_n1 = -0.209337622;
a_51_1 = -1*a_15_n1; 
a_15_0 = -1.084875676; 
a_51_0 = -1*a_15_0; 
a_15_1 = 2.147776050; 
a_51_n1 = -1*a_15_1; 
a_15_2 = -1.388928322; 
a_51_n2 = -1*a_15_2; 
a_15_3 = 0.768949766; 
a_51_n3 = -1*a_15_3; 
a_15_4 = -0.281814650; 
a_51_n4 = -1*a_15_4; 
a_15_5 = 0.048230454; 
a_51_n5 = -1*a_15_5; 

a_24_n2 = 0.049041958; 
a_42_2 = -1*a_24_n2; 
a_24_n1 = -0.468840357; 
a_42_1 = -1*a_24_n1; 
a_24_0 = -0.474760914; 
a_42_0 = -1*a_24_0;
a_24_1 = 1.273274737;
a_42_n1 = -1*a_24_1;
a_24_2 = -0.518484526; 
a_42_n2 = -1*a_24_2;
a_24_3 = 0.166138533; 
a_42_n3 = -1*a_24_3; 
a_24_4 = -0.026369431; 
a_42_n4 = -1*a_24_4; 


%matrix of optimized coefficients
optimized7Centered = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
optimized7Centered2 = [-0.208431427703,0.166705904415,-0.770882380518,0,0.770882380518,-0.166705904415,0.208431427703];

back1_opt = (1/dx).*[a_42_n4,a_42_n3,a_42_n2,a_42_n1,a_42_0,a_42_1,a_42_2]; % -4,-3,-2,-1,0,1,2
back2_opt = (1/dx).*[a_51_n5, a_51_n4, a_51_n3, a_51_n2, a_51_n1, a_51_0, a_51_1]; % -5,-4,-3,-2,-1,0,1
back3_opt = (1/dx).*[a_60_n6,a_60_n5,a_60_n4,a_60_n3,a_60_n2,a_60_n1,a_60_0];% -6,-5,-4,-3,-2,-1,0
forward1_opt = (1/dx).*[a_24_n2,a_24_n1,a_24_0,a_24_1,a_24_2,a_24_3,a_24_4]; %-2,-1,0,1,2,3,4
forward2_opt = (1/dx).*[a_15_n1,a_15_0,a_15_1,a_15_2,a_15_3,a_15_4,a_15_5]; %-1,0,1,2,3,4,5
forward3_opt = (1/dx).*[a_06_0,a_06_1,a_06_2,a_06_3,a_06_4,a_06_5,a_06_6];%0,1,2,3,4,5,6

clearvars a_06_0 a_06_1 a_06_2 a_06_3 a_06_4 a_06_5 a_06_6 a_15_0 a_15_1 a_15_2 a_15_3 a_15_4 a_15_5 a_15_n1 a_24_0 ...
    a_24_1 a_24_2 a_24_3 a_24_4 a_24_n1 a_24_n2 a_42_1 a_42_2 a_42_n1 a_42_n2 a_42_n3 a_42_n4 a_51_0 a_51_1 a_51_n1 ...
    a_51_n2 a_51_n3 a_51_n4 a_51_n5 a_60_0 a_60_n1 a_60_n2 a_60_n3 a_60_n4 a_60_n5 a_60_n6

%coeffM = [optimized7Centered; back1;back2;back3;forward1;forward2;forward3]; % for standard scheme, not useful, so lets use the optimized
coeffM = [optimized7Centered; back1_opt;back2_opt;back3_opt;forward1_opt;forward2_opt;forward3_opt]; % for standard scheme, not useful, so lets use the optimized


%% Running the code


%just do like one time step and lets see what we get. 
for n = 1:17
    Un = RK4(Un,n*dt,S,dt,mean_Values, gamma,c,1, xRange,yRange,coeffM);
end 




%% plotting things


close all
%plot density
figure1 = figure();
surf(X(:,:),Y(:,:),Un(:,:,1))
title("Density")
xlabel("xlabel")
ylabel("ylabel")

figure2 = figure();
surf(X(:,:),Y(:,:),Un(:,:,2))
title("U")
xlabel("x")
ylabel("y")

figure3 = figure();
surf(X(:,:),Y(:,:),Un(:,:,3))
title("V")
xlabel("x")
ylabel("y")

figure4 = figure();
surf(X(:,:),Y(:,:),Un(:,:,4))
title("pressure")
xlabel("x")
ylabel("y")

figure5 = figure();
plot(X(24,:),Un(24,:,4));
title("Pressure along y=0")
xlabel("x")
ylabel("y")

figure6 = figure()
contour(X(:,:),Y(:,:),abs(Un(:,:,4)));
title("Pressure contours")



%% Testing Section

    
    % testing Source, looks reasonable
    %{
    xRange = -23:dx:23; 
    yRange = -23:dy:23;
    [X,Y]= meshgrid(xRange,yRange);
    S_Test = STest(X,Y,pi/2);
    surf(X,Y,S_Test)
    %}


    %Testing E,F 7 point schemes, a test case is implemented in excel to
    %double check passed for E. LOOKS REASONABLE
    %{
    E_row= [1,2,3,4,5,6,7]; %use this foreverything
    F_Col = [7;6;5;4;3;2;1];
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    mean_Values = [1.225, 0, 0, 1E5];
    
    Test_E = assembleE7Point(a_Values,gamma,c,mean_Values,E_row,E_row,E_row,E_row);
    Test_E(1)
    Test_E(2)
    Test_E(3)
    Test_E(4)
    
    Test_F = assembleE7Point(a_Values,gamma,c,mean_Values,F_Col,F_Col,F_Col,F_Col);
    Test_F(1)
    Test_F(2)
    Test_F(3)
    Test_F(4)
    %}

    % testing spatial function to see if it runs.
    %{
     
    xRange = -4:4; 
    yRange = -4:4;
    
    % initialize empty data place holders
    

    [X,Y]= meshgrid(xRange,yRange);
    totalElement = size(X,1)*size(X,2);
    xElemNum = length(xRange);
    yElemNum = length(yRange);
    
    rng(1,'philox')
    U(:,:,1) = rand(xElemNum, yElemNum);
    rng(2,'philox')
    U(:,:,2) = rand(xElemNum, yElemNum);
    rng(3,'philox')
    U(:,:,3) = rand(xElemNum, yElemNum);
    rng(4,'philox')
    U(:,:,4) = rand(xElemNum, yElemNum);
    
    R = sqrt(X.^2+Y.^2);
    cosTheta = Y./R;
    sinTheta = X./R;

    
    mean_Values = [1.0, 0, 0, 1];
    K_Test = assembleK(0,S,mean_Values,gamma,c,U,xRange,yRange,coeffM);7
    
    %expected center value
    centerIndex = (xElemNum-1)/2+1;
    E_center = -1*dot(coeffM(1,:),U(centerIndex, centerIndex-3:centerIndex+3,2));
    F_center = -1*dot(coeffM(1,:),U(centerIndex-3:centerIndex+3,centerIndex,3));
    K_center = E_center+F_center;
    K_center == K_Test(centerIndex,centerIndex,1)
    
    %expected node valueK
    nodeI = 7;
    nodeJ = 3; 
    r= R(nodeI,nodeJ);
    cosT = cosTheta(nodeI,nodeJ);
    sinT = sinTheta(nodeI,nodeJ);
    E_node =1*dot(coeffM(5,:),U(nodeI, nodeJ-2:nodeJ+4,1))*cosT*c; 
    F_node = 1*dot(coeffM(2,:),U(nodeI-4:nodeI+2, nodeJ,1))*sinT*c; 
    U_node = c/(2*r)*U(nodeI,nodeJ,1);
    K_node =-1*(E_node+F_node+U_node);
    K_node == K_Test(nodeI,nodeJ)
    
    %expected bottom left corner 
    BL_I = 9;
    BL_J = 1; 
    r_BL= R(BL_I,BL_J);
    cosT_BL = cosTheta(BL_I,BL_J);
    sinT_BL = sinTheta(BL_I,BL_J);
    
    E_BL =1*dot(coeffM(7,:),U(BL_I, BL_J:BL_J+6,1))*cosT_BL*c; 
    F_BL = 1*dot(coeffM(4,:),U(BL_I-6:BL_I, BL_J,1))*sinT_BL*c; 
    U_BL = c/(2*r_BL)*U(BL_I,BL_J,1);
    K_BL =-1*(E_BL+F_BL+U_BL);
    K_BL == K_Test( BL_I, BL_J,1)
    
    

    %expected top right corner
    TR_I = 1;
    TR_J = 9;
    r_TR= R(TR_I,TR_J);
    cosT_TR = cosTheta(TR_I,TR_J);
    sinT_TR = sinTheta(TR_I,TR_J);
    
    E_TR =1*dot(coeffM(4,:),U(TR_I, TR_J-6:TR_J,1))*cosT_TR*c; 
    F_TR = 1*dot(coeffM(7,:),U(TR_I:TR_I+6, TR_J,1))*sinT_TR*c; 
    U_TR = c/(2*r_TR)*U(TR_I,TR_J,1);
    K_TR =-1*(E_TR+F_TR+U_TR);
    K_TR == K_Test( TR_I, TR_J,1)
    
%}
    
    %testing the damping term implementation:
    %{
    vector = 1:49; 
    matrix_2d = reshape(vector, [7,7]);
    matrix_3d(:,:,1) = matrix_2d;
    matrix_3d(:,:,2) = matrix_2d;
    matrix_3d(:,:,3) = matrix_2d;
    matrix_3d(:,:,4) = matrix_2d; 
    damped = dampingTerm(5,matrix_3d);
    
    d_Values = [0.023853048191,-0.106303578770,0.226146951809,0.287392842460,-0.226146951809,0.106303578770,-0.023853048191];
    row = [4,11,18,25,32,39,46];
    col = [22,23,24,25,26,27,28];
    expectedCenter = (-1/5).*dot(d_Values, row+col);
    
%}
    

    % Testing for Writing file options for debugging purposes.
    %{
    xRange = -4:4; 
    yRange = -4:4;
    
    % initialize empty data place holders
    [X,Y]= meshgrid(xRange,yRange);
    totalElement = size(X,1)*size(X,2);
    xElemNum = length(xRange);
    yElemNum = length(yRange);
    
    rng(1,'philox')
    U(:,:,1) = rand(xElemNum, yElemNum);
    rng(2,'philox')
    U(:,:,2) = rand(xElemNum, yElemNum);
    rng(3,'philox')
    U(:,:,3) = rand(xElemNum, yElemNum);
    rng(4,'philox')
    U(:,:,4) = rand(xElemNum, yElemNum);
    
    writematrix(U,'U.csv','Delimiter',',') %this will save it as 2D matrix, piled together.
 %}

    % Check for spatial du/dt computation again 
    %what does not make sense to me is why we are getting a skewed towards
    %the corners in the U, so lets just plot everything
    % turns out to be a mistake with the case in the interioirregion, i am
    % computing the interior equation insstead of the proper boundary when
    % one of the coordinate is in the interior region.
     
    %{
    xRange = -23:dx:23; 
    yRange = -23:dy:23;
    [X,Y]= meshgrid(xRange,yRange);
    close all
    
    figure1 = figure();
    U_17 = readmatrix('U_17.csv');
    U_17 = reshape(U_17, [47,47,4]);
    surf(X(:,:),Y(:,:),U_17(:,:,1))
    title("density at time step 17")
    xlabel("x")
    ylabel("y")
    
    figure2 = figure();
    U_18 = readmatrix('U_18.csv');
    U_18 = reshape(U_18, [47,47,4]);
    surf(X(:,:),Y(:,:),U_18(:,:,1))
    title("density at time step 18")
    xlabel("x")
    ylabel("y")
    
    figure3 = figure();
    K1_17 = readmatrix('K1_17.csv');
    K1_17 = reshape(K1_17, [47,47,4]);
    surf(X(:,:),Y(:,:),K1_17(:,:,1))
    title("K1 density at time step 17")
    xlabel("x")
    ylabel("y")
    
    % let us check a couple points manually for the update. 
    % if this is correctly implemented, then it is a question of the scheme
   
  
    R = sqrt(X.^2+Y.^2);
    cosTheta = Y./R;
    sinTheta = X./R;
    Source = STest(X,Y,17);
    
    %{ check more points tomorrow.
    % check x = 0, y = -23 
    j_Index = 24; 
    i_Index = 1;
    E_0_N23 = dot(coeffM(1,:),U_17(i_Index,j_Index-3:j_Index+3,1)).*(c*cosTheta(i_Index,j_Index));
    F_0_N23 = dot(coeffM(7,:),U_17(i_Index:i_Index+6,j_Index,1)).*(c*sinTheta(i_Index,j_Index));
    U_0_N23 = U_17(i_Index,j_Index).*(c/(2*R(i_Index,j_Index)));
    K1_0_N23 = -1.*(E_0_N23+F_0_N23+U_0_N23);
   
    K1_17_Computed = assembleK(17,S,mean_Values,gamma,c,U_17,xRange,yRange,coeffM);
    K1_0_N23 == K1_17_Computed(i_Index,j_Index)
    %}
%% Functions



% RK4: time marching scheme for 4 stage rung-kutta
% Inputs: 
%   Un: the solution of Un (x by y by 4), require non-dimesionalization. 
%   t: current time
%   S: source function handle 
%   dt: the timestep 
%   mean_Values: the ambient values
%   gamma: specific heat ratio 
%   c: speed of sound
%   Rs: reynolds number for the mesh
%   xRange:
%   xRange:
%   coeffM: matrix of spatial coefficients
% Outputs:
%   Unp1: the next time step values, non-dimensionalized here. 

function Unp1 = RK4(Un,t,S,dt,mean_Values, gamma, c,Rs,xRange,yRange,coeffM)
    %Standard RK 4
    t_vec = [t, t+1/2*dt, t+1/2*dt,t+dt];
    K1 = assembleK(t_vec(1),S,mean_Values,gamma,c,Un,xRange,yRange,coeffM);
    K2 = assembleK(t_vec(2),S,mean_Values,gamma,c,Un+(0.5*dt).*K1,xRange,yRange,coeffM);
    K3 = assembleK(t_vec(3),S,mean_Values,gamma,c,Un+(0.5*dt).*K2,xRange,yRange,coeffM);
    K4 = assembleK(t_vec(4),S,mean_Values,gamma,c,Un+(1*dt).*K3,xRange,yRange,coeffM);
    Unp1 = Un+ (1/6).*(K1+2.*K2+2.*K3+K4);
    
    
    %lets save all the stuff I want.
    %{
    if (t == 17)
        writematrix(Un,'U_17.csv','Delimiter',',')
        writematrix(K1,'K1_17.csv','Delimiter',',')
        writematrix(K2,'K2_17.csv','Delimiter',',')
        writematrix(K3,'K3_17.csv','Delimiter',',')
        writematrix(K4,'K4_17.csv','Delimiter',',')
        writematrix(Unp1,'U_18.csv','Delimiter',',')
              
    elseif(t == 18)
        
    end
    %}
    
    % paper version
    %{
    c_List = [1,0.5,0.162997,0.0407574];
    
    %calculate beta values
    beta_4 = c_List(2); 
    beta_3 = c_List(3)/beta_4;
    beta_2 = c_List(4)/beta_4/beta_3;
    beta_1 = 0;
    beta_List = [beta_1, beta_2, beta_3, beta_4];
    %one concern is that the time adds up over 1. 
    
    %other values 
    dx = xRange(2)-xRange(1);
    K = assembleK(t,S,mean_Values,gamma,c,Un,xRange,yRange,coeffM);
    %need to check CHECK HERE
    for i = 1:4
       U_Inter = Un +beta_List(i).*K;
       K = dt.*assembleK(t,S,mean_Values,gamma,c,U_Inter,xRange,yRange,coeffM);
    end
    
    Damping = dampingTerm(Rs, Un,dx);
    Unp1 = Un+ K+dt.*beta_List(4).*Damping;
    Unp1 = Un+ K;
    %}
    
end


% assembleK:assemble spatial discretization for estimation of time gradient
% Input:
%   t: current REAL LIFE time for source computation. 
%   S: function handle for the source term
%   mean_Values: list of mean values following same order as U(perturbation
%   vector)
%   gamma: specific heat ratio
%   c: speed of sound
%   U: the matricies storing the values of acoustic perturbation
%              variables at the current time step
%   dx: cell spacing in x 
%   dy: cellspacing in y
%   xRange: x-domain
%   yRange: y-domain
%   coeffM: matrix of spatial coefficients
% Output:
%   K: matrix of [X,Y,4] for the dU/dt of the solution. 
function K = assembleK(t,S,mean_Values,gamma,c,U,xRange,yRange,coeffM)
    %list of coefficients a-3,a-2,a-1,a0,a1,a2,a3
    a_Values = coeffM(1,:);
    
    %other constants
    dx = abs(xRange(2)-xRange(1));
    dy = abs(yRange(2)-yRange(1));
    interiorEndXDir = length(xRange)-3;
    interiorEndYDir = length(yRange)-3;
    jEnd = length(xRange);
    iEnd = length(yRange);
    
    %PlaceHolder
    K = zeros(length(yRange),length(xRange),4);


    %Better Method
     for i = 1:iEnd
        for j = 1:jEnd
            %i = 1; 
            %j = 24;
            
            r = sqrt(xRange(j)^2 + yRange(i)^2);
            cosTheta = yRange(i)/r;
            sinTheta=  xRange(j)/r;
            %initialize the terms
            E = zeros(4,1);
            F = zeros(4,1);
            
         
            
            % Lets compute the E(x derivatives)
            if( j>= 4 && j<=interiorEndXDir)
                
                %need a couple more cases here to check if we are in
                %interior or outside
                if(i >=  4 && i<= interiorEndYDir)
                     %we are in the interior x domain   
                     E = assembleE7Point(coeffM(1,:),gamma,c,mean_Values,U(i,j-3:j+3,1),U(i,j-3:j+3,2),U(i,j-3:j+3,3),U(i,j-3:j+3,4));
                     E = (1/dx).*E;
                else
                     rowIndices = j-3:j+3;
                     for z = 1:4
                        test = dot(coeffM(1,:),U(i,rowIndices,z) ).*(c*cosTheta);
                        E(z,1) = (1/dx).*test;
                     end
                end
            elseif(j<=3)
                %we are in the left boundary x domain (spatially)
                %need forward differencing here
                %E = c*cos*theta*dU/dx
                
                rowIndices = 1:7;
                for z = 1:4
                   test = dot(coeffM(7-j+1,:),U(i,rowIndices,z) ).*(c*cosTheta);
                   E(z,1) = test;
                end
                %E = assembleE7Point(coeffM(7-j+1,:),gamma,c,mean_Values,U(i,rowIndices ,1),U(i,rowIndices,2),U(i,rowIndices,3),U(i,rowIndices,4));
                
            else  
                %we are in the right boundary x domain (spatially)
                %need backward differencing here
                coeffIndex = 4-(jEnd-j);
                rowIndices = jEnd-6: jEnd;
                
                %E = c*cos*theta*dU/dx
                for z = 1:4
                   E(z,1) = dot(coeffM(coeffIndex,:), U(i,rowIndices,z) ).*(c*cosTheta);
                end
                
                %E = assembleE7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(i,rowIndices,1),U(i,rowIndices,2),U(i,rowIndices,3),U(i,rowIndices,4));
        
            end 
            
            %lets compute the F(y Derivatives)
            if( i>= 4 && i<=interiorEndYDir)
                if( j>= 4 && j<=interiorEndXDir)
                    %interior y domain
                    F = assembleF7Point(coeffM(1,:),gamma,c,mean_Values,U(i-3:i+3,j,1),U(i-3:i+3,j,2),U(i-3:i+3,j,3),U(i-3:i+3,j,4));
                    F = (1/dy).*F;
                else
                    coeffIndex = 1;
                    colIndices = i-3:i+3;
                
                    %F = c0*sinTheta*dU/dy
                    for z = 1:4
                       F(z,1) = dot(coeffM(coeffIndex,:), U(colIndices,j,z) ).*(c*sinTheta);
                    end
                end 
            elseif(i<=3)
                % we are at the bottom boundary region (spatially)
                % need forward differencing 
                coeffIndex = 7-i+1;
                colIndices = 1:7;
                
                %F = c0*sinTheta*dU/dy
                for z = 1:4
                   F(z,1) = dot(coeffM(coeffIndex,:), U(colIndices,j,z) ).*(c*sinTheta);
                end
                %F = assembleF7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(colIndices,j,1),U(colIndices,j,2),U(colIndices,j,3),U(colIndices,j,4));
                
            else
                % we are at the top Boundary region(spatially)
                % need backward differencing. 
                coeffIndex = 4 - (iEnd - i);
                colIndices = iEnd-6:iEnd;
                
                %F = c0*sinTheta*dU/dy
                for z = 1:4
                   F(z,1) = dot(coeffM(coeffIndex,:), U(colIndices,j,z) ).*(c*sinTheta);
                end
                
                %F = assembleF7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(colIndices,j,1),U(colIndices,j,2),U(colIndices,j,3),U(colIndices,j,4));
            end
            
            %K matrix ASSEMbly
            if( j>= 4 && j<=interiorEndXDir && i>= 4 && i<=interiorEndYDir )
                %inside intereior, one formulation for du/dt
                S_ij = S(xRange(j),yRange(i),t); 
                K(i,j,1:4) = -1.*(E+F)+S_ij;
            else
                %boundary region another expression for du/dt
                %du/dt = -c0*cosTheta*du/dx - c0*sinTheta*du/dy - c/2r*U
                %=-(E+F+c/2r*U)
                U_ij = reshape(U(i,j,:),[4,1]);
                K(i,j,1:4)  = -1.*(E+F+(c/(2*r)).*U_ij);
            end 
                
        end 
     end
   
end 

%dampingTerm: Implements damping for our scheme. 
%Input: 
%   Rs: reynolds number of mesh, larger means less damping
%   U: The solution matrix U
%   dx:
%   dy:
%Outputs:
%   D: the damping terms. 

function D = dampingTerm(Rs, U,dx)
 
    CD_7 =(1/dx).*[-0.023853048191,0.106303578770,-0.226146951809,0.287392842460,-0.226146951809,0.106303578770,-0.023853048191];
    %CD_5 = [0.0625,-0.25,0.0625,0.375,-0.25,0.0625];
    %CD_3 = [-0.25,0.5,-0.25];
    d_Matrix(1,:,1) = CD_7;
    d_Matrix(1,:,2) = CD_7;
    d_Matrix(1,:,3) = CD_7;
    d_Matrix(1,:,4) = CD_7;

    %other constants
    x_IndexRange = 4: size(U,2)-3; %the actual value cells, not including ghost cells;
    y_IndexRange = 4:size(U,1)-3;
    
    %PlaceHolder
    D = zeros(size(U,2),size(U,1),4);
    % construct the E,F matrices through for loop for each i,j
    for i = y_IndexRange
        for j = x_IndexRange
            %row values
            U_Row = U(i,j-3:j+3,:);
            U_Col = U(i-3:i+3,j,:);
            U_Col = permute(U_Col, [2 1 3]);
            
            D(i,j,:) = dot(-1/Rs.*d_Matrix,U_Row+U_Col);
        end
        
    end
end 

%assembleE7Point: Assemble the 7 point stencil for Ki,j 
%Input: 
%   a_Values: coefficients for the scheme vector of size (1,7)
%   gamma: specific heat ratio
%   c: speed of sound
%   mean_Values: list of mean values size (1,4)
%   
%   rho_Row: 7 points in the row for rho size (1,7)
%   u_Row: 7 points in the row for u size (1,7)
%   v_Row: 7 points in the row for v size (1,7)
%   p_Row: 7 points in the row for p size (1,7)
%Output:
%   E: the 7 point scheme for E_ij dimensions size (4,1)
function E= assembleE7Point(a_Values,gamma,c,mean_Values,rho_Row,u_Row,v_Row,p_Row)
    %place holder and constants
    E = zeros(4,1);
    rho_Mean = mean_Values(1);
    u_Mean = mean_Values(2);
    v_Mean = mean_Values(3);
    p_Mean = mean_Values(4);
    
    u_Row = u_Row;
    v_Row = v_Row; 
    
    % put together E
    E(1) = dot(a_Values, u_Row);
    E(2) = dot(a_Values, p_Row);
    E(3) = 0;
    E(4) = dot(a_Values, (gamma.*(p_Mean/rho_Mean).*u_Row));

end 
    

%assembleF7Point: Assemble the 7 point stencil for Ki,j 
%Input: 
%   a_Values: coefficients for the scheme vector of size (1,7)
%   gamma: specific heat ratio
%   c: speed of sound
%   mean_Values: list of mean values size (1,4)
%   rho_Col: 7 points in the column for rho size (7,1)
%   u_Col: 7 points in the column for u size (7,1)
%   v_Col: 7 points in the column for v size (7,1)
%   p_Col: 7 points in the column for p size (7,1)
%Output:
% F: the 7 point scheme for F at i,j dimensions size (4,1)
function F = assembleF7Point(a_Values,gamma,c,mean_Values,rho_Col,u_Col,v_Col,p_Col)
    %place holder and constants
    F = zeros(4,1);
    a_Values = a_Values';
    rho_Mean = mean_Values(1);
    u_Mean = mean_Values(2);
    v_Mean = mean_Values(3);
    p_Mean = mean_Values(4);
    
    % put together F
    F(1) = dot(a_Values, v_Col);
    F(2) = 0;
    F(3) = dot(a_Values, p_Col);
    F(4) = dot(a_Values, gamma.*(p_Mean/rho_Mean).*v_Col);

end 


