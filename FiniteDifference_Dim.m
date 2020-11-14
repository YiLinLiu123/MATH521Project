%% set up the domain including 3 ghost cells
clear
close all
dx = 1; 
dy = 1; 
c = 10; %speed of sound in air
gamma = 1.4; %specific heat ratio
dt = 1;

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

%matrix of coefficients
optimized7Centered = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    %website: https://bellaard.com/tools/Finite%20difference%20coefficient%20calculator/
back1 = (1/dx^4)*[-1/6,1,-3/2,-2/3,7/2,-3,5/6]; % -4,-3,-2,-1,0,1,2
back2 = (1/dx^4)*[5/6,-6,37/2,-92/3,57/2,-14,17/6]; % -5,-4,-3,-2,-1,0,1
back3 = (1/dx^4)*[17/6,-19,107/2,-242/3,137/2,-31,35/6];% -6,-5,-4,-3,-2,-1,0
forward1 = (1/dx^4)*[5/6,-3,7/2,-2/3,-3/2,1,-1/6]; %-2,-1,0,1,2,3,4
forward2 = (1/dx^4)*[17/6,-14,57/2,-92/3,37/2,-6,5/6]; %-1,0,1,2,3,4,5
forward3 = (1/dx^4)*[35/6,-31,137/2,-242/3,107/2,-19,17/6];%0,1,2,3,4,5,6

coeffM = [optimized7Centered; back1;back2;back3;forward1;forward2;forward3];


%% Running the code

%{
%just do like one time step and lets see what we get. 
for n = 1:5
    Un = RK4(Un,n*dt,S,dt,mean_Values, gamma,c,5, xRange,yRange,coeffM);
end 
%}


%% plotting things

%{
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
surf(X(4:end-3,4:end-3),Y(4:end-3,4:end-3),Un(4:end-3,4:end-3,4))
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
%}
%% Testing Section

    %{
    % testing Source, looks reasonable
    S_Test = STest(X,Y,0);
    surf(X,Y,S_Test)
    %}

    %{
    %Testing E,F 7 point schemes, a test case is implemented in excel to
    %double check passed for E. LOOKS REASONABLE
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
    xRange = -3:3; 
    yRange = -3:3;

    % initialize empty data place holders
    [X,Y]= meshgrid(xRange,yRange);
    vector = 1:49; 
    matrix = reshape(vector, [7,7]);
    U(:,:,1) = matrix; 
    U(:,:,2) = matrix./5;
    U(:,:,3) = matrix./3;
    U(:,:,4) = matrix./2;
    
    mean_Values = [1.0, 0, 0, 1];
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    K_Test = assembleK(0,S,mean_Values,gamma,c,U,xRange,yRange,coeffM);
    
    
    %U(1) Values
    E_Vec_C = U(4,:,2);
    F_Vec_C = U(:,4,3)';
    E = dot(coeffM(1,:),E_Vec_C);
    F = dot(coeffM(1,:),F_Vec_C');
    K_1_C = -E-F;
    K_1_1 = -1.*(dot(coeffM(1,:),E_Vec_C) + dot(coeffM(7,:),F_Vec_C));
    K_1_2 = -1.*(dot(coeffM(1,:),E_Vec_C) + dot(coeffM(6,:),F_Vec_C));
    K_1_3 = -1.*(dot(coeffM(1,:),E_Vec_C) + dot(coeffM(5,:),F_Vec_C));
    
    %center node U(2)
    E_Vec_C2 = U(4,:,4);
    K_2_C = -1.*dot(coeffM(1,:),E_Vec_C2); %the rest should be the same since the stencils are the same and the same weights
    
    %U(3)
    F_Vec_C3 =U(:,4,4);
    K_3_C = -1.*dot(coeffM(1,:),F_Vec_C3);
    K_3_1 = -1.* dot(coeffM(7,:),F_Vec_C3);
    K_3_2 = -1.* dot(coeffM(6,:),F_Vec_C3);
    K_3_3 = -1.* dot(coeffM(5,:),F_Vec_C3);
    
    %U(4) values are just scaled by gamma, easy to check
%{
    %Testing the boundary formulation. OLD, du/dt = 0 and wrong 
    %looks good at first glance, might have to come back and do a proper
    %debug. 
    list = 1:121;
    matrix = reshape(list, [11,11]);
    U(:,:,1) = matrix;
    U(:,:,2) = matrix./5;
    U(:,:,3) = matrix./3;
    U(:,:,4) = matrix./2;
    
    xRange = -5:5; 
    yRange = -5:5;
    c=1;
    gamma = 1.4;
    
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    K =  assembleK(0,S,mean_Values,gamma,c,U,xRange,yRange);
    U_Boundary = NoMeanFlowBoundary(c,U,K,xRange,yRange);
    K_Boundary = assembleK(0,S,mean_Values,gamma,c,U_Boundary,xRange,yRange);
    topLeftU = 4+-12/a_Values(1) + c/(2*sqrt(8)*a_Values(1))*37;
    %}

    %{
    %testing the damping term implementation:
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
    c_List = [1,0.5,0.162997,0.0407574];
    
    %calculate beta values
    beta_4 = c_List(2); 
    beta_3 = c_List(3)/beta_4;
    beta_2 = c_List(4)/beta_4/beta_3;
    beta_1 = 0;
    beta_List = [beta_1, beta_2, beta_3, beta_4];
    
    %list of coefficients a-3,a-2,a-1,a0,a1,a2,a3
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    
    %calculate the stages
    K = assembleK(t,S,mean_Values,gamma,c,Un,xRange,yRange,coeffM);

    
    %need to check CHECK HERE
    for i = 2:4
       
       U_Inter = Un +beta_List(i).*K;
       %boundary conditions patched here.
       U_Inter = NoMeanFlowBoundary(c,U_Inter,K,xRange,yRange);
       K = dt.*assembleK(t,S,mean_Values,gamma,c,U_Inter,xRange,yRange,coeffM);
       %K = U_Inter;
    end
   
    Damping = dampingTerm(Rs, Un);
    Unp1 = Un+ K+dt.*beta_List(4).*Damping;
    Unp1 = NoMeanFlowBoundary(c,Unp1,xRange,yRange,coeffM);
    
    
end


%NoMeanFlowBoundary: populates ghost cell value for no mean flow value
% Input:
%   c: speed of sound
%   U: the matricies storing the values of acoustic perturbation
%              variables at the current time step
%   dx: cell spacing in x 
%   dy: cellspacing in y
%   xRange: x-domain
%   yRange: y-domain
%   coeffM: coefficient for spatial discretizations
% Output:
%   U_boundary: boundary have been modified
function U_Boundary = NoMeanFlowBoundary(c,U,K,xRange,yRange,coeffM)
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    U_Boundary = U;
    U_Boundary;

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
    x_IndexRange = 4: length(xRange)-3; %the actual value cells, not including ghost cells;
    y_IndexRange = 4:length(yRange)-3;
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
            %initialize the terms
            E = zeros(4,1);
            F = zeros(4,1);
            S_ij = S(xRange(j),yRange(i),t); 
            
            % Lets compute the E(x derivatives)
            if( j>= 4 && j<=interiorEndXDir)
                %we are in the interior x domain
                E = assembleE7Point(coeffM(1,:),gamma,c,mean_Values,U(i,j-3:j+3,1),U(i,j-3:j+3,2),U(i,j-3:j+3,3),U(i,j-3:j+3,4));
                E = (1/dx).*E;
            elseif(j<=3)
                %we are in the left boundary x domain (spatially)
                %need forward differencing here
                rowIndices = 1:7;
                E = assembleE7Point(coeffM(7-j+1,:),gamma,c,mean_Values,U(i,rowIndices ,1),U(i,rowIndices,2),U(i,rowIndices,3),U(i,rowIndices,4));
            else  
                %we are in the right boundary x domain (spatially)
                %need backward differencing here
                coeffIndex = 4-(jEnd-j);
                rowIndices = jEnd-6: jEnd;
                E = assembleE7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(i,rowIndices,1),U(i,rowIndices,2),U(i,rowIndices,3),U(i,rowIndices,4));
        
            end 
            
            %lets compute the F(y Derivatives)
            if( i>= 4 && i<=interiorEndYDir)
                %we are in the interior y domain
                F = assembleF7Point(coeffM(1,:),gamma,c,mean_Values,U(i-3:i+3,j,1),U(i-3:i+3,j,2),U(i-3:i+3,j,3),U(i-3:i+3,j,4));
                F = (1/dx).*F;
            elseif(i<=3)
                % we are at the bottom boundary region (spatially)
                % need forward differencing 
                coeffIndex = 7-i+1;
                colIndices = 1:7;
                F = assembleF7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(colIndices,j,1),U(colIndices,j,2),U(colIndices,j,3),U(colIndices,j,4));
                
            else
                % we are at the top Boundary region(spatially)
                coeffIndex = 4 - (iEnd - i);
                colIndices = iEnd-6:iEnd;
                F = assembleF7Point(coeffM(coeffIndex,:),gamma,c,mean_Values,U(colIndices,j,1),U(colIndices,j,2),U(colIndices,j,3),U(colIndices,j,4));
            end
            
            %K matrix ASSEMbly
            K(i,j,1:4) = -1.*(E+F)+S_ij;
        end 
     end
    
    
    %{
    % construct the E,F matrices through for loop for all of I,J, just do
    % cases
    for i = 1:length(yRange)
        for j = 1:length(xRange)
            E = zeros(4,1);
            F = zeros(4,1);
            if( i<= length(yRange)-3 && i>= 4 && j<= length(xRange)-3 && j>=4)
                rho_Row = U(i,j-3:j+3,1);
                rho_Col = U(i-3:i+3,j,1);
                u_Row = U(i,j-3:j+3,2);
                u_Col = U(i-3:i+3,j,2);
                v_Row = U(i,j-3:j+3,3);
                v_Col = U(i-3:i+3,j,3);
                p_Row = U(i,j-3:j+3,4);
                p_Col = U(i-3:i+3,j,4);
                % Construct E,F, call it Z
                E = assembleE7Point(a_Values,gamma,c,mean_Values,rho_Row,u_Row,v_Row,p_Row);
                F = assembleF7Point(a_Values,gamma,c,mean_Values,rho_Col,u_Col,v_Col,p_Col);
            
            %{    
            elseif( i >= 1 && i<=3 && j>= 4 && j<= length(xRange)-3)
                %bottom portion 
                
                %derivative in Y is all forward differencing
                F = assembleF7Point(coeffM(7-i+1,:),gamma,c,mean_Values,U(1:7,j,1),U(1:7,j,2),U(1:7,j,3),U(1:7,j,3));

                %derivative in X is the same everywhere because we have enough
                %nodes
                E = assembleE7Point(coeffM(1,:),gamma,c,mean_Values,U(i,j-3:j+3,1),U(i,j-3:j+3,2),U(i,j-3:j+3,3),U(i,j-3:j+3,4));

            end
            %}
                
            Z_ij = E.*(1/dx)+F.*(1/dy);
            % Source 
            S_ij = S(xRange(j),yRange(i),t); 
            
            %K matrix ASSEMbly
            K(i,j,1:4) = -1.*Z_ij+S_ij;
        end 
    end
    %}
    %{
    %Now do the bottom portion of the boundary excluding the 3x3
    for i = 1:3
        for j = 4:length(xRange)-3
            %derivative in Y is straight forward. 
            F = assembleF7Point(coeffM(7-i+1,:),gamma,c,mean_Values,U(1:7,j,1),U(1:7,j,2),U(1:7,j,3),U(1:7,j,3));
            
            %derivative in X is the same everywhere because we have enough
            %nodes
            E = assembleE7Point(coeffM(1,:),gamma,c,mean_Values,U(i,j-3:j+3,1),U(i,j-3:j+3,2),U(i,j-3:j+3,3),U(i,j-3:j+3,4));
            
            %source computation
            S_ij = S(xRange(j),yRange(i),t);
            
            K(i,j,1:4) = -1.*F+-1.*E+S_ij;
        end
    end 
    %}
    %}
end 

%dampingTerm: Implements damping for our scheme. 
%Input: 
%   Rs: reynolds number of mesh, larger means less damping
%   U: The solution matrix U
%Outputs:
%   D: the damping terms. 

function D = dampingTerm(Rs, U)

    d_Values = [0.023853048191,-0.106303578770,0.226146951809,0.287392842460,-0.226146951809,0.106303578770,-0.023853048191];
    d_Matrix(1,:,1) = d_Values;
    d_Matrix(1,:,2) = d_Values;
    d_Matrix(1,:,3) = d_Values;
    d_Matrix(1,:,4) = d_Values;
    
    
    %other constants
    x_IndexRange = 4: size(U,2)-3; %the actual value cells, not including ghost cells;
    y_IndexRange = 4:size(U,1)-3;
    
    %PlaceHolder
    D = zeros(size(U,2),size(U,1),4);
    % construct the E,F matrices through for loop for each i,j
    for i = x_IndexRange
        for j = y_IndexRange
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
