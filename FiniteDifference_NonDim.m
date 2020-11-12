%% set up the domain including 3 ghost cells
clear
close all
dx = 1; 
dy = 1; 
c = 343; %speed of sound in air
gamma = 1.4; %specific heat ratio
dt = 0.5;

xRange = -23:dx:23; 
yRange = -23:dy:23;

% constants
u_Mean = 0;
v_Mean = 0; 
rho_Mean = 1.225; 
p_Mean = c^2*rho_Mean; 
mean_Values = [rho_Mean, u_Mean, v_Mean, p_Mean];

% initialize the matrix
[X,Y]= meshgrid(xRange,yRange);
Un(:,:,1) = ones(size(X,1),size(X,2)).*rho_Mean;
Un(:,:,2) = zeros(size(X,1),size(X,2));
Un(:,:,3) = zeros(size(X,1),size(X,2));
Un(:,:,4) = ones(size(X,1),size(X,2)).*p_Mean./(rho_Mean*c^2);


%source function
epi= 0.5; 
alpha = log(2)/5;
w = 2*pi/30;
S = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t).*[1/(mean_Values(1)/dx*c);0;0;1/(mean_Values(1)*c^3/dx)]; %monopole only  density perturbation
STest = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t); %monopole only for density perturbation



%% Running the code


%just do like one time step and lets see what we get. 
for n = 1:1
    Un = RK4(Un,n*dt,S,dt,mean_Values, gamma,c,10, xRange,yRange);
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

    %{
    % testing spatial function to see if it runs. 
    xRange = -3:3; 
    yRange = -3:3;

    % initialize empty data place holders
    [X,Y]= meshgrid(xRange,yRange);
    all = zeros(size(X,1),size(X,2));
    all(4,:) = [1,2,3,4,5,6,7];
    all(:,4) = [1,2,3,4,5,6,7];
    K_test = assembleK(0,S,meanValues,gamma,all,all,all,all,xRange,yRange)
    %}

    %{
    %Testing the boundary formulation. OLD, du/dt = 0 and wrong 
    list = 1:121;
    list_Matrix = reshape(list, [11,11]);
    K_Test(:,:,1) = list_Matrix;
    K_Test(:,:,2) = list_Matrix;
    K_Test(:,:,3) = list_Matrix;
    K_Boundary = NoMeanFlowBoundary(K_Test);
    %}

%{
    %Testing the boundary formulation. OLD, du/dt = 0 and wrong 
    %looks good at first glance, might have to come back and do a proper
    %debug. 
    list = 1:121;
    list_Matrix = reshape(list, [11,11]);
    xRange = ones(1,11);
    yRange = xRange;
    K_Test(:,:,1) = list_Matrix;
    K_Test(:,:,2) = list_Matrix;
    K_Test(:,:,3) = list_Matrix;
    K_Test(:,:,4) = list_Matrix;
    K_Boundary = NoMeanFlowBoundary(K_Test,xRange, yRange,[1,1,1,1,1,1,1],1);

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
% Outputs:
%   Unp1: the next time step values, non-dimensionalized here. 

function Unp1 = RK4(Un,t,S,dt,mean_Values, gamma, c,Rs,xRange,yRange)
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
    K = assembleK(t,S,mean_Values,gamma,c,Un,xRange,yRange);
    
    %need to check CHECK HERE
    for i = 2:4
       
       U_Inter = Un +beta_List(i).*K;
       %boundary conditions patched here.
       U_Inter = NoMeanFlowBoundary(U_Inter,xRange,yRange,a_Values,c);
       K = dt.*assembleK(t,S,mean_Values,gamma,c,U_Inter,xRange,yRange);
       %K = U_Inter;
    end
   
    Damping = dampingTerm(Rs, Un);
    Unp1 = Un+ K+dt.*beta_List(4).*Damping;
    Unp1 = NoMeanFlowBoundary(Unp1,xRange,yRange,a_Values,c);
    
    
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
% Output:
%   K: matrix of [X,Y,4] for the dU/dt of the solution. 
function K = assembleK(t,S,mean_Values,gamma,c,U,xRange,yRange)
    %list of coefficients a-3,a-2,a-1,a0,a1,a2,a3
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
   
    rho = U(:,:,1);
    u = U(:,:,2);
    v = U(:,:,3);
    p = U(:,:,4);
    
    %other constants
    x_IndexRange = 4: size(rho,2)-3; %the actual value cells, not including ghost cells;
    y_IndexRange = 4:size(rho,1)-3;
    dx = abs(xRange(2)-xRange(1));
    dy = abs(yRange(2)-yRange(1));
    
    %PlaceHolder
    K = zeros(size(rho,2),size(rho,1),4);
    % construct the E,F matrices through for loop for each i,j
    for i = x_IndexRange
        for j = y_IndexRange
            
            %7 point things
            rho_Row = rho(i,j-3:j+3);
            rho_Col = u(i-3:i+3,j);
            u_Row = u(i,j-3:j+3);
            u_Col = u(i-3:i+3,j);
            v_Row = v(i,j-3:j+3);
            v_Col = v(i-3:i+3,j);
            p_Row = p(i,j-3:j+3);
            p_Col = p(i-3:i+3,j);
            
            
            % Construct E,F, call it Z
            E = assembleE7Point(a_Values,gamma,c,mean_Values,rho_Row,u_Row,v_Row,p_Row);
            F = assembleF7Point(a_Values,gamma,c,mean_Values,rho_Col,u_Col,v_Col,p_Col);
            Z_ij = E+F;
            
            % Source 
            S_ij = S(xRange(i),yRange(j),t);
            
            %K matrix ASSEMbly
            K(i,j,1:4) = (dx/c).*(-1.*Z_ij+S_ij);
        end 
    end
    
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
    E(4) = dot(a_Values, (gamma.*(p_Mean/c^2/rho_Mean).*u_Row));

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
    
    u_Col = u_Col;
    v_Col = v_Col; 
    
    % put together F
    F(1) = dot(a_Values, v_Col);
    F(2) = 0;
    F(3) = dot(a_Values, p_Col);
    F(4) = dot(a_Values, gamma.*(p_Mean/(c^2*rho_Mean)).*v_Col);

end 

%NoMeanFlowBoundary: populates ghost cell value for no mean flow value
% Inputs:
%   U: the cell values at the current time step 
%   xRange: the range of X-values
%   yRange: the range of y-values
%   a_Values: 
%   c: speed of sound
% Output:
%   U_boundary: boundary have been modified
function U_Boundary = NoMeanFlowBoundary(U,xRange,yRange,a_Values,c)
%function K_Boundary = NoMeanFlowBoundary(K)
 
    % Re-implementation
    U_Boundary = U; 
    r_Values =sqrt(max(xRange(4:end-3)).^2+ yRange(4:end-3).*yRange(4:end-3)); %assume square domain and origin source for easy computation,ie r is the same   
    coefficients = ((c)/a_Values(7))./(2.*r_Values);
  
    for i = 1:4
        U_4Col = U(4:end-3,4,i);
        U_endM3Col = U(4:end-3,end-3,i);
        U_4Row = U(4,4:end-3,i);
        U_EndM3Row = U(end-3,4:end-3,i);
        
        U_Boundary(4:end-3,1,i) =coefficients'.*U_4Col+ U_Boundary(4:end-3,7,i);
        U_Boundary(4:end-3,end,i) = coefficients'.*U_endM3Col+U_Boundary(4:end-3, end-6,i);
        U_Boundary(1,4:end-3,i) = coefficients.*U_4Row+U_Boundary(7,4:end-3,i);
        U_Boundary(end,4:end-3,i) = coefficients.* U_EndM3Row+U_Boundary(end-6,4:end-3,i);
    end 
    % kept the same. 
    U_Boundary(4:end-3,3,:) = U_Boundary(4:end-3,5,:);
    U_Boundary(4:end-3,2,:) = U_Boundary(4:end-3,6,:);
    U_Boundary(4:end-3,end-2,:) = U_Boundary(4:end-3,end-4,:);
    U_Boundary(4:end-3,end-1,:) = U_Boundary(4:end-3,end-5,:);
    
    U_Boundary(3,4:end-3,:) = U_Boundary(5,4:end-3,:);
    U_Boundary(2,4:end-3,:) = U_Boundary(6,4:end-3,:);
    U_Boundary(end-2,4:end-3,:) = U_Boundary(end-4, 4:end-3,:);
    U_Boundary(end-1,4:end-3,:) = U_Boundary(end-5, 4:end-3,:);
    
    
    
    %{
    
    % This part is incorrect, so I will comment it out for now and re-do
    % this really quick
    %just set the values to be equal to eachother because of symmetric but
    %opposite sign coefficients for spatial scheme.
    K_Boundary = K; 
    
    K_Boundary(4:end-3,1,:) = K_Boundary(4:end-3,7,:);
    K_Boundary(4:end-3,2,:) = K_Boundary(4:end-3,6,:);
    K_Boundary(4:end-3,3,:) = K_Boundary(4:end-3,5,:);
    K_Boundary(4:end-3,end,:) = K_Boundary(4:end-3, end-6,:);
    K_Boundary(4:end-3,end-1,:) = K_Boundary(4:end-3,end-5,:);
    K_Boundary(4:end-3,end-2,:) = K_Boundary(4:end-3,end-4,:);
    
    
    K_Boundary(1,4:end-3,:) = K_Boundary(7,4:end-3,:);
    K_Boundary(2,4:end-3,:) = K_Boundary(6, 4:end-3,:);
    K_Boundary(3,4:end-3,:) = K_Boundary(5, 4:end-3,:);
    K_Boundary(end,4:end-3,:) = K_Boundary(end-6,4:end-3,:);
    K_Boundary(end-1,4:end-3,:) = K_Boundary(end-5, 4:end-3,:);
    K_Boundary(end-2,4:end-3,:) = K_Boundary(end-4, 4:end-3,:);
    
    %}
end


