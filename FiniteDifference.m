%% set up the domain including 3 ghost cells
clear
dx = 1; 
dy = 1; 
c = 1; %speed of sound in air
gamma = 1.4; %specific heat ratio
dt = dx/c;

xRange = -5:dx:5; 
yRange = -5:dy:5;

% initialize empty data place holders
[X,Y]= meshgrid(xRange,yRange);
Un(:,:,1) = zeros(size(X,1),size(X,2));
Un(:,:,2) = zeros(size(X,1),size(X,2));
Un(:,:,3) = zeros(size(X,1),size(X,2));
Un(:,:,4) = zeros(size(X,1),size(X,2));

% constants
u_Mean = 0;
v_Mean = 0; 
rho_Mean = 1.225; 
p_Mean = 100E3; 
mean_Values = [rho_Mean, u_Mean, v_Mean, p_Mean];

%source function
epi= 0.5; 
alpha = log(2)/5;
w = 2*pi/30;
S = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t).*[1;0;0;1]; %monopole only  density perturbation
STest = @(x1,x2,t) epi.*exp(-alpha.*(x1.^2+x2.^2)).*sin(w.*t); %monopole only for density perturbation

%% Running the code

%just do like one time step and lets see what we get. 
for n = 1:2
    Un = RK4(Un,n*dt,S,dt,mean_Values, gamma, xRange,yRange);
    
end 

%plot density
surf(X,Y,Un(:,:,1))

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
    
    Test_E = assembleE7Point(a_Values,gamma,mean_Values,E_row,E_row,E_row,E_row);
    Test_E(1)
    Test_E(2)
    Test_E(3)
    Test_E(4)
    
    Test_F = assembleE7Point(a_Values,gamma,mean_Values,F_Col,F_Col,F_Col,F_Col);
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
    %Testing the boundary formulation. Looks Good
    list = 1:121;
    list_Matrix = reshape(list, [11,11]);
    K_Test(:,:,1) = list_Matrix;
    K_Test(:,:,2) = list_Matrix;
    K_Test(:,:,3) = list_Matrix;
    K_Boundary = NoMeanFlowBoundary(K_Test);
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
%   xRange:
%   xRange:
% Outputs:
%   Unp1: the next time step values, non-dimensionalized here. 

function Unp1 = RK4(Un,t,S,dt,mean_Values, gamma, xRange,yRange)
    c_List = [1,0.5,0.162997,0.0407574];
    
    %calculate beta values
    beta_4 = c_List(2); 
    beta_3 = c_List(3)/beta_4;
    beta_2 = c_List(4)/beta_4/beta_3;
    beta_1 = 0;
    beta_List = [beta_1, beta_2, beta_3, beta_4];
    
    %calculate the stages
    K = assembleK(t,S,mean_Values,gamma,Un(:,:,1),Un(:,:,2),Un(:,:,3),Un(:,:,4),xRange,yRange);
    for i = 2:4
       U_Inter = Un +beta_List(i).*K;
       K = dt.*assembleK(t,S,mean_Values,gamma,U_Inter(:,:,1),U_Inter(:,:,2),U_Inter(:,:,3),U_Inter(:,:,4),xRange,yRange);
    end
    Unp1 = Un+ K;
    
end


% assembleK:assemble spatial discretization for estimation of time gradient
% Input:
%   t: current REAL LIFE time for source computation. 
%   S: function handle for the source term
%   mean_Values: list of mean values following same order as U(perturbation
%   vector)
%   gamma: specific heat ratio
%   rho,u,v,p: the matricies storing the values of acoustic perturbation
%              variables at the current time step
%   dx: cell spacing in x 
%   dy: cellspacing in y
%   xRange: x-domain
%   yRange: y-domain
% Output:
%   K: matrix of [X,Y,4] for the dU/dt of the solution. 
function K = assembleK(t,S,mean_Values,gamma,rho,u,v,p,xRange,yRange)
    %list of coefficients a-3,a-2,a-1,a0,a1,a2,a3
    a_Values = [-0.02651995,0.18941314,-0.79926643,0,0.79926643,-0.18941314,0.02651995];
    
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
            E = assembleE7Point(a_Values,gamma,mean_Values,rho_Row,u_Row,v_Row,p_Row);
            F = assembleF7Point(a_Values,gamma,mean_Values,rho_Col,u_Col,v_Col,p_Col);
            Z_ij = (1/dx).*E+ (1/dy).*F;
            
            % Source 
            S_ij = S(xRange(i),yRange(j),t);
            
            %K matrix ASSEMbly
            K(i,j,1:4) = Z_ij+S_ij;
        end 
    end
    %boundary conditions
    K = NoMeanFlowBoundary(K);
end 


%assembleE7Point: Assemble the 7 point stencil for Ki,j 
%Input: 
%   a_Values: coefficients for the scheme vector of size (1,7)
%   gamma: specific heat ratio
%   mean_Values: list of mean values size (1,4)
%   rho_Row: 7 points in the row for rho size (1,7)
%   u_Row: 7 points in the row for u size (1,7)
%   v_Row: 7 points in the row for v size (1,7)
%   p_Row: 7 points in the row for p size (1,7)
%Output:
%   E: the 7 point scheme for E_ij dimensions size (4,1)
function E= assembleE7Point(a_Values,gamma,mean_Values,rho_Row,u_Row,v_Row,p_Row)
    %place holder and constants
    E = zeros(4,1);
    rho_Mean = mean_Values(1);
    u_Mean = mean_Values(2);
    v_Mean = mean_Values(3);
    p_Mean = mean_Values(4);
    
    u_Row = u_Row./rho_Mean; 
    v_Row = v_Row./rho_Mean; 
    
    % put together E
    E(1) = dot(a_Values, (u_Mean.*rho_Row + rho_Mean.*u_Row));
    E(2) = dot(a_Values, (u_Mean.*rho_Mean.*u_Row + p_Row));
    E(3) = dot(a_Values, (u_Mean.*rho_Mean.*v_Row));
    E(4) = dot(a_Values, (u_Mean.*p_Row + gamma.*p_Mean.*u_Row));

end 
    

%assembleF7Point: Assemble the 7 point stencil for Ki,j 
%Input: 
%   a_Values: coefficients for the scheme vector of size (1,7)
%   gamma: specific heat ratio
%   mean_Values: list of mean values size (1,4)
%   rho_Col: 7 points in the column for rho size (7,1)
%   u_Col: 7 points in the column for u size (7,1)
%   v_Col: 7 points in the column for v size (7,1)
%   p_Col: 7 points in the column for p size (7,1)
%Output:
% F: the 7 point scheme for F at i,j dimensions size (4,1)
function F = assembleF7Point(a_Values,gamma,mean_Values,rho_Col,u_Col,v_Col,p_Col)
    %place holder and constants
    F = zeros(4,1);
    a_Values = a_Values';
    rho_Mean = mean_Values(1);
    u_Mean = mean_Values(2);
    v_Mean = mean_Values(3);
    p_Mean = mean_Values(4);
    
    u_Col = u_Col./rho_Mean; 
    v_Col = v_Col./rho_Mean; 
    
    % put together F
    F(1) = dot(a_Values, (rho_Col.*v_Mean + rho_Mean.*v_Col));
    F(2) = dot(a_Values, (v_Mean.*rho_Mean.*u_Col));
    F(3) = dot(a_Values, (v_Mean.*rho_Mean.*v_Col+p_Col));
    F(4) = dot(a_Values, (v_Mean.*p_Col+gamma.*p_Mean.*v_Col));

end 

%NoMeanFlowBoundary: populates ghost cell value for no mean flow value%
% Inputs:
%   K: the spatial scheme output. 
% Output:
%   K_boundary: boundary have been modified
function K_Boundary = NoMeanFlowBoundary(K)

    %just set the values to be equal to eachother because of symmetric but
    %opposite sign coefficients for spatial scheme.
    K_Boundary = K; 
    
    K_Boundary(2:end-1,1,:) = K_Boundary(2:end-1,7,:);
    K_Boundary(2:end-1,2,:) = K_Boundary(2:end-1,6,:);
    K_Boundary(2:end-1,3,:) = K_Boundary(2:end-1,5,:);
    K_Boundary(2:end-1,end,:) = K_Boundary(2:end-1, end-6,:);
    K_Boundary(2:end-1,end-1,:) = K_Boundary(2:end-1,end-5,:);
    K_Boundary(2:end-1,end-2,:) = K_Boundary(2:end-1,end-4,:);
    
    
    K_Boundary(1,4:end-4,:) = K_Boundary(7,4:end-4,:);
    K_Boundary(2,4:end-4,:) = K_Boundary(6, 4:end-4,:);
    K_Boundary(3,4:end-4,:) = K_Boundary(5, 4:end-4,:);
    K_Boundary(end,4:end-4,:) = K_Boundary(end-6,4:end-4,:);
    K_Boundary(end-1,4:end-4,:) = K_Boundary(end-5, 4:end-4,:);
    K_Boundary(end-2,4:end-4,:) = K_Boundary(end-4, 4:end-4,:);
    
    
end


