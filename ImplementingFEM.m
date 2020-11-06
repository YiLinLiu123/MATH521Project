%% Lets just start to try and generate a mesh object for our domain
%model = createpde(4);
%% litearlly just gives me a box
model = createpde(4);
%R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
R1 = [3,4,-200,200,200,-200,200,200,-200,-200]';
gm = R1;
sf = 'R1';

ns = char('R1');
ns = ns';
g = decsg(gm,sf,ns);

geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
axis equal
xlim([-250,250])
ylim([-250,250])

%% Now playing with meshing
%for mesh object documentation search up FEMesh object
%mesh data is also useful thing to search up
% General PDE: https://www.mathworks.com/help/pde/pde-problem-setup.html

%Ok so lets try specifyin some coefficients and stuff and lets see if i can
%get it working tomorrow real quick
generateMesh(model,"Hmax", 15)
pdeplot(model)