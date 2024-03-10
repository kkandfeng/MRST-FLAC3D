%% This file only calculates the Liujiagou group, which is convenient for different hole pressures
clc; clear;
mrstModule add ad-core ad-props co2lab coarsegrid mrst-gui

%% Fluid parameters
gravity reset on;
g       = gravity;
rhow    = 1000;                                 % water density (kg/m^3)
co2     = CO2props();                           % CO2 property functions
% p_ref   = 30 *mega*Pascal;                      % reference pressure
% t_ref   = 94+273.15;                            % reference temperature

p_ref   = 28 *mega*Pascal;                      % reference pressure
t_ref   = 124+273.15;                            % reference temperature

co2_rho = co2.rho(p_ref, t_ref);                % CO2 density
co2_c   = co2.rhoDP(p_ref, t_ref) / co2_rho;    % CO2 compressibility
wat_c   = 0;                                    % water compressibility
c_rock  = 4.35e-5 / barsa;                      % rock compressibility
srw     = 0.27;                                 % residual water
src     = 0.20;                                 % residual CO2
pe      = 5 * kilo * Pascal;                    % capillary entry pressure
muw     = 8e-4 * Pascal * second;               % brine viscosity
muco2   = co2.mu(p_ref, t_ref) * Pascal * second; % co2 viscosity


%% Grid, petrophysical data, well, and initial state
% Grid and open boundary faces
[x_nnodes, y_nnodes] = deal(71, 71);
% x = linspace(0, 5000, x_nnodes);
% y = linspace(0, 5000, y_nnodes);
%load('D:\\temp\\mrst_ve_vtk\\ratio_coord.mat');
load('flac_gp_xy_out.mat');
x = new_x;
y = new_y;
%z = [0:100:400, 440:40:600, 640:40:800, 870:70:1080, 1120:40:1200, 1300:100:1500, 1550];
z = [0:17.5:280, 290:10:400];
depth = 3000.0; % 埋深
G = tensorGrid(x, y, z, 'depthz', repmat(depth, 1, x_nnodes * y_nnodes));
G = computeGeometry(G);



n_onelayer = 70*70;
poro = ones(G.cartDims)*0.15;
% poro(             1 : n_onelayer*4) = 0.015;
% poro(n_onelayer*4+1 : n_onelayer*9) = 0.15;
% poro(n_onelayer*9+1 : n_onelayer*14) = 0.014;
% poro(n_onelayer*14+1 : n_onelayer*18) = 0.014;
% poro(n_onelayer*18+1 : n_onelayer*21) = 0.14;
poro(              1 : n_onelayer*16) = 0.014;
poro(n_onelayer*4+1  : n_onelayer*28) = 0.14;
% poro(n_onelayer*21+1 : n_onelayer*24) = 0.010;
% poro(n_onelayer*24+1 : n_onelayer*25) = 0.1;
perm = ones(G.cartDims)*5.0;
% perm(             1 : n_onelayer*4) = 0.2;
% perm(n_onelayer*4+1 : n_onelayer*9) = 20;
% perm(n_onelayer*9+1 : n_onelayer*14) = 0.2;
% perm(n_onelayer*14+1 : n_onelayer*18) = 0.1;
% perm(n_onelayer*18+1 : n_onelayer*21) = 10;
perm(             1 : n_onelayer*16) = 0.05;
perm(n_onelayer*4+1 : n_onelayer*28) = 5;
% perm(n_onelayer*21+1 : n_onelayer*24) = 0.04;
% perm(n_onelayer*24+1 : n_onelayer*25) = 4;
ntg = ones(G.cartDims);
ntg(             1 : n_onelayer*12) = 0.0;
rock = makeRock(G, perm(:) .* milli * darcy, poro(:), 'ntg', ntg(:));

% % boundary 3D
% nx = G.cartDims(1); ny=G.cartDims(2); nz=G.cartDims(3);
% ix1 = searchForBoundaryFaces(G, 'BACK');
% ix2 = searchForBoundaryFaces(G, 'LEFT');
% ix3 = searchForBoundaryFaces(G, 'RIGHT');
% ix4 = searchForBoundaryFaces(G, 'FRONT');
% bcIx = [ix1; ix2; ix3; ix4];
%% %% Confirm the index of the well and locate the location of the well
precision = 10; % Custom precision (number of digits after the decimal point)

A = G.cells.centroids(:, 1);
B = G.cells.centroids(:, 2);
C = G.cells.centroids(:, 3);
rounded_A = round(A, precision);
A = unique(rounded_A);
rounded_B = round(B, precision);
B = unique(rounded_B);
rounded_C = round(C, precision);
C = unique(rounded_C);
% The coordinate point to be queried
x_query = 15529; % Example x-coordinates
y_query = 16107; % Example y-coordinates
z_query = 3300; % Example z-coordinates

% Use geometry information to find the index of the grid point closest to a given coordinate
[~, i] = min(abs(A - x_query));
[~, j] = min(abs(B - y_query));
[~, k] = min(abs(C - z_query));

% Outputs the index of the grid point closest to the coordinates found
fprintf('Closest coordinates (%f, %f, %f) The grid point index of is ：(%d, %d, %d)\n', x_query, y_query, z_query, i, j, k);
%% 

% Setup the well
% wc1 = 70 * 70 * 8 + 70*30+31;      
% inj_rate1 = 1.6e8 / year / co2_rho; % 纸坊组 1.6e5 t /a
W = [];
% Add a well to the set
% W = addWell(W, G, rock, wc, ...
%             'refDepth', G.cells.centroids(wc, 3), ... % BHP reference depth
%             'type', 'rate', ...  % inject at constant rate
%             'val', inj_rate, ... % volumetric injection rate
%             'comp_i', [0 1]);    % inject CO2, not water
% W = addWell(W, G, rock, wc1, ...            
%             'type', 'rate', ...  % inject at constant rate
%             'val', inj_rate1, ... % volumetric injection rate
%             'comp_i', [0 1]);    % inject CO2, not water         
        
wc2 = 70 * 70 * 24 + 70*30+31;      
inj_rate2 = 1.0e8 / year / co2_rho; 
W = addWell(W, G, rock, wc2, 'type', 'rate', ...  
            'val', inj_rate2, 'comp_i', [0 1]);  
        
% wc3 = 70 * 70 * 24 + 70*30+31;      
% inj_rate3 = 1.0e8 / year / co2_rho; 
% W = addWell(W, G, rock, wc3, 'type', 'rate', ...  
%             'val', inj_rate3, 'comp_i', [0 1]);            
        
% Top surface grid, petrophysical data, well, and initial state
[Gt, G, transMult] = topSurfaceGrid(G);
rock2D             = averageRock(rock, Gt);
W2D                = convertwellsVE(W, G, Gt, rock2D);
%initState.pressure = rhow * g(3) * Gt.cells.z;

initState.pressure = repmat(28e6 , Gt.cells.num, 1);
%initState.pressure(n_onelayer*21+1 : n_onelayer*25) = 32e6;

initState.s        = repmat([1, 0], Gt.cells.num, 1);
initState.sGmax    = initState.s(:,2);
% boundary 2D
nx = Gt.cartDims(1); ny=Gt.cartDims(2);
ix1 = searchForBoundaryFaces(Gt, 'BACK');
ix2 = searchForBoundaryFaces(Gt, 'LEFT');
ix3 = searchForBoundaryFaces(Gt, 'RIGHT');
ix4 = searchForBoundaryFaces(Gt, 'FRONT');
bcIxVE = [ix1; ix2; ix3; ix4];
%% 
% To avoid plotting artifacts when visualizing the volumetric and the
% top-surface grids, we shift the volumetric grid 100 meters down,
GG = G; 
GG.nodes.coords(:,3) = GG.nodes.coords(:,3) + 100;
screensize = get(0,'screensize'); 
figure('position',[screensize(3:4)-[845 565] 840 480]);
plotGrid(GG, 'facecolor', [1 1 .7]); 
plotGrid(Gt, 'facecolor', [.4 .5 1]);
[~,ht]=plotWell(G,W); set(ht,'FontSize',10, 'BackgroundColor',[.8 .8 .8]);
view(-65,33); clear GG;


%% Fluid model
% The PVT behavior of the injected CO2 is assumed to be given by an
% equation state, whereas the brine is incompressible. We also include a
% capillary fringe model based on upscaled, sampled capillary pressure.
% Please consult the documentation for the makeVEFluid routine for a
% description of parameters you can use to set up the various fluid models
% implemented in MRST-co2lab.
invPc3D = @(pc) (1-srw) .* (pe./max(pc, pe)).^2 + srw;
kr3D    = @(s) max((s-src)./(1-src), 0).^2; % uses CO2 saturation
fluid   = makeVEFluid(Gt, rock, 'P-scaled table'             , ...
               'co2_mu_ref'  , muco2, ...%6e-5 * Pascal * second , ...
               'wat_mu_ref'  , muw, ...%8e-4 * Pascal * second , ...
               'co2_rho_ref' , co2_rho                , ...
               'wat_rho_ref' , rhow                   , ...
               'co2_rho_pvt' , [co2_c, p_ref]         , ...
               'wat_rho_pvt' , [wat_c, p_ref]         , ...
               'residual'    , [srw, src]             , ...
               'pvMult_p_ref', p_ref                  , ...
               'pvMult_fac'  , c_rock                 , ...
               'invPc3D'     , invPc3D                , ...
               'kr3D'        , kr3D                   , ...
               'transMult'   , transMult);

%% Set up simulation schedule
% The simulation will consist of two periods: during the first 50 years,
% CO2 is injected as constant rate from the single injector. This period is
% simulated with a time step of 1 year. We then simulate a post-injection
% period of 950 years using time steps of 10 years.

% hydrostatic pressure conditions for open boundary faces
%p_bc     = Gt.faces.z(bcIxVE) * rhow * g(3);
p_bc = 28e6;
bc2D     = addBC([], bcIxVE, 'pressure', p_bc); 
bc2D.sat = repmat([1 0], numel(bcIxVE), 1);

% Setting up two copies of the well and boundary specifications. 
% Modifying the well in the second copy to have a zero flow rate.
schedule.control    = struct('W', W2D, 'bc', bc2D);
% schedule.control(2) = struct('W', W2D, 'bc', bc2D);
% schedule.control(2).W.val = 0;

% Specifying length of simulation timesteps
schedule.step.val = [repmat(year/12,    60, 1)];

% Specifying which control to use for each timestep.
% The first 100 timesteps will use control 1, the last 100
% timesteps will use control 2.
schedule.step.control = [ones(60, 1)];

%% Create and simulate model
model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);
[wellSol, states] = simulateScheduleAD(initState, model, schedule);
states = [{initState} states(:)'];

%% Animate the plume migration over the whole simulation period
clf
oG = generateCoarseGrid(Gt.parent, ones(Gt.parent.cells.num,1));
plotFaces(oG, 1:oG.faces.num,'FaceColor','none');
plotWell(Gt.parent, W,'FontSize',10);
view(-63, 50); axis tight; colorbar, clim([0 1-srw]); colormap(parula.^2);
hs     = [];
time   = cumsum([0; schedule.step.val])/year;
period = [1; schedule.step.control];
ptxt   = {'injection','migration'};

for i=1:numel(states)
    delete(hs)
    [h, h_max] = upscaledSat2height(states{i}.s(:,2), states{i}.sGmax, Gt, ...
                                    'pcWG', fluid.pcWG, ...
                                    'rhoW', fluid.rhoW, ...
                                    'rhoG', fluid.rhoG, ...
                                    'p', states{50}.pressure);
    sat = height2Sat(struct('h', h, 'h_max', h_max), Gt, fluid);
    title(sprintf('Time: %4d yrs (%s)', time(i),ptxt{period(i)}));
    ix = sat>0; if ~any(ix), continue; end
    %hs = plotCellData(Gt.parent, sat, ix); drawnow
end
xlim([0 42000]); ylim([0 42000]);
plotToolbar(Gt,states);colormap;colorbar;
save('LiuJiaGou_mrst_sol.mat','Gt','states');
% %% Trapping inventory
% % The result is a more detailed inventory that accounts for six different categories of CO2:
% %
% % # Structural residual - CO2 residually trapped inside a structural trap
% % # Residual - CO2 residually trapped outside any structural traps
% % # Residual in plume - fraction of the CO2 plume outside any structural
% %   traps that will be left behind as residually trapped droplets when the
% %   plume migrates away from its current position
% % # Structural plume - mobile CO2 volume that is currently contained within
% %   a residual trap;  if the containing structure is breached, this volume
% %   is free to migrate upward
% % # Free plume - the fraction of the CO2 plume outside of structural traps
% %   that is free to migrate upward and/or be displaced by imbibing brine.
% % # Exited - volume of CO2 that has migrated out of the domain through its
% %   lateral boundaries
% %
% % This model only has very small structural traps and residual trapping is
% % therefore the main mechanism.
% ta = trapAnalysis(Gt, false);
% reports = makeReports(Gt, states, model.rock, model.fluid, ...
%                       schedule, [srw, src], ta, []);
% 
% h1 = figure; plot(1); ax = get(h1, 'currentaxes');
% plotTrappingDistribution(ax, reports, 'legend_location', 'northwest');
