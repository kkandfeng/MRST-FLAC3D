%% 
data = readtable("LiuJiaGou_tmp/LiuJiaGouRes_zone_center_out");
stress = readtable("flac to mrst/LiuJiaGou_stress_data");
pore_presssure = readtable("flac to mrst/LiuJiaGou_pore_pressure_data");
% Target coordinates
target_x = 15529;
target_y = 16107;
target_z = -3300;
tolerance = 4;
% Extract x, y, z coordinates
x = data.Var2;
y = data.Var3;
z = data.Var4;
% Calculate the Euclidean distance from each point to the target point
distances = sqrt((x - target_x).^2 + (y - target_y).^2 + (z - target_z).^2);
% Finds the index with the minimum distance
[~, closestIndex] = min(distances);
% Get information about the closest point
closestPoint = data(closestIndex, :);
% Use logical indexes to access the corresponding row data directly
esData = stress(stress.Var1 == closestPoint.Var1, :);
espp = pore_presssure(pore_presssure.Var1 == closestPoint.Var1, :);
%% 
% Extract the effective stress values
if ~isempty(esData)
es.max = esData.Var2;
es.int = esData.Var3;
es.min = esData.Var4;
fprintf('The total stress value found：\nmax: %f\nint: %f\nmin: %f\n', es.max, es.int, es.min);
else
fprintf('The total stress value corresponding to the nearest point was not found。\n');
end
if ~isempty(espp)
es.pp = espp.Var2;
fprintf('The pore pressure value found：\npp:  %f\n', es.pp);
else
fprintf('No pore pressure value was found for the nearest point。\n');
end
%% Calculate the average effective stress
es_mean = (es.max + es.int +es.min) / 3 - es.pp;
fprintf('The total stress value found：\nmean:  %f\n', es_mean);

%%  
% Define the parameters
phi_0 = ; % Porosity at zero stress
phi_r = ; % Residual porosity at high stress
sigma_M_prime = es_mean; % Mean effective stress
a = ; % Experimentally determined exponent for porosity
k_0 = ; % Permeability at zero stress
c = ; % Experimentally determined exponent for permeability
P_c0_Se = ; % Initial capillary pressure as a function of effective saturation

% Equation 19 - Calculate porosity
phi = phi_r + (phi_0 - phi_r) * exp(a * sigma_M_prime);

% Equation 20 - Calculate permeability
k = k_0 * exp(c * (phi / phi_0 - 1));

% Equation 21 - Calculate capillary pressure
P_c = P_c0_Se * sqrt(k_0 / phi_0) / sqrt(k / phi);

% Display the results
fprintf('Porosity (phi): %f\n', phi);
fprintf('Permeability (k): %f\n', k);
fprintf('Capillary Pressure (P_c): %f\n', P_c);