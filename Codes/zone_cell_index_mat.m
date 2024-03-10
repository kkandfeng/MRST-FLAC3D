% Clear workspace and close all figures
clear; close all;

% Load data from text files
mrst_cell_cen = readmatrix('LiuJiaGou_tmp/cell_center.txt');
flac_zone_cen = readmatrix('LiuJiaGou_tmp/LiuJiaGouRes_zone_center_out.txt');
save_fname = 'LiuJiaGou_Res_zone_cell_index_out.txt';

ncell_0 = size(mrst_cell_cen, 1); 
ncell_1 = size(flac_zone_cen, 1); 

zone_cell_index = zeros(ncell_1, 2); %Initialize the index

% Add an index column to mrst_cell_cen
mrst_cell_cen = [(1:ncell_0)', mrst_cell_cen];
% disp(mrst_cell_cen)
% Sort flac_zone_cen by columns in the order: X, Y, Z
[~, idx] = sort(flac_zone_cen(:, 2)); %Sort the X columns and find the index
flac_sort_0 = flac_zone_cen(idx, :);

[~, idx] = sortrows(flac_sort_0, [3 2]);
flac_sort_1 = flac_sort_0(idx, :);

[~, idx] = sortrows(flac_sort_1, [4 3 2]);
flac_sort_2 = flac_sort_1(idx, :);

zone_cell_index(:, 1) = flac_sort_2(:, 1); % flac zone id

for i = 1:ncell_1
    zone_cell_index(i, 2) = mod(i-1, ncell_0) + 1; % corresponding mrst cell id
end

% Save the zone_cell_index to a text file
writematrix(zone_cell_index, save_fname, 'Delimiter', ' ');