% Read flac_gp.txt file
data = load('flac_gp.txt');
x = data(:, 2); % X
y = data(:, 3); % Y
z = data(:, 4); % Z

% Select the node with the specified depth
target_z = 0;
tolerance = 1;
mask = abs(z - target_z) < tolerance;
x_selected = x(mask);
y_selected = y(mask);

% Filter approximate X and Y coordinates
x_unique = unique(x_selected);
y_unique = unique(y_selected);

new_x = x_unique;
new_y = y_unique;

% Filter for X and Y coordinates that are too close
for i = length(new_x)-1:-1:1
    if abs((new_x(i+1) - new_x(i)) / (new_x(i+1) + 1e-6)) < 1e-3
        new_x(i+1) = [];
    end 
end

for i = length(new_y)-1:-1:1
    if abs((new_y(i+1) - new_y(i)) / (new_y(i+1) + 1e-6)) < 1e-3
        new_y(i+1) = [];
    end
end
save('flac_gp_xy_out.mat', 'new_x','new_y');
