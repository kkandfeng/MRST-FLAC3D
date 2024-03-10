function mrst_ve2vtk()

    % Ensure the directories exist or create them if they don't

    % Load MRST data
    matlab_fname = 'LiuJiaGou_mrst_sol.mat';
    load(matlab_fname, 'Gt', 'states');

    % Extract grid data
    nodes_coords = Gt.nodes.coords;
    nodes_num = size(nodes_coords, 1);
    nodes_z = -Gt.nodes.z; % Convert depth to z-coordinate

    faces_nodePos = Gt.faces.nodePos;
    faces_nodes = Gt.faces.nodes;
    faces_num = Gt.faces.num;

    cells_faces = Gt.cells.faces;
    cells_facePos = Gt.cells.facePos;
    cells_num = Gt.cells.num;
    cells_centroids = Gt.cells.centroids;

    % Save cell centroids
    % dlmwrite('cell_center.txt', cells_centroids, ' ');
    writematrix(cells_centroids, 'cell_center.txt', 'Delimiter', ' ');

    % Save nodes
    nodes_xyz = [nodes_coords, nodes_z];
    % dlmwrite('nodes.txt', nodes_xyz, ' ');
    writematrix(nodes_xyz, 'nodes.txt', 'Delimiter', ' ');

    % Compute nodes for each cell
    nodes_each_cell = cellNodesList(cells_faces, cells_facePos, faces_nodes, faces_nodePos);
    % dlmwrite('elems.txt', nodes_each_cell, ' ');
    writematrix(nodes_each_cell, 'elems.txt', 'Delimiter', ' ');

    % Write VTK grid
    writeVTKGrid('grid.vtk', nodes_num, cells_num);

    % Read grid text for later use
    grid_txt = fileread('grid.vtk');

    % Extract and save solution data
    for i = 1:length(states)
        state = states{i};
        pressure = state.pressure;
        s = state.s;

        % dlmwrite(['pressure_', num2str(i-1), '.txt'], pressure, ' ');
        % dlmwrite(['sat_', num2str(i-1), '.txt'], s, ' ');
        writematrix(pressure, ['pressure_', num2str(i-1), '.txt'], 'Delimiter', ' ');
        writematrix(s, ['sat_', num2str(i-1), '.txt'], 'Delimiter', ' ');

        % Write VTK cell data
        writeVTKCellData(['step_', num2str(i-1), '.vtk'], cells_num, grid_txt, i-1);
    end
end