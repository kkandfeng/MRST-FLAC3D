# MRST-FLAC3D
The program couples MRST and FLAC3D

The flac command is for version 7.0.

First, establish the flac3d model with init_flac.txt (generate the mesh with geomodel.txt), and compute until equilibrium.

Export the mesh as geo_model.flac3d, process to save only the node information to flac_gp.txt (node ID, X, Y, Z).

Process flac_gp.txt with flac_gp_xy.m to export the XY coordinates of nodes at a specified depth, obtaining flac_gp_xy_out.mat, for mesh generation in MRST.

Run ve_LiuJiaGou.m in the MRST environment, read in flac_gp_xy_out.mat, to generate the VE model of a single reservoir, lasting for 5 years.

After MRST computation, save Gt and states as mrst_sol.mat, and create a separate tmp folder.

In mrst_ve2vtk.m, specify mrst_sol.mat to write MRST results in vtk format stored in the tmp folder, saving nodes, elements, pressure, and saturation for easy later use.

Repeat the first three steps for the simulation of seepage in three reservoirs.

flac3d calls flac_zone_center.txt to export the element ID and center points of specified groups (LiuJiaGou, ShiQianFeng, and so on), saving to flac_zone_center_out.txt.

Run zone_cell_index.m to establish the correspondence between flac and MRST meshes, outputting zone_cell_index_out.txt, with the first column being flac element IDs, and the second column being MRST element IDs (starting from 1, for easier flac invocation later).

Invoke flac_import_mrst_ve_sol.txt to import the pressure and saturation calculated by MRST into the flac model for the next step of computation.

init_ming for mining analysis, using LiuJiaGou_pore_pressure_data.txt and LiuJiaGou_stress_data.txt to export pore pressure and stress, employing the flac_to_mrst.m file to update porosity, permeability, and capillary pressure.
