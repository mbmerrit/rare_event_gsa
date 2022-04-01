function [centroid_coordinates] = get_centroid_coordinates(nx, ny, xa, xb, ya, yb)
% Given a rectangular domain and structured FEM mesh, this code will return
% the x, y coordinates of all FE centroids. For nx, ny elements in each
% direction, you will get 2*nx*ny centroid objects

% start with a linspaced mesh
X = linspace(xa, xb, nx+1)'; Y = linspace(ya, yb, ny+1)';
% Here are the main formulas for the upper and lower centroids of each
% square-shaped FE cell
m1x = (3*X(1:end-1) + X(2:end))/4; m1y = (Y(1:end-1) + 3*Y(2:end))/4;
m2x = (X(1:end-1) + 3*X(2:end))/4; m2y = (3*Y(1:end-1) + Y(2:end))/4;
% This is just some bookkeping to make sure things are duplicated correctly
M1 = [repmat(m1x, ny, 1), repelem(m1y, nx)];
M2 = [repmat(m2x, ny, 1), repelem(m2y, nx)];
centroid_coordinates = [M1; M2];