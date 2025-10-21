% === Clear workspace ===
clear; clc; close all;

% === Load Data ===
data = readtable('7_2.csv', 'VariableNamingRule', 'preserve');
x = data{:,'x(nm)'};
y = data{:,'y(nm)'};
z = data{:,'z(nm)'};
labels = data{:,'cluster label (-1 is noise)'};

% === Filter out noise
validIdx = labels ~= -1;
xValid = -x(validIdx);  % Flip X
yValid = y(validIdx);
zValid = z(validIdx);
xyz = [xValid, yValid, zValid];

% === K-means clustering on XY only
numClusters = 7;
xy = [xValid, yValid];
[idx_kmeans, xyCenters] = kmeans(xy, numClusters, 'Replicates', 5);

% === Estimate Z center and STD for each cluster
zCenters = zeros(numClusters, 1);
clusterSTD = zeros(numClusters, 3);
for c = 1:numClusters
    pts = xyz(idx_kmeans == c, :);
    zCenters(c) = mean(pts(:,3));
    clusterSTD(c,:) = std(pts, 0, 1);  % [stdX, stdY, stdZ]
end
clusterCenters = [xyCenters, zCenters];

% === Center around centroid mean
clusterCOM = mean(clusterCenters, 1);
xyzCentered = xyz - clusterCOM;
centroidsCentered = clusterCenters - clusterCOM;

% === Manual rotation and translation
angleZ = 25; shiftX = 0; shiftY = 0;
thetaZ = deg2rad(angleZ);
Rz = [cos(thetaZ), -sin(thetaZ), 0;
      sin(thetaZ),  cos(thetaZ), 0;
           0     ,       0     , 1];
xyzRotated = (Rz * xyzCentered')';
centroidsRotated = (Rz * centroidsCentered')';
xyzTransformed = xyzRotated + [shiftX, shiftY, 0];
centroidsTransformed = centroidsRotated + [shiftX, shiftY, 0];

% === Voxel Grid for Cloud Rendering
voxelSize = 1.5;
margin = 40;
xRange = min(xyzTransformed(:,1))-margin : voxelSize : max(xyzTransformed(:,1))+margin;
yRange = min(xyzTransformed(:,2))-margin : voxelSize : max(xyzTransformed(:,2))+margin;
zRange = min(xyzTransformed(:,3))-margin : voxelSize : max(xyzTransformed(:,3))+margin;
[Y, X, Z] = ndgrid(yRange, xRange, zRange);
vol = zeros(size(X));

% === Add Gaussian blobs at each cluster center
for c = 1:numClusters
    mu = centroidsTransformed(c,:);
    sigma = clusterSTD(c,:);
    sigma(sigma == 0) = 1;
    G = exp( -((X - mu(1)).^2 / (2*sigma(1)^2) + ...
               (Y - mu(2)).^2 / (2*sigma(2)^2) + ...
               (Z - mu(3)).^2 / (2*sigma(3)^2)) );
    vol = vol + G;
end
vol = vol / max(vol(:));  % Normalize

% === Save Orthographic Views (vector EPS)
plotPointCloudWithVoxel(xyzTransformed, centroidsTransformed, X, Y, Z, vol, [0, 90],  'view_XY_ortho');
plotPointCloudWithVoxel(xyzTransformed, centroidsTransformed, X, Y, Z, vol, [0, 0],   'view_XZ_ortho');
plotPointCloudWithVoxel(xyzTransformed, centroidsTransformed, X, Y, Z, vol, [90, 0],  'view_YZ_ortho');
plotPointCloudWithVoxel(xyzTransformed, centroidsTransformed, X, Y, Z, vol, [45, 25], 'view_3D_ortho');

% === Show Interactive 3D View (no export)
figure('Color', 'k'); hold on;
level = 0.3;
p = patch(isosurface(X, Y, Z, vol, level));
isonormals(X, Y, Z, vol, p);
p.FaceColor = [1 1 1];
p.FaceAlpha = 1;
p.EdgeColor = 'none';

scatter3(xyzTransformed(:,1), xyzTransformed(:,2), xyzTransformed(:,3), ...
    15, 'r', 'filled');
scatter3(centroidsTransformed(:,1), centroidsTransformed(:,2), centroidsTransformed(:,3), ...
    200, 'blue', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

xlabel('x', 'Color', 'w'); ylabel('y', 'Color', 'w'); zlabel('z', 'Color', 'w');
axis equal; pbaspect([1 1 1]); xlim padded; ylim padded; zlim padded;
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', ...
    'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out');
view(45, 25); rotate3d on;

% === Subfunction: Save vector EPS without transparency
function plotPointCloudWithVoxel(xyzPoints, centroids, X, Y, Z, vol, viewAngles, filenameBase)
    figure('Color', 'k'); hold on;

    % --- Solid Voxel Surface (no transparency)
    level = 0.3;
    p = patch(isosurface(X, Y, Z, vol, level));
    isonormals(X, Y, Z, vol, p);
    p.FaceColor = [1 1 1];
    p.FaceAlpha = 1;
    p.EdgeColor = 'none';

    scatter3(xyzPoints(:,1), xyzPoints(:,2), xyzPoints(:,3), 15, 'r', 'filled');
    scatter3(centroids(:,1), centroids(:,2), centroids(:,3), 200, 'blue', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

    xlabel('x', 'Color', 'w'); ylabel('y', 'Color', 'w'); zlabel('z', 'Color', 'w');
    axis equal; pbaspect([1 1 1]);
    xlim([-75 75]); ylim([-75 75]); zlim([-100 100]);
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', ...
        'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'out');
    camproj('orthographic'); view(viewAngles(1), viewAngles(2));
    set(gcf, 'InvertHardcopy', 'off');

    % --- Export to EPS (vector)
    epsFilename = [filenameBase '.eps'];
    print(gcf, epsFilename, '-depsc2', '-painters');
    close;
end
