clear
clc
close all

load('C:\Users\kyler\OneDrive\School\Colorado State\Research\Dr. Mueller\FEM\OOEIT_FEM_Driver\Results\Sbj009_Comp\Equal_Spacing\case142473-Cond-z0.133-BreathCycle-LBelt.mat')

top_nodes = nodes(sigma_GT{1},:);
top_sigma = sigma(sigma_GT{1},:);
bot_nodes = nodes(sigma_GT{2},:);
bot_sigma = sigma(sigma_GT{2},:);

frame = 1;
xt = top_nodes(:,1);
yt = top_nodes(:,2);

top_DT = delaunayTriangulation(xt, yt);
top_tri = top_DT.ConnectivityList;   % triangle vertex indices
top_pts = top_DT.Points;             % [x y] coordinates

xb = bot_nodes(:,1);
yb = bot_nodes(:,2);

bot_DT = delaunayTriangulation(xb, yb);
bot_tri = bot_DT.ConnectivityList;   % triangle vertex indices
bot_pts = bot_DT.Points;             % [x y] coordinates


for frame = [15, 20, 37]
    It = top_sigma(:,frame);
    Ib = bot_sigma(:,frame);
    figure()
        subplot(2,1,1)
            trisurf(top_tri, top_pts(:,1), top_pts(:,2), It, ...
                    'EdgeColor', 'k', ...      % no triangle edges
                    'EdgeAlpha', 0.5, ...
                    'FaceColor', 'interp');       % interpolate across vertices
            view(2);                 % look straight down
            axis equal off;
            set(gca, 'XDir', 'reverse');
            colorbar;
            clim([0, 0.8])
            colormap('jet')
            title(sprintf("Frame %d Ground Truth\nConductivity", frame))
        subplot(2,1,2)
            trisurf(bot_tri, bot_pts(:,1), bot_pts(:,2), Ib, ...
                    'EdgeColor', 'k', ...      % no triangle edges
                    'EdgeAlpha', 0.5, ...
                    'FaceColor', 'interp');       % interpolate across vertices
            view(2);                 % look straight down
            axis equal off;
            set(gca, 'XDir', 'reverse');
            colorbar;
            colormap('jet')
            clim([0, 0.8])
end

