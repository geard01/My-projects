clear
clc

%% Loading data

circ = load('C:\Users\geard\OneDrive\Download\Computational\HW2\Circle.mat');
spiral = load('C:\Users\geard\OneDrive\Download\Computational\HW2\Spiral.mat');
circ = circ.X;
spiral = spiral.X;
sigma = 1;

%% Point 1
%Create the graph from the data

Sim_circ = exp((-squareform(pdist(circ).^2))/(2*sigma));
Sim_spiral = exp((-squareform(pdist(spiral(:,1:2)).^2))/(2*sigma));

k = 40;
W_circ = knn_sim_graph(Sim_circ, k);
W_spiral = knn_sim_graph(Sim_spiral, k);

graph_circ = graph(W_circ);
graph_spiral = graph(W_spiral);

plot(graph_circ,'XData',circ(:,1),'YData',circ(:,2))
figure; 
plot(graph_spiral,'XData',spiral(:,1),'YData',spiral(:,2))

%% Point 2
%Compute D, W and L in sparse form

D_circ = sparse(diag(sum(W_circ, 2)));
D_spiral = sparse(diag(sum(W_spiral, 2)));

W_circ = sparse(W_circ);
W_spiral = sparse(W_spiral);

L_circ = sparse(D_circ - W_circ);
L_spiral = sparse(D_spiral - W_spiral);

%% Point 3
%Count the number of connected components

comp_circ = conncomp(graph_circ);
n_comp_circ = max(comp_circ);

comp_spiral = conncomp(graph_spiral);
n_comp_spiral = max(comp_spiral);

%% Point 4
%Coumpute eigenvalues with inverse power function and choose how many cluster by looking at the smallest eigenvalues

tol = 1e-6;       % Tolerance for convergence
max_iter = 1000;  % Maximum number of iterations per eigenvalue

m_circ = 10;
[eig_val_circ, ~] = inverse_power_method(L_circ,m_circ,tol,max_iter);

m_spiral = 10;
[eig_val_spiral,~] = inverse_power_method(L_spiral,m_spiral,tol,max_iter);

m_circ = 3;
m_spiral = 3;

%% Point 5
%Compute the eigenvectors of the smallest eigenvalues

[~, eig_vec_circ] = inverse_power_method(L_circ,m_circ,tol,max_iter);
[~, eig_vec_spiral] = inverse_power_method(L_spiral,m_spiral,tol,max_iter);
U_circ = eig_vec_circ;
U_spiral = eig_vec_spiral;

%% Point 6
%Normalize the U matrixes and apply the kmean algorithm

row_norms_circ = sqrt(sum(U_circ.^2, 2)); 
U_circ = U_circ ./ row_norms_circ; 

row_norms_spiral = sqrt(sum(U_spiral.^2, 2)); 
U_spiral = U_spiral ./ row_norms_spiral; 

idx_circ = kmeans(U_circ,m_circ,'Replicates',10);
idx_spiral = kmeans(U_spiral,m_spiral, 'Replicates',10);

%% Point 7
%Assign the points of the datasets to cell arrays

Clu_circ = cell(m_circ, 1);  
Clu_spiral = cell(m_spiral, 1); 

for i = 1:size(circ, 1)
    cluster_idx = idx_circ(i);
    Clu_circ{cluster_idx} = [Clu_circ{cluster_idx}; circ(i, :)];  % Add the original point circ(i, :) to the corresponding cluster
end

for i = 1:size(spiral, 1)
    cluster_idx = idx_spiral(i); 
    Clu_spiral{cluster_idx} = [Clu_spiral{cluster_idx}; spiral(i, 1:2)];  % Add the original point spir(i, :) to the corresponding cluster
end

%% Point 8
% Plot both sets of point colored by the cluster they belong to

figure;
subplot(1, 2, 1); 
hold on;

colormap_circ = lines(m_circ); 

for cluster_idx = 1:m_circ
    cluster_points = Clu_circ{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_circ(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Circular Data Clusters');
legend('show');
hold off;

subplot(1, 2, 2); 
hold on;

colormap_spiral = lines(m_spiral); 

for cluster_idx = 1:m_spiral
    cluster_points = Clu_spiral{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_spiral(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Spiral Data Clusters');
legend('show');
hold off;

%% Point 9
%Applying the kmeans directly to the dataset without the spectral
%clustering understanding the differences

idx_circ_bas = kmeans(circ,m_circ,'Replicates',10);
idx_spiral_bas = kmeans(spiral,m_spiral, 'Replicates',10);


Clu_circ_bas = cell(m_circ, 1);  
Clu_spiral_bas = cell(m_spiral, 1); 

for i = 1:size(circ, 1)
    cluster_idx = idx_circ_bas(i);
    Clu_circ_bas{cluster_idx} = [Clu_circ_bas{cluster_idx}; circ(i, :)];  % Add the original point circ(i, :) to the corresponding cluster
end

for i = 1:size(spiral, 1)
    cluster_idx = idx_spiral_bas(i); 
    Clu_spiral_bas{cluster_idx} = [Clu_spiral_bas{cluster_idx}; spiral(i, 1:2)];  % Add the original point spir(i, :) to the corresponding cluster
end

figure;
subplot(1, 2, 1); 
hold on;

colormap_circ = lines(m_circ); 

for cluster_idx = 1:m_circ
    cluster_points = Clu_circ_bas{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_circ(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Circular without spectral clustering');
legend('show');
hold off;

subplot(1, 2, 2); 
hold on;

colormap_spiral = lines(m_spiral); 

for cluster_idx = 1:m_spiral
    cluster_points = Clu_spiral_bas{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_spiral(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Spiral without spectral clustering');
legend('show');
hold off;

%% Point optional 1
% Generate a synthetic 3D data and perform the same process

rng(43);
num_clusters = 3;
data1 = mvnrnd([1, 1, 1], eye(3), 100); % Cluster 1 (3D Gaussian)
data2 = mvnrnd([5, 5, 5], eye(3), 100); % Cluster 2 (3D Gaussian)
data3 = mvnrnd([9, 1, 1], eye(3), 100); % Cluster 3 (3D Gaussian)
data = [data1; data2; data3];

figure;
scatter3(data(:,1), data(:,2), data(:,3), 36, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Synthetic 3D Data');

%Same process
Sim_3d = exp((-squareform(pdist(data).^2))/(2*sigma));
W_3d = knn_sim_graph(Sim_3d, k);
graph_3d = graph(W_3d);
D_3d = sparse(diag(sum(W_3d, 2)));
W_3d = sparse(W_3d);
L_3d = sparse(D_3d - W_3d);
m_3d=3;
[~, eig_vec_3d] = inverse_power_method(L_3d,m_3d,tol,max_iter);
U_3d = eig_vec_3d;
row_norms_3d = sqrt(sum(U_3d.^2, 2)); 
U_3d = U_3d ./ row_norms_3d; 
[idx_3d, C_3d] = kmeans(U_3d, m_3d, 'Replicates', 10);
Clu_3d = cell(m_3d, 1); 

for i = 1:size(data, 1)
    cluster_idx = idx_3d(i);
    Clu_3d{cluster_idx} = [Clu_3d{cluster_idx}; data(i, :)];  
end

figure;
hold on;

colors = lines(m_3d);

for cluster_idx = 1:m_3d
    cluster_points = Clu_3d{cluster_idx};
    scatter3(cluster_points(:, 1), cluster_points(:, 2), cluster_points(:, 3), 36, colors(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Spectral Clustering');
legend('show');
hold off;

%% Point optional 2
%Same process again but with the 2 L_sym matrixes

L_sym_circ = eye(length(L_circ(1,:))) - (D_circ^(-1/2)) * W_circ * (D_circ^(-1/2)); % Normalized symmetric Laplacian
[~, eig_vec_circ_sym] = inverse_power_method(L_sym_circ, m_circ, tol, max_iter);
U_circ_sym = eig_vec_circ_sym;
row_norms_circ_sym = sqrt(sum(U_circ_sym.^2, 2)); 
U_circ_sym = U_circ_sym ./ row_norms_circ_sym; 
idx_circ_sym= kmeans(U_circ_sym, m_circ, 'Replicates', 10); 

Clu_circ_sym = cell(m_circ, 1);
for i = 1:size(circ, 1)
    cluster_idx = idx_circ_sym(i);
    Clu_circ_sym{cluster_idx} = [Clu_circ_sym{cluster_idx}; circ(i, :)];
end

L_sym_spiral = eye(length(L_spiral(1,:))) - (D_spiral^(-1/2)) * W_spiral * (D_spiral^(-1/2)); % Normalized symmetric Laplacian
[~, eig_vec_spiral_sym] = inverse_power_method(L_sym_spiral, m_spiral, tol, max_iter); 
U_spiral_sym = eig_vec_spiral_sym;
row_norms_spiral_sym = sqrt(sum(U_spiral_sym.^2, 2)); 
U_spiral_sym = U_spiral_sym ./ row_norms_spiral_sym; 
idx_spiral_sym = kmeans(U_spiral_sym, m_spiral, 'Replicates', 10); 

Clu_spiral_sym = cell(m_spiral, 1);
for i = 1:size(spiral, 1)
    cluster_idx = idx_spiral_sym(i);
    Clu_spiral_sym{cluster_idx} = [Clu_spiral_sym{cluster_idx}; spiral(i, 1:2)];
end

figure;
subplot(2, 2, 1); 
hold on;

colormap_circ = lines(m_circ);

for cluster_idx = 1:m_circ
    cluster_points = Clu_circ{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_circ(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Circular Data Clusters');
legend('show');
hold off;

subplot(2, 2, 2);
hold on;

colormap_spiral = lines(m_spiral);

for cluster_idx = 1:m_spiral
    cluster_points = Clu_spiral{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colormap_spiral(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Spiral Data Clusters');
legend('show');
hold off;

subplot(2, 2, 3); 
hold on;

colors = lines(m_circ);

for cluster_idx = 1:m_circ
    cluster_points = Clu_circ_sym{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Circular Data Clustering with L_{sym}');
legend('show');
hold off;

subplot(2, 2, 4);
hold on;

colors = lines(m_spiral);

for cluster_idx = 1:m_spiral
    cluster_points = Clu_spiral_sym{cluster_idx};
    scatter(cluster_points(:, 1), cluster_points(:, 2), 36, colors(cluster_idx, :), 'filled', 'DisplayName', sprintf('Cluster %d', cluster_idx));
end

xlabel('X');
ylabel('Y');
title('Spiral Data Clustering with L_{sym}');
legend('show');
hold off;

%% Point Optional 3
% Create the inverse power method function
function [eigval, eigvec] = inverse_power_method(B, M, tol, max_iter)

    n = size(B, 1);  
    eigval = zeros(M, 1);    
    eigvec = zeros(n, M);

    for i = 1:M
        x = rand(n, 1);  
        x = x / norm(x);  
        lambda_old = 0;   

        % Inverse Power Method to find smallest eigenvalue
        for iter = 1:max_iter
            v = B \ x; 
            lambda = x' * v;  
            x = v / norm(v);      
            if abs(lambda - lambda_old) < tol
                break;
            end           
            lambda_old = lambda;  
        end

        eigval(i) = 1 / lambda;  
        eigvec(:, i) = x;       
        % Deflation: Remove the contribution of the found eigenvalue/eigenvector
        B = B \ eye(size(B))  - lambda * (x * x');  % Update the matrix by subtracting rank-1 update
        B= B \ eye(size(B));
    end
end


%% Functions

function W_knn = knn_sim_graph(Sim, k)
    n = size(Sim, 1);
    W_knn = zeros(n,n);
    for i = 1:n
        [~, indices] = sort(Sim(i, :), 'descend');
        neighbors = indices(2:k+1);
        W_knn(i, neighbors) = Sim(i, neighbors);
    end
    W_knn = max(W_knn,W_knn');
end
