clear
clc

circ = load('Circle.mat');
spi = load('Spiral.mat');
circ = circ.X; % 2 columns (x ; y)
spi = spi.X; % 3 columns (x ; y ; label) 
w_circ = exp(-squareform(pdist(circ)).^2 / 2); % Sigma=1
w_spi = exp(-squareform(pdist(spi(:,1:2))).^2 / 2); % Sigma=1

%% Graphs Analysis

k = input('Insert k for KNN similarity matrix [Recommended: 10 - 20 - 40]: '); % Choose a parameter for KNN similarity matrix building
tol = input('Insert tolerance (just as a negative power of ten, for example 12 means 1e-12) to estimate eigenvalues [Recommended 10]: ');
tol = 10.^(-tol);
maxIter = input('Insert maximum number of iterations to operate for estimating eigenvalues [Recommended 500]: ');

tic;

w_circ = knn_sim(w_circ,k);
w_spi = knn_sim(w_spi,k);
w_circ = triu(w_circ) + triu(w_circ, 1)'; % Using 'upper' to force simmetry
w_spi = triu(w_spi) + triu(w_spi, 1)'; % Using 'upper' to force simmetry
w_circ(logical(eye(size(w_circ)))) = 0; % Setting 0 on the diagonal
w_spi(logical(eye(size(w_spi)))) = 0; % Setting 0 on the diagonal

graph_circ = graph(w_circ); 
graph_spi = graph(w_spi); 
figure(1);
subplot(1,2,1)
plot(graph_circ,'XData',circ(:,1),'YData',circ(:,2))
subplot(1,2,2)
plot(graph_spi,'XData',spi(:,1),'YData',spi(:,2))

%% Spectral Clustering

D_circ = diag(sum(w_circ,2));
D_spi = diag(sum(w_spi,2));

L_circ = D_circ - w_circ;
L_spi = D_spi - w_spi;

aux_circ = D_circ^(-1/2);
aux_spi = D_spi^(-1/2);
L_circ_sym = eye(length(L_circ(1,:))) - aux_circ * w_circ * aux_circ;
L_spi_sym = eye(length(L_spi(1,:))) - aux_spi * w_spi * aux_spi;

clc
n_eigenval_to_find = input('Insert number of singular values to estimate [Recommended 5]: '); % Number of smallest eigenvalues to find with inverse method power

[egvec_circ_ex,egval_circ_ex] = eig(L_circ);
egval_circ_ex = diag(egval_circ_ex);
[egvec_spi_ex,egval_spi_ex] = eig(L_spi);
egval_spi_ex = diag(egval_spi_ex);

[egvec_circ_sym_ex,egval_circ_sym_ex] = eig(L_circ_sym);
egval_circ_sym_ex = diag(egval_circ_sym_ex);
[egvec_spi_sym_ex,egval_spi_sym_ex] = eig(L_spi_sym);
egval_spi_sym_ex = diag(egval_spi_sym_ex);

[egvec_circ,egval_circ] = inv_pow(L_circ,n_eigenval_to_find,tol,maxIter);
[egvec_spi,egval_spi] = inv_pow(L_spi,n_eigenval_to_find,tol,maxIter);

[egvec_circ_sym,egval_circ_sym] = inv_pow(L_circ_sym,n_eigenval_to_find,tol,maxIter);
[egvec_spi_sym,egval_spi_sym] = inv_pow(L_spi_sym,n_eigenval_to_find,tol,maxIter);

egval_circ_ex_s = sort(egval_circ_ex);
egval_spi_ex_s = sort(egval_spi_ex);
egval_spi_sym_ex_s = sort(egval_spi_sym_ex);
egval_circ_sym_ex_s = sort(egval_circ_sym_ex);

figure(2);
subplot(2,2,1);
plot(1:length(egval_circ),egval_circ,LineWidth=2,LineStyle="--");
hold on
plot(1:length(egval_circ),egval_circ_ex_s(1:n_eigenval_to_find));
title('Circle')
hold off
subplot(2,2,2);
plot(1:length(egval_circ),egval_spi,LineWidth=2,LineStyle="--");
hold on
plot(1:length(egval_circ),egval_spi_ex_s(1:n_eigenval_to_find));
title('Spiral')
hold off
subplot(2,2,3);
plot(1:length(egval_circ),egval_circ_sym,LineWidth=2,LineStyle="--");
hold on
plot(1:length(egval_circ),egval_circ_sym_ex_s(1:n_eigenval_to_find));
title('Circle - L Sym.')
hold off
subplot(2,2,4);
plot(1:length(egval_circ),egval_spi_sym,LineWidth=2,LineStyle="--");
hold on
plot(1:length(egval_circ),egval_spi_sym_ex_s(1:n_eigenval_to_find));
title('Spiral - L Sym.')
hold off

clc
m_circ = input('Observe the plot and choose a number of clusters for Circle [Recommended 3]: ');
m_spir = input('Observe the plot and choose a number of clusters for Spiral [Recommended 3]: ');
m_circ_sym = input('Observe the plot and choose a number of clusters for Circle - L Sym. [Recommended 3]: ');
m_spir_sym = input('Observe the plot and choose a number of clusters for Spiral - L Sym. [Recommended 3]: ');


U_circ = extr_smallest_eigenv(egvec_circ,egval_circ,m_circ);
U_spir = extr_smallest_eigenv(egvec_spi,egval_spi,m_spir);

U_circ_sym = extr_smallest_eigenv(egvec_circ_sym,egval_circ_sym,m_circ_sym);
U_spir_sym = extr_smallest_eigenv(egvec_spi_sym,egval_spi_sym,m_spir_sym);
row_norms_circ = sqrt(sum(U_circ_sym.^2, 2));
U_circ_sym = U_circ_sym ./ row_norms_circ;
row_norms_spir = sqrt(sum(U_spir_sym.^2, 2));
U_spir_sym = U_spir_sym ./ row_norms_spir;

iteropt = statset('MaxIter', 2000);

cluster_circ = kmeans(U_circ,m_circ,'Options',iteropt,'Replicates',10);
cluster_spi = kmeans(U_spir,m_spir,'Options',iteropt,'Replicates',10);

cluster_circ_sym = kmeans(U_circ_sym,m_circ_sym,'Options',iteropt,'Replicates',10);
cluster_spi_sym = kmeans(U_spir_sym,m_spir_sym,'Options',iteropt,'Replicates',10);

U_circ_ex = extr_smallest_eigenv(egvec_circ_ex,egval_circ_ex,m_circ);
U_spir_ex = extr_smallest_eigenv(egvec_spi_ex,egval_spi_ex,m_spir);

U_circ_sym_ex = extr_smallest_eigenv(egvec_circ_sym_ex,egval_circ_sym_ex,m_circ_sym);
U_spir_sym_ex = extr_smallest_eigenv(egvec_spi_sym_ex,egval_spi_sym_ex,m_spir_sym);
row_norms_circ = sqrt(sum(U_circ_sym_ex.^2, 2));
U_circ_sym_ex = U_circ_sym_ex ./ row_norms_circ;
row_norms_spir = sqrt(sum(U_spir_sym_ex.^2, 2));
U_spir_sym_ex = U_spir_sym_ex ./ row_norms_spir;

cluster_circ_ex = kmeans(U_circ_ex,m_circ,'Replicates',2);
cluster_spi_ex = kmeans(U_spir_ex,m_spir,'Replicates',2);

cluster_circ_sym_ex = kmeans(U_circ_sym_ex,m_circ_sym,'Replicates',2);
cluster_spi_sym_ex = kmeans(U_spir_sym_ex,m_spir_sym,'Replicates',2);

cluster_circ_control = kmeans(circ,m_circ_sym,'Replicates',2);

time = toc;

%% Results Representation
clc 

disp(['Using KNN with k: ', num2str(k)])

comp = conncomp(graph_circ);
n_comp = max(comp);
disp(['Number of Circle connected components: ', num2str(n_comp)]);

comp = conncomp(graph_spi);
n_comp = max(comp);
disp(['Number of Spiral connected components: ', num2str(n_comp)]);

disp(['Time Elapsed: ', num2str(time)]);

figure(3);
subplot(2,4,1);
scatter(circ(:,1), circ(:,2), [], cluster_circ, 'filled')
title('Circle - Est. Eig.')
subplot(2,4,2);
scatter(spi(:,1), spi(:,2), [], cluster_spi, 'filled')
title('Spiral - Est. Eig.')

subplot(2,4,3);
scatter(circ(:,1), circ(:,2), [], cluster_circ_sym, 'filled')
title('Circle - Norm. L - Est. Eig.')
subplot(2,4,4);
scatter(spi(:,1), spi(:,2), [], cluster_spi_sym, 'filled')
title('Spiral - Norm. L - Est. Eig.')

subplot(2,4,5);
scatter(circ(:,1), circ(:,2), [], cluster_circ_ex, 'filled')
title('Circle - Ex. Eig.')
subplot(2,4,6);
scatter(spi(:,1), spi(:,2), [], cluster_spi_ex, 'filled')
title('Spiral - Ex. Eig.')

subplot(2,4,7);
scatter(circ(:,1), circ(:,2), [], cluster_circ_sym_ex, 'filled')
title('Circle - Norm. L - Ex. Eig.')
subplot(2,4,8);
scatter(spi(:,1), spi(:,2), [], cluster_spi_sym_ex, 'filled')
title('Spiral - Norm. L - Ex. Eig.')

figure(4);
subplot(1,2,1);
scatter(circ(:,1), circ(:,2), [], cluster_circ_control, 'filled')
title('Circle - KMeans')
subplot(1,2,2);
scatter(spi(:,1), spi(:,2), [], spi(:,3), 'filled')
title('Spiral - KMeans')

%% Functions
function w = knn_sim(mat,k) % Create a matrix with the k-nearest neighbours of each node, putting 0 on others
    k = k+1; % To not consider the point itself
    [nRows, nCols] = size(mat);
    w = zeros(nRows, nCols);
    for i = 1:nRows
        [~, idx] = maxk(mat(i, :), k); 
        w(i, idx) = mat(i, idx);
    end 
end

function U = extr_smallest_eigenv(eigenvec,eigenval,m) % Extracts the m eigenvectors related to the m smallest eigenvalues
    non_zero_indices = find(eigenval ~= 0);
    non_zero_values = eigenval(non_zero_indices);
    [~, sorted_indices] = mink(non_zero_values, m);
    original_indices = non_zero_indices(sorted_indices);
    U = eigenvec(:,original_indices);
end

function [eigvec,eigval] = inv_pow(A,k,tol,maxIter) % Estimates the k smallest eigenvalues and their eigenvectors of a matrix A
    n = size(A, 1); % Size of the matrix
    eigval = zeros(k, 1); % To store eigenvalues
    eigvec = zeros(n, k); % To store eigenvectors
    aux_A = A;
    for i = 1:k
        v = rand(n, 1); % Random initial vector
        v = v ./ norm(v); % Normalize
        lambda_old = inf; % Initialize previous eigenvalue
        for iter = 1:maxIter
            v_k1 = aux_A \ v; % Compute new A^k * v
            lambda = v' * v_k1; 
            v = v_k1 / norm(v_k1); % Normalize eigenvector
            if abs(lambda - lambda_old) < tol % Check convergence
                break;
            end
            lambda_old = lambda; % Update eigenvalue
        end
        eigval(i) = abs(1/lambda); % Compute actual eigenvalue
        eigvec(:, i) = v; % Store eigenvector
        aux_A = aux_A + (v * v'); % Sliding the eigenvalue found up to a random relatively big number to find new smallest eigenvalue in next iteration
    end
end
