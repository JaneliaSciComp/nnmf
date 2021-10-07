function cluster_sizes = get_cluster_sizes(spatial)
%GET_CLUSTER_SIZES Returns size of clusters

num_clusters = size(spatial,3);
cluster_sizes = zeros(num_clusters,1);

for i=1:num_clusters
    cluster_sizes(i) = nnz(spatial(:,:,i));
end

