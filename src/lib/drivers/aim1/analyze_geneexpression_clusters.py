analyze_geneexpression_clusters (compose with combination of functions/variables){
	load gtex data
	if apply transformation? nil/ normalize / batch effect
	Quantify tissue association
	Store results
}

quantify tissue association
{
    Load clustering results
    silhouette
    method
    ARI
    SMC
    jaccard index
    https: // stats.stackexchange.com / questions / 95782 / what - are - the - most - common - metrics -
for -comparing - two - clustering - algorithms - especi
    https: //
    hdbscan.readthedocs.io / en / latest / comparing_clustering_algorithms.html
    Freeze an object of results and store
}
