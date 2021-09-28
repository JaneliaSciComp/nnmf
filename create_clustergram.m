function cg = create_clustergram(DFF)
%CREATE_CLUSTERGRAM Create a clustergram of DFF

cg = clustergram(DFF,'Cluster','Column','RowPDist','Correlation','ColorMap',jet,'OptimalLeafOrder','True');
cgAxis = plot(cg);
set(cgAxis,'Clim',[0,0.25])
end

