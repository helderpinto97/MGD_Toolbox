function dist=mgd_dist_qtl(y)
    qtl=quantile(y,[0.95 0.05]);
    dist=[qtl(1,:)-mean(y); mean(y)-qtl(2,:)];
end