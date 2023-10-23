## plotting function
plotting_func<-function(df_data,plot_name,sim){
  png(file=paste0(plot_name,sim,'.png'),width = 500,height = 300)
  p<-ggplot(data=df_data, aes(channel1, channel2, colour = factor(group)))+ 
    geom_point(size=4,show.legend = FALSE) +labs(x = "Green Channel",y='Red Channel')+theme(text = element_text(size = 50),
                                                                                            panel.grid.major = element_blank(),
                                                                                              panel.grid.minor = element_blank(),
                                                                                              panel.background = element_blank(),
                                                                                              axis.line = element_line(colour = "black"),
                                                                                              plot.margin = margin(t = 20,  # Top margin
                                                                                                                   r = 20,  # Right margin
                                                                                                                   b = 20,  # Bottom margin
                                                                                                                   l = 20))
  print(p)
  dev.off()
}


eval_func<-function(class_group,true_group){
  adj_randindex<-ARI(class_group,true_group)
  
  return(adj_randindex)
}


## gridding function
gridding_density<-function(data,s1=15,s2=15){
  range_G<-seq(min(data[,1])-abs(min(data[,1])/10),max(data[,1])+abs(min(data[,1])/10),length.out = s1)
  range_R<-seq(min(data[,2])-abs(min(data[,2])/10),max(data[,2])+abs(min(data[,2])/10),length.out = s2)
  
  intensity_mat<-matrix(0,ncol=s1,nrow=s2)
  for (i in 1:(s1-1)){
    for (j in 1:(s2-1)){
      intensity_mat[j,i]<-sum((data[,1]>=range_G[i]& data[,1]<range_G[i+1]) & (data[,2]>=range_R[j]& data[,2]<range_R[j+1])) 
    }
  }
  return (list(intensity_mat,range_G,range_R))
}


## search for each central points in each cluster
# central_points_func<-function(data,method='neighbourhood',s1=15,s2=15,perc=0.1,plotname=NULL,plot=FALSE){
#   unique_group<-unique(data$group)
#   data_index<-1:nrow(data)
#   outliers<-NULL
#   outliers_index<-NULL
#   non_outliers<-NULL
#   non_outliers_index<-NULL
#   
#   if (method=='neighbourhood'){
#     for (i in unique_group) {
#       data_group<-data[data$group==i,]
#       data_group_n<-nrow(data_group)
#       data_group_index<-data_index[data$group==i]
#       intensity_cal<-gridding_density(data_group,s1=s1,s2=s1)
#       range_G<-intensity_cal[[2]]
#       range_R<-intensity_cal[[3]]
#       unit_G<-range_G[2]-range_G[1]
#       unit_R<-range_R[2]-range_R[1]
#       neighbour_count<-NULL
#       
#       for (j in 1:data_group_n){
#         neighbour_count<-c(neighbour_count,sum((data_group[,1]>=(data_group[j,1]-unit_R))&(data_group[,1]<(data_group[j,1]+unit_R))&(data_group[,2]>=(data_group[j,2]-unit_G))&(data_group[,2]<(data_group[j,2]+unit_G))))
#       }
#       
#       data_group_reordered<-data_group[order(neighbour_count),]
#       data_group_index_reordered<-data_group_index[order(neighbour_count)]
#       
#       outliers<-rbind(outliers,data_group_reordered[1:round(data_group_n*perc),])
#       outliers_index<-c(outliers_index,data_group_index_reordered[1:round(data_group_n*perc)])
#       non_outliers<-rbind(non_outliers,data_group_reordered[(round(data_group_n*perc)+1):data_group_n,])
#       non_outliers_index<-c(non_outliers_index,data_group_index_reordered[(round(data_group_n*perc)+1):data_group_n])
#       
#       if(plot){
#         png(file=paste0(plotname,method,'.png'),width = 500,height = 300)
#         plot(outliers,col='red',xlim=c(min(data[,1])-10,max(data[,1])+10),ylim=c(min(data[,2])-10,max(data[,2])+10),xlab='Channel2', ylab='Channel3')
#         points(non_outliers)
#         dev.off()
#       }
#     }
#     
#   }
#   
#   return(list(outliers,outliers_index,non_outliers,non_outliers_index))
# }

central_points_func<-function(data,perc=0.1,plotname=NULL,plot=FALSE){
  unique_group<-unique(data$group)
  data_index<-1:nrow(data)
  outliers<-NULL
  outliers_index<-NULL
  non_outliers<-NULL
  non_outliers_index<-NULL
  
  data_sil<-silhouette_coef(data,data$group)[[1]]
  
  for (i in unique_group) {
    data_group<-data[data$group==i,]
    data_group_n<-nrow(data_group)
    data_group_index<-data_index[data$group==i]
    
    group_sil<-data_sil[data_sil$cluster==i,'width']
    data_group_reordered<-data_group[order(group_sil),]
    data_group_index_reordered<-data_group_index[order(group_sil)]
      
    outliers<-rbind(outliers,data_group_reordered[1:round(data_group_n*perc),])
    outliers_index<-c(outliers_index,data_group_index_reordered[1:round(data_group_n*perc)])
    non_outliers<-rbind(non_outliers,data_group_reordered[(round(data_group_n*perc)+1):data_group_n,])
    non_outliers_index<-c(non_outliers_index,data_group_index_reordered[(round(data_group_n*perc)+1):data_group_n])
      
    if(plot){
      png(file=paste0(plotname,method,'.png'),width = 500,height = 300)
      plot(outliers,col='red',xlim=c(min(data[,1])-10,max(data[,1])+10),ylim=c(min(data[,2])-10,max(data[,2])+10),xlab='Channel2', ylab='Channel3')
      points(non_outliers)
      dev.off()
    }
  }
    
  return(list(outliers,outliers_index,non_outliers,non_outliers_index))
}


euclidean.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

misclass_index<-function(data,class_group,true_centers){
  ## wrongly classified data points
  unique_data_group<-unique(class_group)
  unique_data_group<-unique_data_group[unique_data_group!=0]
  dist_mat<-NULL
  for (i in unique_data_group){
    if (!is.null(nrow(data[class_group==i,]))){
      group_mean<-apply(data[class_group==i,],2,mean)
    } else {group_mean<-data[class_group==i,]}
    dist<-NULL
    for (j in 1:nrow(true_centers)){
      dist<-c(dist,euclidean.dist(group_mean,true_centers[j,]))
    }
    dist_mat<-rbind(dist_mat,dist)
    
  }
  
  data_group_new<-class_group+100
  
  if (length(unique_data_group)<=nrow(true_centers)){
    assign_res<-solve_LSAP(dist_mat)
    
    for (k in 1:length(unique_data_group)){
       data_group_new[class_group==unique_data_group[k]]<-assign_res[k]
    }
    
  } else {
    assign_res<-solve_LSAP(t(dist_mat))

    for (k in 1:length(assign_res)){
      data_group_new[class_group==unique_data_group[assign_res[k]]]<-k
    }
  }
  
  return(data_group_new)
}


# silhouette_coef<-function(data,clustering,dist_method='euclidean',plot=FALSE,plot_name='orig',sim='orig'){
silhouette_coef<-function(data,clustering,plot=FALSE,plot_name='orig',sim='orig'){
  
  # si <- silhouette(clustering, dist(data, "euclidean"))
  ndim<-ncol(data)
  si <- approxSilhouette(data[,-ndim], clustering)
  
  sil_tab<-as.data.frame(si) %>% group_by(cluster) %>% summarise(mean_sil=round(mean(width),2))
  
  # find the mean of each cluster
  
  clust_pos<- data %>% group_by(group) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
  colnames(clust_pos)[1]<-'cluster'
  df_merged <- merge(x=sil_tab,y=clust_pos, by="cluster", all.x=TRUE)
  
  if(plot){
    png(file=paste0(plot_name,sim,'.png'),width = 500,height = 300)
    p<-ggplot(data=data, aes(channel1, channel2, colour = factor(group)))+ 
      geom_point(size=0.7,show.legend = FALSE) +labs(x = "Green Channel",y='Red Channel')+theme(panel.grid.major = element_blank(),
                                                                                                panel.grid.minor = element_blank(),
                                                                                                panel.background = element_blank(),
                                                                                                axis.line = element_line(colour = "black"),
                                                                                                plot.margin = margin(t = 20,  # Top margin
                                                                                                                     r = 20,  # Right margin
                                                                                                                     b = 20,  # Bottom margin
                                                                                                                     l = 20)) + annotate("text", x=df_merged$mean_x, y=df_merged$mean_y, label= df_merged$mean_sil)
    
    print(p)
    dev.off()
  }
  
  return(list(si,df_merged))
  
}

mahal_dist<-function(data){
  index_sel<-NULL
  index_all<-1:nrow(data)
  group_num<-length(unique(data$group))
  clust_pos<- data %>% group_by(group) %>% summarise(mean_x=mean(channel1),mean_y=mean(channel2))
  # Cutoff value for distances from Chi-Sqaure Dist. 
  cutoff <- qchisq(p = 0.99, df = ncol(data))
  
  for (i in 1:group_num){
    index_group<-index_all[data$group==i]
    data_group<-data[data$group==i,c(1,2)]
    data_group_cov<-cov(data_group)
    distances <- mahalanobis(x = data_group, center = as.numeric(clust_pos[clust_pos$group==i,c(2,3)]), cov = as.matrix(data_group_cov))
    
    ## Display observation whose distance greater than cutoff value
    index_sel<-c(index_sel,index_group[distances <= cutoff])
  } 
  return(index_sel)
}

simulation_plotcheck<-function(data,plot_name,sim){
  #png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(data[[1]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p2 <- ggplot(data[[2]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p3 <- ggplot(data[[3]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic() 
  
  p4 <- ggplot(data[[4]], aes(x,y)) +
    geom_point() + 
    #geom_density_2d() +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon") +
    theme_classic()
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  #print(p)
  #dev.off()
  ggsave(paste0(plot_name,sim,".png"),plot = p,width = 8,
         height = 8, dpi = 350)
}


simulation_plotcheck_depth<-function(data1,data2,plot_name,sim,method='Projection'){
  cluster1_mahadist<-mahalanobis(data1[[1]], colMeans(data1[[1]]), cov(data1[[1]]))
  cluster1_depth1<-DepthProc::depth(data1[[1]],data1[[1]],method = method)
  cluster1_depth2<-DepthProc::depth(data1[[1]],data2[[1]],method = method)
  cluster1_depth<-data.frame(OrigDepth=cluster1_depth1,SimDepth=cluster1_depth2,MahaDist=cluster1_mahadist)
  
  cluster2_mahadist<-mahalanobis(data1[[2]], colMeans(data1[[2]]), cov(data1[[2]]))
  cluster2_depth1<-DepthProc::depth(data1[[2]],data1[[2]],method = method)
  cluster2_depth2<-DepthProc::depth(data1[[2]],data2[[2]],method = method)
  cluster2_depth<-data.frame(OrigDepth=cluster2_depth1,SimDepth=cluster2_depth2,MahaDist=cluster2_mahadist)
  
  cluster3_mahadist<-mahalanobis(data1[[3]], colMeans(data1[[3]]), cov(data1[[3]]))
  cluster3_depth1<-DepthProc::depth(data1[[3]],data1[[3]],method = method)
  cluster3_depth2<-DepthProc::depth(data1[[3]],data2[[3]],method = method)
  cluster3_depth<-data.frame(OrigDepth=cluster3_depth1,SimDepth=cluster3_depth2,MahaDist=cluster3_mahadist)
  
  cluster4_mahadist<-mahalanobis(data1[[4]], colMeans(data1[[4]]), cov(data1[[4]]))
  cluster4_depth1<-DepthProc::depth(data1[[4]],data1[[4]],method = method)
  cluster4_depth2<-DepthProc::depth(data1[[4]],data2[[4]],method = method)
  cluster4_depth<-data.frame(OrigDepth=cluster4_depth1,SimDepth=cluster4_depth2,MahaDist=cluster4_mahadist)
  
  p1 <- ggplot(cluster1_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p2 <- ggplot(cluster2_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p3 <- ggplot(cluster3_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  p4 <- ggplot(cluster4_depth, aes(x = OrigDepth, y = SimDepth, colour = MahaDist)) +
    geom_point()+ geom_abline(intercept = 0, slope = 1)+
    theme_classic() 
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  
  ggsave(paste0(plot_name,sim,".png"),plot = p,width = 10,
         height = 8, dpi = 350)
  
}


approxSilhouette <- function(x, clusters) {
  x <- as.matrix(x)
  uclust <- sort(unique(clusters))
  averaged <- list(length(uclust))
  clust.var <- numeric(length(uclust))
  
  for (i in seq_along(uclust)) {
    current <- uclust[i]==clusters
    xcurrent <- x[current,,drop=FALSE]
    centroid <- colMeans(xcurrent)
    averaged[[i]] <- centroid
    clust.var[i] <- sum(colMeans(sweep(xcurrent, 2, centroid)^2))
  }
  
  self.dist <- other.dist <- rep(Inf, nrow(x))
  other.clust <- integer(nrow(x))
  tx <- t(x)
  
  for (i in seq_along(uclust)) {
    D <- sqrt(colSums((tx - averaged[[i]])^2) + clust.var[i])
    
    is.self <- uclust[i]==clusters
    self.dist[is.self] <- D[is.self]
    
    is.other <- !is.self
    other.D <- D[is.other]
    better <- other.D < other.dist[is.other]
    other.dist[is.other][better] <- other.D[better]
    other.clust[is.other][better] <- i
  }
  
  result<-data.frame(
    cluster=clusters,
    other=uclust[other.clust],
    width=(other.dist - self.dist)/pmax(other.dist, self.dist),
    row.names=rownames(x)
  )
  
  return(result)
}


## silhouette coefficient plots
silhouette_plot<-function(sil_data,plot_name,sim){
  cluster_sil<-list()
  for (i in 1:4){
    cluster_sil_tmp<-NULL
    for (j in 1:100){
      cluster_sil_tmp<-c(cluster_sil_tmp,sil_data[[j]][i,2])
    }
    cluster_sil[[i]]<-data.frame(silhouette_coef=cluster_sil_tmp)
  }
  
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(cluster_sil[[4]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  
  p2 <- ggplot(cluster_sil[[3]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p3 <- ggplot(cluster_sil[[1]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p4 <- ggplot(cluster_sil[[2]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  print(p)
  dev.off()
}

silhouette_single_plot<-function(sil_data,plot_name,sim){
  cluster_sil<-list()
  for (i in 1:4){
    cluster_sil[[i]]<-data.frame(silhouette_coef=sil_data[sil_data[,1]==i,3])
  }
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 600)
  p1 <- ggplot(cluster_sil[[4]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  
  p2 <- ggplot(cluster_sil[[3]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p3 <- ggplot(cluster_sil[[1]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  p4 <- ggplot(cluster_sil[[2]],aes(x=silhouette_coef)) + geom_histogram() + theme_classic()+labs(x = "silhouette coefficient")
  
  prow <- plot_grid (p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",plot.margin = margin(l=0,r = 0)),
                     p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p3 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     p4 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position="none",axis.title.y = element_blank(),plot.margin = margin(r = 0)),
                     align = 'vh',
                     labels = c("A","B","C","D"),
                     hjust = 0,
                     nrow = 2
  )
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, ncol = 1, rel_heights = c(1.5, 0.2))
  print(p)
  dev.off()
}

eval_plot<-function(result_list,plot_name,sim){
  ### argument: 
  ### result_list should be a list of evaluation results, such as rand index
  
  methods_name<-c('kmeans','kmeans_initials','cmeans','cmeans_initials','dbscan', 'fsom','flowpeaks','flowclust',
                  'flowclust_auto','flowmerge_moreclust','samspectral','samspectral_auto')
  
  result_df<-NULL
  for (i in 1:length(methods_name)){
    result_df<-rbind(result_df,cbind(result_list[[methods_name[i]]],methods_name[i]))
  }
  
  result_df<-data.frame(result_df)
  colnames(result_df)<-c('ri','adj_ri','ami','methods')
  result_df$adj_ri<-as.numeric(result_df$adj_ri)
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 400)
  
  p<-ggplot(result_df, aes(x=methods, y=adj_ri)) +
    geom_boxplot()+theme_classic()+labs(x = "Methods",y='Adjusted Rand Index')
  
  print(p)
  dev.off()
}


sens_plot<-function(result_list,true_lambda,plot_name,sim){
  ### argument: 
  ### result_list should be a list of 
  
  methods_name<-c('kmeans','kmeans_initials','cmeans','cmeans_initials','dbscan', 'fsom','flowpeaks','flowclust',
                  'flowclust_auto','flowmerge_moreclust','samspectral','samspectral_auto')
  
  result_df<-NULL
  for (i in 1:length(methods_name)){
    result_df<-rbind(result_df,cbind((result_list[[methods_name[i]]]-true_lambda)/true_lambda,methods_name[i]))
  }
  
  result_df<-data.frame(result_df)
  colnames(result_df)<-c('rel_bias','methods')
  result_df$rel_bias<-as.numeric(result_df$rel_bias)
  result_df<-result_df[!is.na(result_df$rel_bias),]
  
  png(file=paste0(plot_name,sim,'.png'),width = 1000,height = 400)
  
  p<-ggplot(result_df, aes(x=methods, y=rel_bias)) +
    geom_boxplot()+theme_classic()+labs(x = "methods",y='relative bias')+ylim(-1, 1)+
    geom_hline(yintercept=0, linetype="dashed", color = "red")
  
  print(p)
  dev.off()
}

