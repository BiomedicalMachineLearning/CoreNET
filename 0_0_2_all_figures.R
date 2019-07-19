
library(zinbwave)
library(DESeq2)
library(edgeR)
library(UpSetR)
library(gplots)
library(ggraph)
library(igraph)
library(RColorBrewer)
require(VennDiagram)

run_deseq2<-function(data,clusters,N,list_cells_1,list_cells_2){
  
  data_1<-do.call('cbind',lapply(list_cells_1,function(X){
    data[,which(clusters==unique(clusters)[X])]
  }))
  colnames(data_1)<-rep("one",dim(data_1)[2])
  
  data_2<-do.call('cbind',lapply(list_cells_2,function(X){
    data[,which(clusters==unique(clusters)[X])]
  }))
  colnames(data_2)<-rep("two",dim(data_2)[2])
  
  cts<-cbind(data_1,data_2)
  coldata=DataFrame(as.integer(as.factor(colnames(cts))),row.names = as.factor(colnames(cts)))
  
  zinb<-SummarizedExperiment(assays=list(counts=cts), colData=coldata,rowData=row.names(cts))
  zinb$condition<-factor(colnames(cts))
  
  coldata_zinb=data.frame(as.factor(colnames(cts)))
  colnames(coldata_zinb)<-"condition"
  zinb <- zinbwave(zinb, K=0, BPPARAM=SerialParam(), epsilon=1e12)
  dds<-DESeqDataSetFromMatrix(zinb@assays$data$counts,colData = coldata_zinb,
                              design=~condition)
  dds <- DESeq(dds, test="LRT", reduced=~1,
               sfType="poscounts", minmu=1e-6, minRep=Inf)
  res <- results(dds)
  return(res)
}



run_edger<-function(data,clusters,N,list_cells_1,list_cells_2){
  
  data_1<-do.call('cbind',lapply(list_cells_1,function(X){
    data[,which(clusters==unique(clusters)[X])]
  }))
  colnames(data_1)<-rep("one",dim(data_1)[2])
  
  data_2<-do.call('cbind',lapply(list_cells_2,function(X){
    data[,which(clusters==unique(clusters)[X])]
  }))
  colnames(data_2)<-rep("two",dim(data_2)[2])

  x<-cbind(data_1,data_2)
  
  group <- as.factor(colnames(x))
  y <- DGEList(counts=x,group=group)
  y <- calcNormFactors(y)
  design <- model.matrix(~0+group)
  y <- estimateDisp(y,design)
  
  list_df<-c(1,-1)
  
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,contrast =list_df)
  test1<-topTags(qlf,n = N)$table
  return(test1)
}

diffrk_analysis<-function(ADJ1,ADJ2){
  
  adjMtrxSample1=ADJ1
  adjMtrxSample2=ADJ2
  
  n_nodes=dim(adjMtrxSample1)[1]
  
  eps=10E-30
  lambda=0.5
  DBC=rep(1/n_nodes, n_nodes) 
  N=dim(adjMtrxSample1)[1];
  
  
  Delta_C_i_store=Matrix::Matrix((adjMtrxSample1-adjMtrxSample2),sparse=T)
  Delta_C_i_store@x=abs(Delta_C_i_store@x)
  Delta_C_i<-Matrix::Matrix(Delta_C_i_store/(10E-100+Matrix::rowSums(Delta_C_i_store)),sparse=T)
  
  
  error=100
  count=1
  
  solution=rep(1/N,N)
  
  while(error > eps)
  {
    count=count+1;
    formerSoulution=solution;
    
    s = as.vector(t(Delta_C_i)%*%solution)
    s[which(s==0)]<-sum(solution/N)
    solution = (1-lambda) * DBC + lambda * s
    
    error=sum((as.matrix(formerSoulution) - as.matrix(solution))^2)		
  }
  
  
  return(as.matrix(unlist(solution)))
}





for_figures<-function(data,select_number=1000,network_gamma=2.93,sparse=F){
  

gene_names<-row.names(data)
unique_cell_names<-unique(colnames(data))

print("CoreNet started")
test_edge<-ceiling((select_number)/(select_number^(1/(network_gamma-1))))+1
order_genes<-order(apply(data,1,var),decreasing=T)[1:select_number]
gene_names_test<-as.matrix(gene_names)[order_genes]

data<-data[order_genes,]

to_remove<-lapply(unique(colnames(data)),function(X){
  c(which(rowSums(data[,which(colnames(data)==X)])==0),
    which(is.na(apply(data[,which(colnames(data)==X)],1,var))==T)
  )
})
removalist<-unique(c(do.call('c',to_remove)))

print("Cosine Adjacency Matrix construction")
kernel<-lapply(c(1:length(unique_cell_names)),function(ce){
  if (sparse==F){
    S_cell<-as.matrix(proxy::dist((data[-removalist,which(colnames(data)==unique(colnames(data))[ce])]),method = "cosine"))
  }
  if (sparse==T){
    S_cell<-1/as.matrix(qlcMatrix::cosSparse(t(data[-removalist,which(colnames(data)==unique(colnames(data))[ce])])))
  }
  diag(S_cell)<-10E-100
  return((S_cell))
})

print("Cell-Type Network construction")
sub<-lapply(c(1:length(unique_cell_names)),function(ce){
  A<-kernel[[ce]]
  A=do.call('rbind',lapply(c(1:dim(A)[1]),function(X){A[X,]*(A[X,]<=(A[X,order(A[X,],decreasing = F)[test_edge]]))}))
  A_final=((A+t(A)))
  A_final<-t(apply(A_final,1,function(X){X/sum(X)}))
  return(A_final)
})


print("Core Network construction")
average_kernel<-Reduce('+',kernel)/length(kernel)

take<-Reduce('+',lapply(c(1:length(unique_cell_names)),function(X){
  D<-abs(kernel[[X]]-average_kernel)
  return(D)
}))

diag(take)<-0
A<-take
A=do.call('rbind',lapply(c(1:dim(A)[1]),function(X){A[X,]*(A[X,]<=(A[X,order(A[X,],decreasing = F)[test_edge]]))}))
core=A+t(A)

A<-average_kernel
A=do.call('rbind',lapply(c(1:dim(A)[1]),function(X){A[X,]*(A[X,]<=(A[X,order(A[X,],decreasing = F)[test_edge]]))}))
A_final=A+t(A)
core_for_diffrank<-t(apply(A_final,1,function(X){X/sum(X)}))


print("Graphlet counting")
core<-netdist::count_orbits_per_node(igraph::graph_from_adjacency_matrix(core!=0,"undirected",diag = F),4)   

cell_list<-lapply(c(1:length(unique_cell_names)),function(ce){
  return(netdist::count_orbits_per_node(igraph::graph_from_adjacency_matrix(sub[[ce]]!=0,"undirected",diag = F),4))
})

return(list(data=data,
            kernel=kernel,
            core_for_diffrank=core_for_diffrank,
            sub=sub,
            cell_list=cell_list,
            core=core,
            removalist=removalist,
            gene_names_test=gene_names_test,
            unique_cell_names=unique_cell_names))

}

#
#
#
#
#

data<-counts_10x_new
select_number=1000
corenet_for_figures<-for_figures(data,select_number=select_number,network_gamma=2.93,sparse=F)







#
#
#
#
#






print("Run DESeq")

run_deseq_core=1 # Turn on/off


if (run_deseq_core==1){
  
  edger_full_list_all<-as.vector(c(),mode = "list")
  deseq_list_all<-as.vector(c(),mode = "list")
  graphlet_list_all<-as.vector(c(),mode = "list")
  diffrank_list_all<-as.vector(c(),mode="list")
  
  par(mfcol=c(2,8))
  for (i in 1:8){
    chi_test<-(scale(log10((1+corenet_for_figures$cell_list[[i]])/(1+corenet_for_figures$core)))^2)
    ppvals<-pchisq(rowSums((-2*log(apply(chi_test,2,function(X){pchisq(q = X,df = 1,lower.tail = F)})))),2*dim(corenet_for_figures$core)[2],lower.tail = F)
    pvals<-p.adjust(ppvals,"BH")
    
    deseq_run<-run_deseq2(data = corenet_for_figures$data[-corenet_for_figures$removalist,],colnames(corenet_for_figures$data[-corenet_for_figures$removalist,]),N = select_number,i,c(1:8)[-i])
    edger_run<-run_edger(data = corenet_for_figures$data[-corenet_for_figures$removalist,],colnames(corenet_for_figures$data[-corenet_for_figures$removalist,]),N = select_number,i,c(1:8)[-i])
    diffrk_run<-diffrk_analysis(corenet_for_figures$core_for_diffrank,corenet_for_figures$sub[[i]])
    
    graphlet_fdr<-(-log10(10E-300+ppvals))
    deseq_fdr<-(-log10(10E-300+deseq_run$padj[match(row.names(deseq_run),corenet_for_figures$gene_names_test[-corenet_for_figures$removalist])]))
    edger_fdr<-(-log10(10E-300+edger_run$FDR[match(row.names(edger_run),corenet_for_figures$gene_names_test[-corenet_for_figures$removalist])]))
    edger_full_list_all<-c(edger_full_list_all,list(edger_run))
    deseq_list_all<-c(deseq_list_all,list(deseq_fdr))
    graphlet_list_all<-c(graphlet_list_all,list(-log10(10E-300+pvals)))
    diffrank_list_all<-c(diffrank_list_all,list(diffrk_run))
    
    id<-which(graphlet_fdr>(-log10(1)))
    
    smoothScatter(rank(deseq_fdr[id]),rank((graphlet_fdr[id])),main=unique(colnames(corenet_for_figures$data))[i],xlab="DESeq Model Rank",ylab="Graphlet Rank",colramp = colorRampPalette(c("white", blues9)),nbin = 300,bandwidth = length(id)/15,transformation = function(x) x^2)
    lm_mod<-lm(rank((graphlet_fdr[id]))~rank(deseq_fdr[id]))
    abline(lm_mod,col="red",lwd=2,lty=3)
    s_lm_mod<-summary(lm_mod)
    r<-s_lm_mod$r.squared
    r_pval<-pt(q = (r/sqrt(1-r^2))*sqrt(length(graphlet_fdr[id])-2),df =c(length(graphlet_fdr[id])-2),lower.tail = F)
    kendall_deseq<-cor.test(rank(deseq_fdr[id]),rank((graphlet_fdr)[id]),method = "kendall")
    legend("bottomright",legend = c(paste("Kendall: ",kendall_deseq$estimate),paste("p-value: ",kendall_deseq$p.value)),cex = 1.2)
    
    smoothScatter(rank(diffrk_run[id]),rank((graphlet_fdr)[id]),main=unique(colnames(corenet_for_figures$data))[i],xlab="DiffRank Model Rank",ylab="Graphlet Rank",colramp = colorRampPalette(c("white", blues9)),nbin = 300,bandwidth = length(id)/15,transformation = function(x) x^2)
    lm_mod<-lm(rank((graphlet_fdr)[id])~rank(diffrk_run[id]))
    abline(lm_mod,col="red",lwd=2,lty=3)
    s_lm_mod<-summary(lm_mod)
    r<-s_lm_mod$r.squared
    r_pval<-pt(q = (r/sqrt(1-r^2))*sqrt(length(graphlet_fdr[id])-2),df =c(length(graphlet_fdr[id])-2),lower.tail = F)
    kendall_diffrank<-cor.test(rank(diffrk_run[id]),rank((graphlet_fdr)[id]),method = "kendall")
    legend("bottomright",legend = c(paste("Kendall: ",round(kendall_diffrank$estimate,3)),paste("p-value: ",kendall_diffrank$p.value)),cex = 1.2)
  }
  
}



#
#
#
#
#





print("Run Barplot Venn diagram")

run_bar_venn=1 # Turn on/off


if (run_bar_venn==1){
  

lapply(c(1:8),function(X){
  
  graphlet_sig<-which(graphlet_list_all[[X]]>c(-log10(0.05)))
  edger_sig<-which(edger_list_all[[X]]>=c(-log10(0.05)))
  deseq_sig<-which(deseq_list_all[[X]]>=c(-log10(0.05)))
  diffrank_sig<-(order(diffrank_list_all[[X]],decreasing = T)[1:length(graphlet_sig)])
  
  
  listInput<-list("Graphlet"=graphlet_sig,"EdgeR"=edger_sig,"DESeq"=deseq_sig,"DiffRank"=diffrank_sig)
  UpSetR::upset(fromList(listInput),order.by = "freq", nsets = 4,  point.size = 3.5, line.size = 2,
                mainbar.y.label = "Gene Sets", sets.x.label = "Genes Per Method",
                text.scale = c(2, 2, 1.5, 2, 2, 2),empty.intersections = "on", sets.bar.color = "#56B4E9")
  
})

}




#
#
#
#
#




print("Run Circle Venn diagram")

run_circle_venn=1 # Turn on/off

if (run_circle_venn==1){
  
lapply(c(1:8),function(X){
  graphlet_sig<-which(graphlet_list_all[[X]]>=c(-log10(0.05)))
  edger_sig<-which(edger_list_all[[X]]>=c(-log10(0.05)))
  deseq_sig<-which(deseq_list_all[[X]]>=c(-log10(0.05)))
  diffrank_sig<-(order(diffrank_list_all[[X]],decreasing = T)[1:length(graphlet_sig)])
  
  venn<-VennDiagram::venn.diagram(list("EdgeR"=edger_sig,"DESeq"=deseq_sig,"CoreNet"=graphlet_sig,"DiffRank"=diffrank_sig),fill = 2:5, alpha = 0.3, filename = NULL)
  grid.draw(venn)
})
}


#
#
#
#
#


print("Run Network diagram")

run_network_cell=1 # Turn on/off

if (run_network_cell==1){
  
  truth=c(1:8)
  lapply(c(1:8),function(Z){
    if (any(grepl(Z,truth))){
      N=30
      id<-setdiff(which(graphlet_list_all[[Z]]>=(-log10(0.05))),unique(c(which(edger_list_all[[Z]]>=(-log10(0.05))),order(diffrank_list_all[[Z]],decreasing = T)[1:length(which(graphlet_list_all[[Z]]>=(-log10(0.05))))],which((deseq_list_all[[Z]])>=(-log10(0.05))))))
      id<-id[order(graphlet_list_all[[Z]][id],decreasing = T)[1:min(length(id),N)]]
      names<-(corenet_for_figures$gene_names_test[-corenet_for_figures$removalist])[id]
      fdr<-graphlet_list_all[[Z]][id]
      if(length(c(id))>1){
        print(id)
        top_N_genes<-apply(corenet_for_figures$kernel[[Z]][id,id],2,as.numeric)
        
        A<-apply(top_N_genes,2,function(X){
          internal<-rep(0,dim(top_N_genes)[1])
          internal[order(X,decreasing = T)[1:floor(dim(top_N_genes)/2)]]<-1
          return(internal)
        })
        top_N_genes<-A+t(A)
        row.names(top_N_genes)<-colnames(top_N_genes)<-names
        
        # fdr<-fdr[-which(colSums(top_N_genes)==0)]
        # names<-names[-which(colSums(top_N_genes)==0)]
        # top_N_genes<-top_N_genes[-which(colSums(top_N_genes)==0),-which(colSums(top_N_genes)==0)]
        
        
        
        graph <- graph_from_adjacency_matrix(top_N_genes,mode = "undirected",weighted = T,diag = F)
        ddf<-corenet_for_figures$data.frame(id=id,value=edger_full_list_all[[Z]]$logFC[id])
        rbPal <- colorRampPalette(c("light pink","light blue"))
        ddf$cols_f<-rbPal(2)[1+as.numeric(ddf$value>0)]
        print(ddf$cols)
        print(as.numeric(ddf$value))
        
        set.seed(1)
        ggraph(graph, layout = 'graphopt') + 
          geom_edge_link(aes(start_cap = label_rect(node1.name), end_cap = label_rect(node2.name),colour=E(graph)$weight)) + 
          theme_graph()+
          scale_edge_colour_gradient2( low = ("grey80"), mid = "grey80",
                                       high = ("grey80"), midpoint = mean(top_N_genes[top_N_genes>0]), space = "Lab",
                                       na.value = "grey50", guide = "edge_colourbar") + 
          geom_node_point(size = (fdr*2+5)^(1/1.3),color=sapply(ddf$cols_f,col2hex))+
          geom_node_text(aes(label = names)
                         , colour="black", fontface = "bold",cex=5)
      }
    }
  })
  
}




