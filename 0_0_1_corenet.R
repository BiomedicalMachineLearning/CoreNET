#' CoreNet v0.1.3 - 2019-July-18
#' @description CoreNet is a Network Tool to Identify Differential Gene Interactions.
#' It receives a set of networks (here, each network represents a different cell-type), in the form of an adjacency matrix 
#' and constructs a core network - a network that shares common edges across all networks in the set. The idea is to compare 
#' the network of a single cell-type (e.g. B Cell), with the common core network shared across all cell-types. The aim is to 
#' identify a set of differentially interacting genes that are contributors to the difference between the network of a single cell type, 
#' and the broader core network.  
#' 
#' Given that the networks aremade sparse according to a network property - the minimum degree for each node is defined by the 
#' power-law ( a scale-free network ).Computation is quick for a small number  of genes (around a thousand will take about a minute or 
#' under). However, given that networks scales with the size of the adjacency matrix several thousand genes may take several minutes 
#' to run or more. It should be noted that the full run-time will analyse all cell-types, identifying differential interacting genes
#' across all the unique cell-types, from where the core network is generated. The output is a list of genes, for each cell-type 
#' based network. Preferrably, a core network should be constructed from cell-types which have a common biological interpretation 
#' (e.g. immune subtypes, B Cell, Dendritic Cell, CD4 T helper Cell, Natural Killer Cell etc.)
#' 
#' The main component of statistical testing comes from a graphlet analysis, an evaluation of counts of sub-networks identified
#' in the larger, original cell-type (or core) network. Please refer to other documents for a comprehensive review of graphlets.
#' Once graphlets are counted for each gene, the counts for each graphlet type are transformed into Chi-square statistics by
#' first scaling it into a standard normal (provided there are sufficient number of genes (approximately > 1000) - normality assumptions 
#' should hold). 
#' 
#' @param data is a data frame with genes as rows and the cell types as columns. Names should be given for both rows and columns.
#' @param select_number the number of genes to analyse for differential gene interactions
#' @param network_gamma the network parameter to tweak the power exponent of the hypothetical network degree distribution. 
#' It is 1.93 in the literature (https://www.pnas.org/content/101/11/3765 # Universality and flexibility in gene expression from bacteria to human), 
#' but here the network notation follows the Barabasi book Network Science, so 1 is added to it
#' 
#' @examples 
#' diff_gene_interac<-corenet(data=data,select_number=3000,network_gamma=2.93)
#' @export
#' @author David N. Banh, Quan H. Nguyen 

###
###  ***
###

print("install packages before running")
# install.packages("devtools")
# devtools::install_github("alan-turing-institute/network-comparison")
# install.packages("qlcMatrix") 
# 
# library(netdist)  # uses count_orbits_per_node, a function to count orbits (permutation of graphlets)
# library(qlcMatrix)  # uses cosSparse, a function to calculate cosine distance

###
###  ***
###



corenet<-function(data,select_number=1000,network_gamma=2.93,sparse=F){
  gene_names<-row.names(data)
  unique_cell_names<-unique(colnames(data))
  
  print("CoreNet started")
  test_edge<-ceiling((select_number)/(select_number^(1/(network_gamma-1))))+1
  order_genes<-order(apply(data,1,var),decreasing=T)[1:select_number]
  gene_names_test<-as.matrix(gene_names)[order_genes]
  
  data_order<-data[order_genes,]
  
  to_remove<-lapply(unique(colnames(data_order)),function(X){
    c(which(rowSums(data_order[,which(colnames(data_order)==X)])==0),
      which(is.na(apply(data_order[,which(colnames(data_order)==X)],1,var))==T)
    )
  })
  removalist<-unique(c(do.call('c',to_remove)))
  
  print("Cosine Adjacency Matrix construction")
  kernel<-lapply(c(1:length(unique_cell_names)),function(ce){
    if (sparse==F){
      S_cell<-as.matrix(proxy::dist((data_order[-removalist,which(colnames(data_order)==unique(colnames(data_order))[ce])]),method = "cosine"))
    }
    if (sparse==T){
      S_cell<-1/as.matrix(qlcMatrix::cosSparse(t(data_order[-removalist,which(colnames(data_order)==unique(colnames(data_order))[ce])])))
    }
    diag(S_cell)<-10E-100
    return((S_cell))
  })
  
  print("Cell-Type Network construction")
  sub<-lapply(c(1:length(unique_cell_names)),function(ce){
    A<-kernel[[ce]]
    A=do.call('rbind',lapply(c(1:dim(A)[1]),function(X){A[X,]*(A[X,]<=(A[X,order(A[X,],decreasing = F)[test_edge]]))}))
    A_final=((A+t(A)))
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
  
  
  print("Graphlet counting")
  core<-netdist::count_orbits_per_node(igraph::graph_from_adjacency_matrix(core!=0,"undirected",diag = F),4)   

  cell_list<-lapply(c(1:length(unique_cell_names)),function(ce){
    return(netdist::count_orbits_per_node(igraph::graph_from_adjacency_matrix(sub[[ce]]!=0,"undirected",diag = F),4))
  })
  
  print("Graphlet analysis")
  graphlet_test<-lapply(c(1:length(unique_cell_names)),function(X){
    chi_test<-(scale(log10((1+cell_list[[X]])/(1+core)))^2)
    
    ppvals<-pchisq(rowSums((-2*log(apply(chi_test,2,function(X){pchisq(q = X,df = 1,lower.tail = F)})))),2*dim(core)[2],lower.tail = F)
    pvals<-p.adjust(ppvals,"BH")
    sigsbh=sort((gene_names_test[-removalist])[which(pvals<0.05)])
    sigsbf=sort((gene_names_test[-removalist])[which(ppvals<0.05/dim(core)[1])])
    
    return(list(p=pvals,genes_fdr=sigsbh,genes_p=sigsbf))
    
  })
  print("CoreNet finished")
  return(graphlet_test)
}

print("please install packages before running - see top, [CTRL-F ***]")

