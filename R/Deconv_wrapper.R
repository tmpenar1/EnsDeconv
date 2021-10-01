#' This function is a wrapper for ensemble deconvolution
#'
#' @param allfold_res all gene deconvolution results
#' @param allgene_res two folds deconvolution results
#' @param df_true Optional true cell type fraction
#' @param trueMet character 
#' 
#'
#' @export

deconv_wrap <- function(count_bulk,meta_bulk = NULL,ref_list,customed_markers = NULL,markers_range = NULL,true_frac = NULL,params = NULL,ncv_input =2,
                        outpath,parallel_comp = FALSE,ncore,os = "win",myseed = 2020,feature.select = FALSE,clustering = TRUE,rm.duplicated =FALSE,mrkpen = FALSE,dmeths = NULL,cur_mrk,trueMet){

  res = list()
  # Windows user
  for (i  in 1:2) {
    res[[i]] = gen_all_res_list(count_bulk = as.matrix(count_bulk), meta_bulk = meta_bulk, ref_list = ref_list, true_frac =true_frac, ncv_input = i, outpath =outpath, ncore =ncore, parallel_comp = parallel_comp, params = params)
  }
  
  res[[1]] = lapply(res[[1]], function(x){
    gene =  intersect(rownames(count_bulk),unique(unlist(cur_mrk)))
    x[[3]][[1]] = NULL
    x[[3]][[1]] = list()
    x[[3]][[1]][["test_bulk"]] = as.matrix(count_bulk)[gene,]
    return(x)
  })
  
  res[[2]] = lapply(res[[2]], function(x){
    gene = intersect(rownames(x[[3]][[1]][[1]]),unique(unlist(cur_mrk)))
    x[[3]][[1]][[1]] = x[[3]][[1]][[1]][gene,]
    
    gene = intersect(rownames(x[[3]][[2]][[1]]),unique(unlist(cur_mrk)))
    x[[3]][[2]][[1]] = x[[3]][[2]][[1]][gene,]
    return(x)
  })
  
  # Check available scenarios
  ind = sapply(res[[1]], length)
  ind2 = sapply(res[[2]], length)
  ind3 = union(ind,ind2)
  res[[1]] = res[[1]][which(ind3 == 3)]
  res[[2]] = res[[2]][which(ind3 == 3)]

  out_res = get_res_wrap(allfold_res = res[[2]],allgene_res = res[[1]],df_true = true_frac,trueMet = trueMet)
  
  
  return(list(res = res,out_res=out_res))
}
