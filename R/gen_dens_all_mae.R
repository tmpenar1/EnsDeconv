#' This function is used to generate density plot of pearson correlation and MAE
#' @export
gen_dens_all_mae = function(re_org_out_all,tit,customed_mrk = FALSE, cust,notation = T, note = NULL){
  # library(ggridges)
  # library(ggplot2)
  # library(ggpubr)
  
  ########## customed marker name need to be changed
  if(customed_mrk){
    re_org_out_all = re_org_out_all %>% mutate(Marker.Method = ifelse(Marker.Method == "customed",cust,Marker.Method))
  }
  
  re_sub_pe = re_org_out_all %>% ungroup() %>%
    mutate(NormalizationTC = ifelse(str_detect(Method,"fold|Mean|allgene") ,"none",NormalizationTC),
           Scale =ifelse(str_detect(Method,"fold|Mean|allgene"),"none",Scale) )%>%
    select(Method,NormalizationTC,Pearson,CellType,Scale) %>% unique()
  
  re_sub_MAE = re_org_out_all %>% ungroup() %>%
    mutate(NormalizationTC = ifelse(str_detect(Method,"fold|Mean|allgene"),"none",NormalizationTC),
           Scale =ifelse(str_detect(Method,"fold|Mean|allgene"),"none",Scale) ) %>%
    select(Method,NormalizationTC,MAE_Total,Scale) %>% unique()
  
  re_sub_pe_total = re_org_out_all %>% ungroup() %>%
    mutate(NormalizationTC = ifelse(str_detect(Method,"fold|Mean|allgene"),"none",NormalizationTC),
           Scale =ifelse(str_detect(Method,"fold|Mean|allgene"),"none",Scale) ) %>%
    select(Method,NormalizationTC,Pearson_Total,Scale) %>% unique()
  
  cuts <-  c("2fold_nnls","2fold_solnp","2fold_mean","1fold_lar_nnls","1fold_lar_solnp","1fold_lar_mean","2foldrm_nnls","2foldrm_solnp","2foldrm_mean","1foldrm_lar_nnls","1foldrm_lar_solnp","1foldrm_lar_mean","allgenerm_mean","allgene_mean","Mean")
  re_sub_pe1 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[1])
  
  re_sub_pe2 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[2])
  
  re_sub_pe3 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[3])
  
  re_sub_pe4 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[4])
  
  re_sub_pe5 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[5])
  
  re_sub_pe6 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[6])
  
  re_sub_pe7 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[7])
  
  re_sub_pe8 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[8])
  
  re_sub_pe9 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[9])
  
  
  re_sub_pe10 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[10])
  
  re_sub_pe11 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[11])
  
  re_sub_pe12 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[12])
  
  re_sub_pe13 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[13])
  
  re_sub_pe14 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[14])
  
  re_sub_pe15 = re_sub_pe %>%
    group_by(CellType) %>%
    filter(Method==cuts[15])
  
  
  col = get_palette('jco',15)
  
  re_sub_pe_print = re_sub_pe %>% filter(Method %in% cuts) %>% unique()%>%
    select(-Scale,-NormalizationTC)
  re_sub_pe_print <- as.matrix(spread(re_sub_pe_print, Method, Pearson))
  
  rownames(re_sub_pe_print) = re_sub_pe_print[,1]
  re_sub_pe_print = re_sub_pe_print[,cuts]
  colnames(re_sub_pe_print) = toupper(gsub("ensemble_","",gsub("fold","",colnames(re_sub_pe_print))))
  
  tbl1 = re_sub_pe_print
  
  re_sub_pe = re_sub_pe %>% group_by(CellType)%>% mutate(Pearson = ifelse(is.na(Pearson),min(Pearson, na.rm = T),Pearson))
  p = re_sub_pe%>%
    ggplot(aes(x=Pearson,y=CellType)) + # , fill=CellType
    stat_density_ridges(quantile_lines = TRUE,scale = 0.8,alpha = .5,)+
    geom_segment(data = re_sub_pe1, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[1])+
    geom_segment(data = re_sub_pe2, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[2])+
    geom_segment(data = re_sub_pe3, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[3])+
    geom_segment(data = re_sub_pe4, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[4])+
    geom_segment(data = re_sub_pe5, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[5])+
    geom_segment(data = re_sub_pe6, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[6])+
    geom_segment(data = re_sub_pe7, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[7])+
    geom_segment(data = re_sub_pe8, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[8])+
    geom_segment(data = re_sub_pe9, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                        yend = as.numeric(CellType) + .8), size=1, color = col[9])+
    geom_segment(data = re_sub_pe10, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[10])+
    geom_segment(data = re_sub_pe11, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[11])+
    geom_segment(data = re_sub_pe12, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[12])+
    geom_segment(data = re_sub_pe13, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[13])+
    geom_segment(data = re_sub_pe14, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[14])+
    geom_segment(data = re_sub_pe15, aes(x = Pearson, xend = Pearson,y = as.numeric(CellType),
                                         yend = as.numeric(CellType) + .8), size=1, color = col[15])+
    labs(x= "Pearson's correlation" ,title = tit)+
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(-0.5, 1))+
    theme_ridges(grid = FALSE, center = TRUE) + theme(legend.position = "none")+
    scale_colour_manual(name="Ensemble",values=c('2_NNLS'=col[1],'2_SOLNP'=col[2],'2_MEAN'=col[3], '1_LAR_NNLS' = col[4],'1_LAR_SOLNP' = col[5],"1_LAR_MEAN" = col[6],"2RM_NNLS" = col[7],"2RM_SOLNP" = col[8],"2RM_MEAN" = col[9],"1RM_LAR_NNLS" = col[10],"1RM_LAR_SOLNP" = col[11],"1RM_LAR_MEAN" = col[12],"ALLGENERM_MEAN" = col[13],"ALLGENE_MEAN" = col[14],"MEAN" = col[15]))
  
  
  re_sub_MAE2 = re_sub_MAE %>%
    filter(Method%in% cuts)%>%
    arrange(factor(Method, levels = cuts))
  
  print_MAE = re_sub_MAE2 %>%
    mutate(Method = gsub("ensemble_","",Method),
           Method = gsub("fold","",Method),
           Method = toupper(Method)) %>% select(Method,MAE_Total)
  
  print_MAE = t(print_MAE)
  
  
  re_sub_MAE$y = ''
  p_MAE =  re_sub_MAE%>%
    ggplot(aes(x=MAE_Total, y = y, fill = factor(stat(quantile)))) + ylab('Density') +
    # stat_density_ridges(quantile_lines = TRUE,scale = 0.8,alpha = .5)+ # , color="#e9ecef"
    # scale_fill_viridis_d(name = "Quartiles")+
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = c(0.25,0.5, 0.75),alpha = .5
    ) +
    scale_fill_manual(
      name = "Probability", values = c("#D3D3D3",  "#9E9E9E", "#7D7D7D","#696969"),
      labels = c("(0, 0.25]","(0.25,0.50]", "(0.50, 0.75]", "(0.75, 1]")
    )+
    theme_ridges(grid = FALSE, center = TRUE)+
    # theme_bw()+
    # geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[1],color = "2_NNLS"), size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[1], y=2, label=re_sub_MAE2$Method[1],angle = 90,color = col[1], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[2],color = "2_SOLNP"), size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[2], y=2, label=re_sub_MAE2$Method[2],angle = 90,color = col[2], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[3],color = "2_MEAN"),   size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[4], y=2, label=re_sub_MAE2$Method[4],angle = 90,color = col[4], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[4],color = "1_LAR_NNLS"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[5],color = "1_LAR_SOLNP"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[6],color = "1_LAR_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[7],color = "2RM_NNLS"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[8],color = "2RM_SOLNP"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[9],color = "2RM_MEAN"), size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[1], y=2, label=re_sub_MAE2$Method[1],angle = 90,color = col[1], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[10],color = "1RM_LAR_NNLS"), size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[2], y=2, label=re_sub_MAE2$Method[2],angle = 90,color = col[2], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[11],color = "1RM_LAR_SOLNP"),   size = 1) +
    # geom_text(x=re_sub_MAE2$MAE_Total[4], y=2, label=re_sub_MAE2$Method[4],angle = 90,color = col[4], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[12],color = "1RM_LAR_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[13],color = "ALLGENERM_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[14],color = "ALLGENE_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_MAE2$MAE_Total[15],color = "MEAN"),   size = 1) +
    scale_colour_manual(name="Ensemble",values=c('2_NNLS'=col[1],'2_SOLNP'=col[2],'2_MEAN'=col[3], '1_LAR_NNLS' = col[4],'1_LAR_SOLNP' = col[5],"1_LAR_MEAN" = col[6],"2RM_NNLS" = col[7],"2RM_SOLNP" = col[8],"2RM_MEAN" = col[9],"1RM_LAR_NNLS" = col[10],"1RM_LAR_SOLNP" = col[11],"1RM_LAR_MEAN" = col[12],"ALLGENERM_MEAN" = col[13],"ALLGENE_MEAN" = col[14],"MEAN" = col[15]))
  # geom_text(x=re_sub_MAE2$MAE_Total[5], y=2, label=re_sub_MAE2$Method[5],angle = 90,color = col[5], vjust = -0.2)+
  # geom_vline(data = d2, aes(xintercept = lower, linetype="dashed")) +
  # geom_vline(data = d2, aes(xintercept = med, linetype="dashed"))+
  # geom_vline(data = d2, aes(xintercept = upper, linetype="dashed"))+
  
  re_sub_pe_total2 = re_sub_pe_total %>%
    filter(Method%in% cuts)%>%
    arrange(factor(Method, levels = cuts))
  
  print_pe_total = re_sub_pe_total2 %>%
    mutate(Method = gsub("ensemble_","",Method),
           Method = gsub("fold","",Method),
           Method = toupper(Method)) %>% select(Method,Pearson_Total)
  
  print_pe_total = t(print_pe_total)
  tbl = as.data.frame(rbind(print_MAE,print_pe_total))
  colnames(tbl) = tbl[1,]
  tbl = tbl[-c(1,3),]
  #tbl <- tableGrob(tbl)
  re_sub_pe_total$y = ''
  p_pet =  re_sub_pe_total%>%
    ggplot(aes(x=Pearson_Total, y = y, fill = factor(stat(quantile)))) + ylab('Density') +
    # stat_density_ridges(quantile_lines = TRUE,scale = 0.8,alpha = .5)+ # , color="#e9ecef"
    # scale_fill_viridis_d(name = "Quartiles")+
    stat_density_ridges(
      geom = "density_ridges_gradient",
      calc_ecdf = TRUE,
      quantiles = c(0.25,0.5, 0.75),alpha = .5
    ) +
    scale_fill_manual(
      name = "Probability", values = c("#D3D3D3",  "#9E9E9E", "#7D7D7D","#696969")
    )+
    theme_ridges(grid = FALSE, center = TRUE)+
    # theme_bw()+
    # geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[1],color = "2_NNLS"), size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[1], y=2, label=re_sub_pe_total2$Method[1],angle = 90,color = col[1], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[2],color = "2_SOLNP"), size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[2], y=2, label=re_sub_pe_total2$Method[2],angle = 90,color = col[2], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[3],color = "2_MEAN"),   size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[4], y=2, label=re_sub_pe_total2$Method[4],angle = 90,color = col[4], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[4],color = "1_LAR_NNLS"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[5],color = "1_LAR_SOLNP"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[6],color = "1_LAR_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[7],color = "2RM_NNLS"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[8],color = "2RM_SOLNP"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[9],color = "2RM_MEAN"), size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[1], y=2, label=re_sub_pe_total2$Method[1],angle = 90,color = col[1], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[10],color = "1RM_LAR_NNLS"), size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[2], y=2, label=re_sub_pe_total2$Method[2],angle = 90,color = col[2], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[11],color = "1RM_LAR_SOLNP"),   size = 1) +
    # geom_text(x=re_sub_pe_total2$Pearson_Total[4], y=2, label=re_sub_pe_total2$Method[4],angle = 90,color = col[4], vjust = -0.2)+
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[12],color ="1RM_LAR_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[13],color = "ALLGENERM_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[14],color = "ALLGENE_MEAN"),   size = 1) +
    geom_vline(aes(xintercept = re_sub_pe_total2$Pearson_Total[15],color = "MEAN"),   size = 1) +
    scale_colour_manual(name="Ensemble",values=c('2_NNLS'=col[1],'2_SOLNP'=col[2],'2_MEAN'=col[3], '1_LAR_NNLS' = col[4],'1_LAR_SOLNP' = col[5],"1_LAR_MEAN" = col[6],"2RM_NNLS" = col[7],"2RM_SOLNP" = col[8],"2RM_MEAN" = col[9],"1RM_LAR_NNLS" = col[10],"1RM_LAR_SOLNP" = col[11],"1RM_LAR_MEAN" = col[12],"ALLGENERM_MEAN" = col[13],"ALLGENE_MEAN" = col[14],"MEAN" = col[15]))+
    # geom_text(x=re_sub_pe_total2$Pearson_Total[5], y=2, label=re_sub_pe_total2$Method[5],angle = 90,color = col[5], vjust = -0.2)+
    # geom_vline(data = d2, aes(xintercept = lower, linetype="dashed")) +
    # geom_vline(data = d2, aes(xintercept = med, linetype="dashed"))+
    # geom_vline(data = d2, aes(xintercept = upper, linetype="dashed"))+
    theme(legend.position = "none")
  
  
  figure3 = ggarrange(p,tableGrob(tbl1), p_MAE,p_pet,tableGrob(tbl),ncol = 1, heights = c(2,2, 2,2,0.5))
  
  # p1 = annotate_figure(figure3, top = text_grob(tit, color = "red", face = "bold", size = 14))
  
  if(notation){
    figure3 = annotate_figure(figure3,
                              bottom = text_grob(note, color = "blue",
                                                 hjust = 1, x = 1, face = "italic", size = 10)
    )
  }
  
  return(list(p = figure3,tbl1 = tbl1,tbl2 = tbl))
}
