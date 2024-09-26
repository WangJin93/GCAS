load("www/genelist.rda")
themes_list <- list(
  "Base" = ggthemes::theme_few(),
  "cowplot" = cowplot::theme_cowplot(),
  "Light" = theme_light(),
  "Minimal" = theme_minimal(),
  "Classic" = theme_classic(),
  "Gray" = theme_gray(),
  "half_open" = cowplot::theme_half_open(),
  "minimal_grid" = cowplot::theme_minimal_grid()
)


stage_plot<-function(clinic_data,feature="Tumor_Stage_simplify",is.cont = F,cont.cut){
  clinic <- clinic_data %>% dplyr::select(c("ID",feature,colnames(clinic_data)[3]))
  clinic[clinic == ""] = NA
  clinic <- na.omit(clinic)
  if (is.cont){
    clinic[,feature] <- ifelse(clinic[,feature] > cont.cut, paste0("> ",cont.cut),paste0("≤ ",cont.cut))
  }

  clinic[,feature] <- factor(clinic[,feature])
  colnames(clinic)[2] <- "Group"
  group.n <- unique(clinic$Group)
  p<-ggplot(clinic,
            aes(x=Group,y=get(colnames(clinic_data)[3]),fill=Group))+
    geom_bar(stat="summary",fun=mean,position="dodge",color="black")+
    stat_summary(fun.data = 'mean_sd', geom = "errorbar", linewidth = 0.4,size=0.8,width=0.3,position = position_dodge(0.9),color="black")+
    labs(x=NULL,y='Relative protein expression')
  # if  ("None" %in% my.compair){
  #   p <- p
  # }else{
  #   if (group.n>2){
  #     all_compare=c()
  #     groups <- unique(clinic$Group)
  #     for (i in 1:(length(groups)-1)) {
  #       for (j in (i+1):length(groups)) {
  #         all_compare <- c(all_compare, paste(groups[i],"vs",groups[j]))
  #       }
  #     }
  #     my.compair <- strsplit(all_compare," vs ")
  #     p<-p+
  #       stat_compare_means(method = 't.test',label = "p.signif", comparisons = my.compair,size=5)
  #
  #   }else if (group.n==2){
  #     p<-p+ stat_compare_means(method = 't.test',label.y = NULL,label = "p.signif",comparisons = "all",size=5)
  #   }
  #
  # }

  df_p<-list("df"=clinic,"p"=p)
  return(df_p)
}



