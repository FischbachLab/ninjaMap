#!/usr/bin/env Rscript
library(tidyverse)
library(waffle)
library(hrbrthemes)
library(fs)
library(mdthemes)
library(scales)

total_samples = function(seedfile_dir, seedfile_ext,accepted_runs_list){
  # seedfile_dir = study_seedfiles_dir
  # seedfile_ext = "seedfile.csv"
  # accepted_runs_list = runs_list
  #   rm(seedfile_dir,seedfile_ext,accepted_runs_list)
  df = dir_ls(path = seedfile_dir, 
              glob = paste0('*.',seedfile_ext), 
              recurse = TRUE) %>% 
    map_df(read_csv,.id="filePath") %>% 
    rowwise() %>% 
    mutate(run_name=nth(str_split(basename(filePath),'\\.')[[1]], 1)) %>% 
    ungroup() %>% 
    filter(run_name %in% accepted_runs_list) %>% 
    distinct()
  return(nrow(df))
}

aggregate_abundance_data = function(parent_path) {
  abundance.df = dir_ls(path = parent_path, 
                        glob = "*.ninjaMap.abundance.csv", 
                        recurse = TRUE) %>% 
    map_df(read_csv,.id="filePath") %>% 
    rowwise() %>% 
    mutate(sample_id=nth(str_split(filePath,"/")[[1]], -3)) %>% 
    select(sample_id,Strain_Name,Read_Fraction,Percent_Coverage,Coverage_Depth) %>%
    ungroup()
}

aggregate_read_stats = function(parent_path){
  reads_acct.df = dir_ls(path = parent_path, 
                         glob = "*/Stats/read_accounting.csv", 
                         recurse = TRUE) %>% 
    map_df(read_csv,.id="filePath") %>% 
    dplyr::rename(sample_id = "Sample_Name") %>% 
    rowwise() %>%
    mutate(run_name = nth(str_split(filePath,"/")[[1]], -4)) %>%
    select(sample_id, run_name, Total_Fragments, Fragments_After_Trim, Fragments_Aligned) %>% 
    ungroup()
  
  reads_used.df = dir_ls(path = parent_path, 
                         glob = "*.ninjaMap.read_stats.csv", 
                         recurse = TRUE) %>% 
    map_df(read_csv,.id="filePath") %>% 
    rowwise() %>%
    mutate(sample_id = basename(File_Name),
           run_name = nth(str_split(filePath,"/")[[1]], -4)) %>% 
    select(-filePath,-File_Name) %>%
    ungroup()
  
  left_join(reads_acct.df,reads_used.df, by=c("sample_id"="sample_id", "run_name"="run_name"))
}

aggregate_fragment_stats = function(df,how="all"){
  switch(how, "all" = aggregate_fragment_stats_all(df), "post_qc" = aggregate_fragment_stats_post_qc(df))
}

aggregate_fragment_stats_all = function(read_stats_df){
  read_stats_df %>% 
    mutate(Singular_Fragments = Reads_wSingular_Votes/2,
           Escrowed_Fragments = Reads_wEscrowed_Votes/2,
           Discarded_Fragments = Discarded_Reads_w_Perfect_Aln/2,
           Unaligned_Fragments = (Fragments_After_Trim - Fragments_Aligned),
           QC_Fail_Fragments = (Total_Fragments-Fragments_After_Trim),
           All_Fragments=Total_Fragments) %>% 
    select(sample_id,contains('Fragments')) %>% 
    gather(Step,Frags, -c(sample_id,Total_Fragments)) %>% 
    rowwise() %>% 
    mutate(perc_frags = Frags*100/Total_Fragments) %>% 
    select(-Frags,-Total_Fragments) %>% 
    spread(Step,perc_frags)
}

aggregate_fragment_stats_post_qc = function(read_stats_df){
  read_stats_df %>% 
    mutate(Singular_Fragments = Reads_wSingular_Votes/2,
           Escrowed_Fragments = Reads_wEscrowed_Votes/2,
           Discarded_Fragments = Discarded_Reads_w_Perfect_Aln/2,
           Unaligned_Fragments = (Fragments_After_Trim - Fragments_Aligned),
           QC_Fail_Fragments = (Total_Fragments-Fragments_After_Trim),
           Fragments_After_Trim_tmp=Fragments_After_Trim) %>% 
    select(sample_id,contains('Fragments')) %>% 
    gather(Step,Frags, -c(sample_id,Fragments_After_Trim_tmp)) %>% 
    rowwise() %>% 
    mutate(perc_frags = Frags*100/Fragments_After_Trim_tmp) %>% 
    select(-Frags,-Fragments_After_Trim_tmp) %>% 
    spread(Step,perc_frags) %>% 
    select(-QC_Fail_Fragments, -Unaligned_Fragments)
}


download_from_s3 = function(local_path, s3_path,debug = FALSE){
  # check if local path exists; if not create it.
  if (! dir.exists(local_path)){
    dir.create(local_path, recursive = TRUE, mode = "0777")
  }
  
  debug_val = TRUE
  if (debug){
    debug_val = FALSE
  }
  
  system(paste("aws s3 cp",
               s3_path,
               local_path,
               "--recursive --exclude '*' --include '*/ninjaMap/*.ninjaMap.abundance.csv'", sep = " "),
         ignore.stdout = debug_val,
         ignore.stderr = debug_val)
  
  system(paste("aws s3 cp",
               s3_path,
               local_path,
               "--recursive --exclude '*' --include '*/Stats/read_accounting.csv'", sep = " "),
         ignore.stdout = debug_val,
         ignore.stderr = debug_val)
  
  system(paste("aws s3 cp",
               s3_path,
               local_path,
               "--recursive --exclude '*' --include '*/ninjaMap/*.ninjaMap.read_stats.csv'", sep = " "),
         ignore.stdout = debug_val,
         ignore.stderr = debug_val)
}

download_all_runs = function(list_runs,local_raw_output_path,s3_raw_output_path){
  for(run_name in runs_list){
    strains_local_path = paste(local_raw_output_path, run_name, sep = '/')
    strains_s3_path = paste(s3_raw_output_path, run_name, sep = '/')
    print(paste("Processing",strains_s3_path))
    download_from_s3(strains_local_path, strains_s3_path)
  }
}

make_dir = function(dir_path){
  if (! dir.exists(dir_path)){
    dir.create(dir_path, recursive = TRUE, mode = "0777")
  }
  return(dir_path)
}

## PLOTS ##
brian_plot = function(df,title="Brian Plot - Basic"){
  # requires 3 columns: sample_id, Strain_Name (of the strain), Read_Fraction
  p = df %>%
    ungroup() %>% 
    mutate(label=if_else(sample_id == Strain_Name, 15, 0),
           label=na_if(label,0),
           sample_id = str_replace_all(sample_id, "-", " "),
           abundance = Read_Fraction/100,
           # sample_id = fct_relevel(sample_id,sort(unique(sample_id)))
           sample_id = formatted_name
           ) %>% 
    # group_by(sample_id) %>%
    # mutate(norm_abundance = abundance/sum(abundance)) %>%
    # ungroup() %>%
    # ggplot(aes(sample_id,norm_abundance))+
    ggplot(aes(sample_id,abundance))+
    # geom_hline(yintercept = 1e-2, linetype="dashed", color = "red")+
    # geom_hline(yintercept = 1, linetype="dashed", color = "red")+
    # geom_violin(alpha = 0.5,draw_quantiles = TRUE,trim = TRUE, color = "lightgrey")+
    # geom_beeswarm(aes(color=Strain_Name),
    #               size = 1,
    #               # varwidth = TRUE, 
    #               show.legend = FALSE)+
    # geom_point(aes(color=Strain_Name), show.legend = FALSE)+
    geom_point(color="grey", show.legend = FALSE, shape=21, stroke=0.75)+
    geom_point(data = (. %>% filter(!is.na(label))),
               # color = 'black', shape=21, size=2)+
               color = 'black', fill='#66c2a5', shape=21, size=3, stroke=1)+
    # geom_line(aes(group=Strain_Name, color=Strain_Name), show.legend = FALSE, alpha = 0.2)+
    scale_y_continuous(trans='log10', breaks = 10^(seq(-10,10,1)),
                       labels = trans_format("log10", math_format(10^.x)))+
    # scale_x_discrete(aes(labels=formatted_name))+
    # coord_cartesian(ylim=c(1e-5, 100))+
    # theme_bw()+
    # theme(text = element_text(size = 13),
    #       axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10))+
    xlab("Query")+
    ylab("Relative abundance")+
    ggtitle(title)+
    md_theme_bw()
  
  
  p$theme$axis.text.x$angle=90
  p$theme$axis.text.x$hjust=1
  p$theme$axis.text.x$vjust=0.5
  p$theme$axis.text.x$size=10
  p$theme$axis.text.y$size=15
  p$theme$text$size=15
  # p1$theme$axis.title.x="Did this work?"
  p$theme$text$size=15
  
  return(p)
}

plot_fragment_stats = function(fragment_stats_df){
  library(ggbeeswarm)
  si_read_stats_df_fig = fragment_stats_df %>% 
    rowwise() %>% 
    # mutate(Total_Fragments_Lost = sum(QC_Fail_Fragments,Unaligned_Fragments,Discarded_Fragments, na.rm = TRUE)) %>% 
    gather(Step,perc_of_total,-sample_id) # %>% 
  # filter(! Step %in% c('QC_Fail_Fragments','Unaligned_Fragments','Discarded_Fragments'))
  
  si_read_stats_df_fig$Step = factor(si_read_stats_df_fig$Step, 
                                     levels = c('All_Fragments', 
                                                    'QC_Fail_Fragments',
                                                    'Fragments_After_Trim',
                                                    'Fragments_Aligned','Unaligned_Fragments',
                                                    'Singular_Fragments','Escrowed_Fragments',
                                                    'Discarded_Fragments'))
  
  si_read_stats_df_fig %>% 
    filter(! is.na(Step)) %>% 
    ggplot(aes(Step,perc_of_total,fill=Step)) +
    # geom_violin(aes(fill=Step))+
    geom_line(aes(group=sample_id), alpha=0.1)+
    geom_boxplot(aes(fill=Step))+
    scale_fill_viridis_d(option="C")+
    # geom_beeswarm(aes(color=Step))+
    geom_rug(aes(color=Step),sides = "l")+
    scale_color_viridis_d(option="C")+
    scale_y_continuous(breaks = seq(0,100,5))+
    ggtitle("Study Read Fate")+
    # coord_flip()+
    theme_bw()+
    theme(text = element_text(size = 10),
          axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5))
}

#' Generate a waffle plot for the given data frame.
#'
#' @description
#' The waffle_plot function expects a data frame with at least 2 columns named `names` and `values`,
#' returns a proportional 10 x 10 grid waffle plot ggplot object
#' 
#' @param df - a data frame. 
#' @param title - what are we plotting today?
#' @param units - optional, default='%'.
#' @param subtitle - optional.
#'
#' @return ggplot object.
waffle_plot = function(df, title, units = "%", subtitle=NULL){
  df %>% 
    ungroup() %>% 
    # arrange(desc(values)) %>% 
    # mutate(names = fct_relevel(names, names)) %>% 
    mutate(print_name = paste0(names, " (", round(values,2)," ",units,")")) %>% 
    ggplot(aes(fill=print_name, values=values)) +
    geom_waffle(color = "white", size=1, n_rows = 10, make_proportional = TRUE)+
    # facet_wrap(~fct, ncol=1) 
    scale_x_discrete(expand=c(0,0))+
    scale_y_discrete(expand=c(0,0)) +
    ggthemes::scale_fill_tableau(name=NULL) +
    coord_equal() +
    labs(
      title = title, subtitle = subtitle
    ) +
    theme_ipsum_rc(grid="") +
    theme_enhance_waffle()
}

#' Summarize fragment stats for each step
#'
#' @param frag_stats_df 
#'
#' @return
summarize_stats_by_steps = function(frag_stats_df) {
  
  frag_stats_df %>% 
    group_by(Step) %>%
    summarise_if(is.numeric, median, na.rm = TRUE) %>% 
    mutate(names = case_when(Step == "QC_Fail_Fragments" ~ "QC_Fail",
                             Step == "Unaligned_Fragments" ~ "Unaligned",
                             Step == "Singular_Fragments" ~ "Primary",
                             Step == "Escrowed_Fragments" ~ "Escrow",
                             Step == "Discarded_Fragments" ~ "Discard",
                             Step == "Aligned_Missed_Fragments" ~ "Missed")) %>% 
    filter(!is.na(names)) %>% 
    rename(values = perc_frags) %>% 
    select(names, values)
}

process_fragment_stats = function(df){
  df %>%
    mutate(Singular_Fragments = Reads_wSingular_Votes/2,
           Escrowed_Fragments = Reads_wEscrowed_Votes/2,
           Discarded_Fragments = Discarded_Reads_w_Perfect_Aln/2,
           Unaligned_Fragments = (Fragments_After_Trim - Fragments_Aligned),
           Aligned_Missed_Fragments = (Fragments_Aligned - Singular_Fragments - Escrowed_Fragments - Discarded_Fragments),
           QC_Fail_Fragments = (Total_Fragments-Fragments_After_Trim),
           All_Fragments=Total_Fragments) %>% 
    select(sample_id,contains('Fragments')) %>% 
    gather(Step,Frags, -c(sample_id,Total_Fragments)) %>% 
    rowwise() %>% 
    mutate(perc_frags = Frags*100/Total_Fragments) %>% 
    select(-Frags,-Total_Fragments)
}