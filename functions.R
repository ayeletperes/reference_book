## functions

hline <- function(y = 0, color = "red") {
  list(
    type = "line",
    x0 = 0,
    x1 = 1,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, dash = "dot")
  )
}

plot_zygousity <- function(tmp, state, allele_thresh, g){
  
  tmp_plot <- dplyr::filter(tmp, zygousity_state == state) %>% dplyr::rowwise() %>% dplyr::mutate(
    v_alleles_p = v_alleles_abc,
    v_alleles_p = gsub(";", "\n", v_alleles_p),
    text = paste(
      '</br>Project: ',
      project,
      '</br>Subject: ',
      subject,
      '</br>Alleles: ',
      v_alleles,
      '</br># assignments: ',
      count,
      '</br>Relative freq.: ',
      round(freq, 4),
      '</br>Relative Rep. freq.: ',
      round(freq2, 4)
    )
  ) %>% ungroup()
  
  loc2 <-
    setNames(1:length(unique(tmp_plot$v_allele_axis)), sort(unique(tmp_plot$v_allele_axis)))
  
  tmp_plot$loc2 <- loc2[tmp_plot$v_allele_axis]
  
  if(state!=1 & length(unique(tmp_plot$v_alleles_p))!=1){
  loc_jitter <- c()
  for(ii in unique(tmp_plot$loc2)){
    loc_jitter[[ii]] <-
      seq(0, 0.5, length.out = length(unique(tmp_plot$v_alleles_p[tmp_plot$loc2==ii])))
    
    loc_jitter[[ii]]  <-
      setNames(loc_jitter[[ii]] , sort(unique(tmp_plot$v_alleles_p[tmp_plot$loc2==ii])))
  }
  
  
  tmp_plot <-
    tmp_plot %>% dplyr::arrange(loc2, v_alleles_p) %>% dplyr::group_by(loc2) %>% 
   dplyr:: mutate(loc_plot = loc2 + loc_jitter[[unique(loc2)]][v_alleles_p],) %>% ungroup() %>% 
    dplyr::mutate(jitter_offset = (loc_plot))
  }else{
    tmp_plot <-
      tmp_plot %>% dplyr::arrange(loc2, v_alleles_p) %>% dplyr::group_by(loc2) %>% 
      dplyr:: mutate(loc_plot = loc2,) %>% ungroup() %>% 
      dplyr::mutate(jitter_offset = (loc_plot))
  }
  
  tickvals_tmp <-
    tmp_plot %>% dplyr::pull(loc_plot) %>% unique() %>% sort()
  
  tickvals <- c()
  
  for (i in 1:length(loc2)) {
    tickvals <- c(tickvals, mean(tickvals_tmp[floor(tickvals_tmp) == i]))
  }
  
  
  ticktext <-
    tmp_plot %>% dplyr::pull(v_allele_axis) %>% unique() %>% sort()
  
  plotly1 <-
    tmp_plot %>% rowwise() %>% dplyr::mutate(group = paste0(project, "-", v_alleles_p)) %>%
    highlight_key(., ~ subject) %>%
    plot_ly() %>%
    add_trace(
      type = "scatter",
      x = ~ (jitter_offset),
      y = ~ freq,
      text = ~ text,
      symbol = ~ project,
      mode = 'markers',
      marker = list(color = "grey", size = 12),
      showlegend = TRUE,
      opacity = 0.9,
      hoverinfo = 'none',
      legendgroup = ~ project
    ) %>%
    add_trace(
      type = "scatter",
      x = ~ (jitter_offset),
      y = ~ freq,
      text = ~ text,
      color = ~ v_alleles_p,
      mode = 'markers',
      showlegend = FALSE,
      opacity = 0.8,
      hoverinfo = 'text',
      legendgroup = ~ v_alleles_p
    ) %>%
    add_trace(
      x = ~ as.numeric(loc_plot),
      y = ~ freq,
      color = ~ v_alleles_p,
      type = "box",
      hoverinfo = "none",
      fillcolor = "transparent",
      name = ~ v_alleles_p,
      legendgroup = ~ v_alleles_p
    ) %>%
    layout(
      hovermode = 'closest',
      shapes = list(hline(allele_thresh/100)),
      legend = list(
        tracegroupgap = 20,
        title = list(text =
                       '<b>  </b>'),
        orientation = "V"
      ),
      xaxis = list(
        title = paste0(g, " Alleles"),
        autotick = F,
        tickmode = "array",
        tickvals = tickvals,
        ticktext = ticktext
      ),
      yaxis = list(title = "Relative allele frequency",
                   range = c(0,1.05))
    )  %>% plotly::highlight(
      on = "plotly_click",
      selected = attrs_selected(showlegend = FALSE),
      opacityDim = 0.3,
      persistent = TRUE
    ) %>% plotly_build()
  
  return(plotly1)
}

data_cutoff <- function(tmp, func_groups, g, allele_thresh = 0.5, or_allele){
  tmp <- tmp %>%
    dplyr:: filter(v_gene == func_groups[as.character(g)], !is.na(v_allele)) %>% 
    ungroup()
  
  tmp <- tmp %>% dplyr::group_by(subject)
  
  tmp <- tmp %>% dplyr::arrange(desc(freq)) %>%
    dplyr::group_by(subject, v_gene) %>% dplyr::mutate(
      zygousity_state = as.numeric(sum(freq > allele_thresh/100, na.rm = T)),
      v_alleles = paste0(1:unique(zygousity_state), " - ", or_allele[v_allele[1:unique(zygousity_state)]], collapse = ";"),
      v_alleles_abc = paste0(sort(or_allele[v_allele[1:unique(zygousity_state)]]), collapse = ";"),
      v_allele_axis = or_allele[v_allele]
    ) %>% arrange(subject)
  tmp <- tmp %>% dplyr::group_by(subject, zygousity_state) %>% dplyr::mutate(loc_state = loc <= zygousity_state) %>% filter(loc_state) %>% ungroup()
  
  return(tmp)
}


seq_align <- function(v_calls, allele_db, vgerms, chain, mat, g_group){
  alleles <- allele_db %>% dplyr::filter(new_allele %in% v_calls) %>% dplyr::pull(or_allele)
  new_alleles <- setNames(allele_db %>% filter(new_allele %in% v_calls) %>% dplyr::pull(new_allele), alleles)
  sequences <- substr(vgerms[[chain]][alleles],1,318)
  names(sequences) <- new_alleles[names(sequences)]
  
  mat_sub <- mat[alleles,alleles]
  
  colnames(mat_sub) <-  gsub(paste0(g_group,"[*]"),"",new_alleles[colnames(mat_sub)])
  rownames(mat_sub) <-  gsub(paste0(g_group,"[*]"),"",new_alleles[rownames(mat_sub)])
  
  hc <- hclust(as.dist(mat_sub), method = "complete")
  dend <- as.dendrogram(hc)
  dend <- dendextend::set(dend, "labels_cex", 1.2)
  ggd1 <- as.ggdend(dend)
  p_dend <- ggplot(ggd1, theme = bbplot::bbc_style())  + theme(
    axis.line = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.border = element_blank(), panel.background = element_blank(), 
    legend.position = "none" )
  
  
  matrix_sequences <- as.data.frame(sapply(sequences,seqinr::s2c), stringsAsFactors = F)
  matrix_sequences$annot <- apply(matrix_sequences, 1, function(x) length(unique(x)) != 1) 
  matrix_sequences$pos <- 1:318
  matrix_sequences_plot <- reshape2::melt(matrix_sequences, id.vars = c("pos","annot"))
  matrix_sequences_plot$id <- matrix_sequences_plot$pos
  matrix_sequences_plot$allele <- gsub(paste0(g_group,"[*]"),"",matrix_sequences_plot$variable)
  matrix_sequences_plot$allele <- factor(matrix_sequences_plot$allele, levels = unique(matrix_sequences_plot$allele))
  matrix_sequences_plot$value[matrix_sequences_plot$value=="."] <- NA
  matrix_sequences_plot$annot_text <- sapply(1:nrow(matrix_sequences_plot), function(i) ifelse(matrix_sequences_plot$annot[i],matrix_sequences_plot$value[i],""))
  
  p1 <- ggplot(matrix_sequences_plot[matrix_sequences_plot$id<=80, ], aes(x=(pos), y=(allele))) + 
    geom_tile(aes(fill=value),colour="white") + 
    geom_text(aes(label = annot_text), color = "black") +
    coord_equal(expand = F, xlim = c(1, 80)) + bbplot::bbc_style()  +
    scale_fill_manual(values = c("#1380A1", "#FAAB18", "#990000", "#588300", "gray50")) + theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_blank(), panel.background = element_blank(), 
      legend.position = "none" )
  
  p2 <- ggplot(matrix_sequences_plot[matrix_sequences_plot$id>80 & matrix_sequences_plot$id <=160, ], aes(x=(pos), y=(allele))) + 
    geom_tile(aes(fill=value),colour="white") + 
    geom_text(aes(label = annot_text), color = "black") +
    coord_equal(expand = F, xlim = c(81, 160)) + bbplot::bbc_style()  + 
    scale_fill_manual(values = c("#1380A1", "#FAAB18", "#990000", "#588300", "gray50")) + theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_blank(), panel.background = element_blank(), 
      legend.position = "none" ) 
  
  p3 <- ggplot(matrix_sequences_plot[matrix_sequences_plot$id>160 & matrix_sequences_plot$id <=240, ], aes(x=(pos), y=(allele))) + 
    geom_tile(aes(fill=value),colour="white") + 
    geom_text(aes(label = annot_text), color = "black") +
    coord_equal(expand = F, xlim = c(161, 240)) + bbplot::bbc_style()  + 
    scale_fill_manual(values = c("#1380A1", "#FAAB18", "#990000", "#588300", "gray50")) + theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_blank(), panel.background = element_blank(), 
      legend.position = "none" ) 
  
  
  p4 <- ggplot(matrix_sequences_plot[matrix_sequences_plot$id>240, ], aes(x=(pos), y=(allele))) + 
    geom_tile(aes(fill=value),colour="white") + 
    geom_text(aes(label = annot_text), color = "black") +
    coord_equal(expand = F, xlim = c(240, 318)) + bbplot::bbc_style()  + 
    scale_fill_manual(values = c("#1380A1", "#FAAB18", "#990000", "#588300", "gray50")) + theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_blank(), panel.background = element_blank(), 
      legend.position = "none" ) 
  
  
  cowplot::plot_grid(p_dend, p1, p2, p3, p4, nrow=5, ncol = 1, rel_heights = c(0.4, 0.15, 0.15, 0.15, 0.15))
  
}
