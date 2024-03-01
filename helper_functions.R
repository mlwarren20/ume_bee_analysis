###############################################
#   Author:       Magdalena Warren
#   Date created:  12 March 2020
###############################################


# Taxonomy wrangling ------------------------------------------------------
taxa_wrangle = function(ps, taxdf = NULL) {
  if (is.null(taxdf)) {
    taxdf = tax_table(ps)
  }
  if (!identical(taxa_names(ps), rownames(taxdf))) {
    rownames(taxdf) = attributes(rownames(taxdf))$names
  }
  df = tax_table(as.matrix(as.data.frame(taxdf) %>% 
                             mutate(Phylum = coalesce(Phylum, paste0(Kingdom, "_")), 
                                    Class = coalesce(Class, paste0(Phylum, "_")), 
                                    Order = coalesce(Order, paste0(Class, "_")), 
                                    Family = coalesce(Family, paste0(Order, "_")), 
                                    Genus = coalesce(Genus, paste0(Family, "_")), 
                                    Species = coalesce(Species, paste0(Genus, "_"))))) 
  df = sub("_.*", " sp.", df)
  tax_table(ps) = df
  return(ps)
}


# Beta diversity ----------------------------------------------------------
subset_pcoa <- function(physeq, d = "bray") {
  if (d %in% c("bray", "jaccard")) {
    ps = veganify(physeq)
    if (d == "bray") {
      dist = vegdist(otu_table(ps), method = d)
    } else {
      dist = vegdist(otu_table(ps), method = d, binary = T)
    }
  }
  
  if (d == "unifrac") {
    dist = UniFrac(physeq, weighted = F)
  }
  
  if (d == "wunifrac") {
    dist = UniFrac(physeq, weighted = T)
  }
  
  tmp_pcoa <- pcoa(dist)
  tmp_samData <- data.frame(sample_data(physeq))
  num_axes = ifelse (dim(tmp_pcoa$vectors)[2] > 10, 10, dim(tmp_pcoa$vectors)[2]) 
  tmp_samData <- cbind(tmp_samData, tmp_pcoa$vectors[,1:num_axes])
  # tmp_samData <- cbind(tmp_samData, tmp_pcoa$vectors[,1:num_axes], tmp_pcoa$values)
  row.names(tmp_samData) <- row.names(sample_data(physeq))
  sample_data(physeq) <- tmp_samData
  results = list("physeq" = physeq, "values" = tmp_pcoa$values, "dist" = dist)
  return(results)
}

# Make sure otu table has taxa as columns for vegan
veganify <- function(physeq) {
  if (taxa_are_rows(physeq)) {
    otu_table(physeq) = t(otu_table(physeq))
  }
  return(physeq)
}

## To add ellipses 
draw_ellipses <- function(parameter, param2 = NULL, np2 = NULL, physeq, axis1, axis2) {
  centroids <- aggregate(cbind(get(axis1), get(axis2)) ~ get(parameter), data.frame(sample_data(physeq)), mean)
  df <- do.call(rbind, lapply(na.omit(unique(data.frame(sample_data(physeq))[, parameter])), function(t) {
    data.frame(tag = as.character(t), 
               ellipse::ellipse(cov(na.omit(data.frame(sample_data(physeq))[data.frame(sample_data(physeq))[, parameter] == t, c(axis1, axis2)])),
                       centre = as.matrix(centroids[centroids[, "get(parameter)"] == t, 2:3]), 
                       level = 0.95), 
               stringsAsFactors = FALSE)}))
  names(df)[names(df) == "tag"] <- parameter
  if (!is.null(param2)) {
    df = df %>% mutate(np2 = param2)
    names(df)[names(df) == "np2"] <- np2
  }
  return(df)
}

plot_PCoA_w_highlight <- function(phyobj, phyname, values, colorcode, shapecode, linecode,
                                  color_scheme = bee_colors, line_scheme = NULL, shape_scheme = NULL,
                                  catg1, grp1, catg2 = NULL, grp2 = NULL, catg3 = NULL, grp3 = NULL, op = NULL) {
  DATA = data.frame(sample_data(phyobj))
  full_name = phyname
  if (is.null(catg2) & is.null(catg3)) {
    ellipse_ps = prune_samples(row.names(DATA %>% dplyr::filter(.data[[catg1]] == grp1)), phyobj)
  } else if (op == "and") {
    ellipse_ps = prune_samples(row.names(DATA %>% dplyr::filter(.data[[catg1]] == grp1 & .data[[catg2]] == grp2)), phyobj)
  } else if (op == "or") {
    ellipse_ps = prune_samples(row.names(DATA %>% dplyr::filter(.data[[catg1]] == grp1 | .data[[catg2]] == grp2)), phyobj)
  } 
  if(!is.null(catg3)) {
    ellipse_ps = prune_samples(sample_data(phyobj)$catg3 == grp3, ellipse_ps)
  }
  lapply(list(c(1, 2), c(2, 3), c(1, 3)), function(x) {
    xaxis = x[1]
    yaxis = x[2]
    
    x_name = paste0("Axis.", xaxis)
    y_name = paste0("Axis.", yaxis)
    
    ellipse_df = draw_ellipses(shapecode, param2 = grp1, np2 = catg1, physeq = ellipse_ps, axis1 = x_name, axis2 = y_name)
    if (grp1 == "nectar") {
      others = rbind(crossing(c("crop", "mouth"), unique(DATA[,catg2])), c("nectar", "winter"))
    } else {
      others = rbind(crossing(c("crop", "mouth"), unique(DATA[,catg2])))
    }
    
    f_df = do.call(rbind, apply(others, MARGIN = 1, function(x) {
      ps = prune_samples(DATA[,catg1] == x[[1]] & DATA[,catg2] == x[[2]], phyobj)
      df = draw_ellipses(shapecode, physeq = ps, axis1 = x_name, axis2 = y_name)
      return(df)
    }))
    
    total_lims = lapply(2:3, function(x) {
      c = c(DATA[,x_name], ellipse_df[,x], f_df[,x])
      return(list("high" = max(c), "low" = min(c)))
    })
    names(total_lims) = c("x", "y")
    
    axis1_exp <- round(values$Relative_eig[xaxis] * 100, 1)
    axis2_exp <- round(values$Relative_eig[yaxis] * 100, 1)
      
    pcoa_plot = ggplot(DATA, aes(x = .data[[x_name]], y = .data[[y_name]])) + 
      suppressWarnings(geom_point(aes(color = .data[[colorcode]], shape = .data[[shapecode]], text = sample_name), alpha = 0.7, size = 2)) +
      geom_vline(xintercept = 0, color = "darkgrey", linetype = "dashed") + 
      geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed") +
      xlim(total_lims$x$low, total_lims$x$high) +
      ylim(total_lims$y$low, total_lims$y$high) +
      scale_color_manual("Sample type", values = color_scheme) +
      labs(color = colorcode, shape = shapecode, line = shapecode) +
      theme_classic() + 
      theme(legend.title = element_blank()) +
      
      if (is.null(catg2) & is.null(catg3)) {
        gghighlight(.data[[catg1]] == grp1, use_group_by = FALSE) 
      } else if (op == "and") {
        if(is.null(catg3)) {
          gghighlight(.data[[catg1]] == grp1 & .data[[catg2]] == grp2, use_group_by = FALSE) 
        } else {
          gghighlight(.data[[catg1]] == grp1 & .data[[catg2]] == grp2 & .data[[catg3]] == grp3, use_group_by = FALSE) 
        }
      } else if (op == "or") {
        if(is.null(catg3)) {
          gghighlight(.data[[catg1]] == grp1 | .data[[catg2]] == grp2, use_group_by = FALSE) 
        } else {
          gghighlight(.data[[catg1]] == grp1 | .data[[catg2]] == grp2 | .data[[catg3]] == grp3, use_group_by = FALSE) 
        }
      } 
    pcoa_w_axis_lab = pcoa_plot + xlab(paste0("Axis ", xaxis, " (", axis1_exp, "%)")) + 
      ylab(paste0("Axis ", yaxis, " (", axis2_exp, "%)")) + 
      ggtitle(ifelse(!is.null(grp2), paste(grp2, grp1, "PCoA"), paste(grp1, "PCoA"))) +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    pcoa_w_ell = pcoa_w_axis_lab + geom_path(data = ellipse_df, aes(color = .data[[colorcode]], linetype = .data[[linecode]])) +
      scale_linetype_manual("", values = line_scheme) + theme(plot.title = element_blank()) 
    
    pcoas = list("pcoa" = pcoa_w_axis_lab, "pcoa_w_ell" = pcoa_w_ell)
    
    return(pcoas)
  })
}

# Function to do the permanova, betadisper, and tukeyhsd on betadisper real quick
perma_run = function(ps, parameter) {
  dist = UniFrac(ps, weighted = T)
  form = as.formula(paste0("dist ~ ", parameter))
  perma = adonis2(form, data = data.frame(sample_data(ps)))
  return(perma)
}

adonis_pairs = function(ps, parameter) {
  sams = unique(ps@sam_data[[parameter]])
  if (length(sams) <= 2) {
    return(NULL)
  }
  res = lapply(sams, function(not_include) {
    ps2 = prune_samples(ps@sam_data[[parameter]] != not_include, ps)
    perma = perma_run(ps2, parameter)
    included = sams[sams != not_include]
    tag = paste(included, collapse = "_")
    return(list("perma" = perma, "tag" = tag))
  })
  names(res) = unlist(lapply(res, function(x) x$tag))
  return(res)
}

betadiv_stats = function(ps, parameter) {
  dist = UniFrac(ps, weighted = T)
  form = as.formula(paste0("dist ~ ", parameter))
  perma = adonis2(form, data = data.frame(sample_data(ps)))
  perma_pair = adonis_pairs(ps, parameter)
  betadisp = (with(data.frame(sample_data(ps)), betadisper(dist, get(parameter))))
  anova_betadisp = anova(betadisp)
  tukey_betadisp = TukeyHSD(betadisp)
  return(list("permanova" = perma, "pair_permanova" = perma_pair, "betadisper" = betadisp, 
              "anova_betadisper" = anova_betadisp, "TukeyHSD_betadisper" = tukey_betadisp))
}


# Differential abundance --------------------------------------------------
dif_abund_deseq <- function(physeq, variable, t = "Wald", ft = "parametric", 
                            alpha = NULL, physeq_name, ...) {
  otu_table(physeq) <- otu_table(physeq) + 1
  dds = phyloseq_to_deseq2(physeq, reformulate(variable))
  dds = DESeq(dds, test = t, fitType = ft)
  res = results(dds, cooksCutoff = FALSE)
  if (!is.null(alpha)) {
    sigtab = res[which(res$padj < alpha), ]
  } else {
    sigtab = res
  }
  
  if (sigtab@nrows > 0) {
    sigtab = cbind(as(sigtab, "data.frame"), as(physeq@tax_table[rownames(sigtab), ], "matrix"))
    sigtab$Phylum = factor(as.character(sigtab$Phylum))
    sigtab$Genus = factor(as.character(sigtab$Genus))
    p = ggplot(sigtab, aes(x = log2FoldChange, y = Genus, color = Phylum, fill = Phylum)) + 
      geom_point(size = 1.5, shape = 21, alpha = 0.75) +
      ggtitle(paste0("DESeq2 Differential Abundance of '", physeq_name, "' \n(", res@elementMetadata@listData$description[2], ")")) + theme_bw() +
      theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
            plot.title = element_text(hjust = 0.5))
    results = list("dds" = dds, "results" = res, "df" = sigtab, "plot" = p)
    return(results)
  } else {
    return(physeq)
  }
}


# Alpha diversity ---------------------------------------------------------
### Functions for alpha diversity
alpha_wrangle <- function(df, tiss_vec) {
  if (length(tiss_vec) == 2) {
    new_df <- df %>% filter(origin == tiss_vec[1] & season == tiss_vec[2])
  } else {
    new_df <- df %>% filter(origin == tiss_vec)
  }
  return(new_df)
}

boot_res <- function(alpha_measure, bootfun, df) {
  form <- as.formula(paste0(alpha_measure, "~ (1|site_alpha)"))
  alpha_model <- lmer(form, data = df)
  set.seed(1279466746)
  res <- bootMer(alpha_model, bootfun, nsim = 2000, parallel = "multiprocess", ncpus = 6, re.form = ~ (1|site_alpha))
  colnames(res$t) <- c(rownames(ranef(alpha_model)[["site_alpha"]]), "Intercept")
  return(list("alpha_mod" = alpha_model, "results" = res))
}

ss_factored <- function(df) {
  new_df = df %>% 
    mutate(season = factor(ifelse(site %in% c("Ichiigawa", "Kiyokawa 1", "Oshine 1", "Nishi-Honjo_S"), "summer", "winter"), 
                           levels = c("summer", "winter")), 
           site = factor(site, levels = c("Ichiigawa", "Kiyokawa 1", "Oshine 1", "Nishi-Honjo_S","Nishi-Honjo_W", 
                                          "Higashi-honjo", "Hirono", "Kamihaya", "Kiyokawa 2", "Nagano", "Oshine 2", 
                                          "Shimanose", "Yamauchi")))
  return(new_df)
}

df_summary <- function(df) {
  new_df = ss_factored(df) %>% 
    group_by(site, season) %>% 
    summarize(mean = mean(value), med = median(value), lo = quantile(value, 0.025), hi = quantile(value, 0.975))
  return(new_df)
}

random_effect_df <- function(alpha_boot) {
  df = alpha_boot$t %>% 
    as_tibble() %>% 
    pivot_longer(!contains("Intercept"), names_to = "site")
  return(df_summary(df))
}

random_effect_plot <- function(ref_df, tissue) {
  df = ss_factored(ref_df)
  p = df %>% ggplot(aes(med, site)) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray50") + 
    geom_errorbar(aes(xmin = lo, xmax = hi)) + 
    geom_point(color = origin_colors[[tissue[1]]])
  return(p)
}

alpha_results_df <- function(alpha_res) {
  c = ncol(alpha_res$t)
  df = (alpha_res$t[,1:c-1] + alpha_res$t[,c]) %>% 
    as_tibble() %>% 
    pivot_longer(everything(), names_to = "site")
  return(df_summary(df))
}

alpha_results_plot <- function(rich_df, alpha_df, tissue, m) {
  p = rich_df %>% ggplot(aes(med, site)) +
    geom_jitter(data = alpha_df, aes(x = get(m), y = site_alpha, shape = bee_species), 
                alpha = 0.7, width = 0.2, color = origin_colors[[tissue[1]]], size = 2) +
    scale_shape_manual(values = bee_shapes) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2) + 
    geom_point(size = 3, shape = origin_shapes[[tissue[1]]], fill = origin_colors[[tissue[1]]]) +
    xlab(m) + 
    theme_classic() + 
    theme(axis.title.y = element_blank())
}


# Clam test plots ---------------------------------------------------------
clam_plot = function(df1, df2 = NULL, title = NULL, y_lab = T) {
  cn = colnames(df1)
  lim = plyr::round_any(max(df1[2:3]), 1000, f = ceiling)
  p = ggplot(df1, aes(x = get(cn[2]) + 1, y = get(cn[3]) + 1, color = Classes, fill = Classes)) + 
    geom_point(alpha = 0.8, shape = 21) + 
    scale_color_manual(values = class_colors) +
    scale_fill_manual(values = class_fills) +
    scale_x_log10(limits = c(1, 100000), breaks = c(1e+01, 1e+03, 1e+05), labels = c("1e+01", "1e+03", "1e+05")) + 
    scale_y_log10(limits = c(1, 100000), breaks = c(1e+01, 1e+03, 1e+05), labels = c("1e+01", "1e+03", "1e+05")) +
    geom_abline(slope = 1, intercept = 0, color = "darkgrey", linetype = "dashed") +
    xlab(paste0("Log(", str_to_sentence(str_split_i(cn[2], "_", 2)), " OTU abundance + 1)")) +
    ggtitle(title) + theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "bold", 
                                    size = 16, margin = margin(t = 0, r = 0, b = 20, l = 0)),
          legend.position = "none",
          axis.text = element_text(size = 10), axis.title = element_text(size = 12))
  if (y_lab == T) {
    p + ylab(paste0("Log(", str_to_sentence(str_split_i(cn[3], "_", 2)), " OTU abundance + 1)"))
  } else {
    p + theme(axis.title.y = element_blank())
  }
}

clam_barplot = function(df, title = NULL, y_lab = T) {
  p = df %>%  dplyr::count(Classes) %>% ungroup() %>% mutate(freq = n/sum(n)) %>%
    ggplot(aes(x = Classes, y = freq, fill = Classes, color = Classes)) + 
    geom_bar(aes(fill = Classes, color = Classes), stat = "identity") +
    geom_text(aes(label = n), color = "black", vjust = -0.4, size = 2.75) +
    scale_fill_manual(values = class_fills) +
    scale_color_manual(values = class_colors) +
    scale_x_discrete(labels = clam_labs) +
    ylim(0, 1) +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 12), 
          axis.text = element_text(size = 10)) 
  if (y_lab == T) {
    p + ylab("Proportion of OTUs Classified")  
  } else {
    p + theme(axis.title.y = element_blank())
  }
}



