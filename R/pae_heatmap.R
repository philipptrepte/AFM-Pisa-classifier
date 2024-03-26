#' Plot kmeans interface clusters
#'
#' @param import_afm: list in the format from import.afm()
#' @param afm_complex: AFM predicted complex in the format "NSP10-NSP16"
#' @param afm_rank: which afm_rank from AFM predictions to plot
#' @param interface_cluster: cluster PAE values
#' @param plldt: pLLDT cutoff
#'
#' @return ComplexHeatmap
#' @export
#'
#' @examples
pae.heatmap <- function(import_afm, afm_complex = "NSP14-NSP10", afm_rank = 1, interface_cluster = FALSE, plldt = 50) {
  set.seed = 555
  assertthat::assert_that(any(str_detect(names(import_afm$pae), afm_complex)))
  name <- names(import_afm$pae)[str_detect(names(import_afm$pae), afm_complex) & str_detect(names(import_afm$pae), paste0("rank_0*", afm_rank))]
  afm_plldt <- import_afm$plldt
  afm_pae <- import_afm$pae
  
  pae_matrix <- import_afm$pae[[name]]
  length_A <- import_afm$protein %>%
    dplyr::filter(A_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[1]] & B_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[2]]) %>%
    dplyr::pull(A_length)
  protein_A <- import_afm$protein %>%
    dplyr::filter(A_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[1]] & B_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[2]]) %>%
    dplyr::pull(A_protein)
  length_B <- import_afm$protein %>%
    dplyr::filter(A_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[1]] & B_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[2]]) %>%
    dplyr::pull(B_length)
  protein_B <- import_afm$protein %>%
    dplyr::filter(A_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[1]] & B_protein == str_split(afm_complex, pattern = "-|_|\\+")[[1]][[2]]) %>%
    dplyr::pull(B_protein)
  
  if(interface_cluster == FALSE) {
    plot <- ComplexHeatmap::Heatmap(pae_matrix, cluster_rows = FALSE, cluster_columns = FALSE,
                                    row_split = c(rep(protein_A, length_A), rep(protein_B, length_B)),
                                    column_split = c(rep(protein_A, length_A), rep(protein_B, length_B)),
                                    border = TRUE, use_raster = TRUE,
                                    heatmap_legend_param = list(title = "PAE (Å)"))
  } else if(interface_cluster == TRUE) {
    plldt_filtered <- as.numeric(which(afm_plldt[[name]] > plldt))
    plldt_filteredB <- plldt_filtered-length_A
    plldt_filteredB <- plldt_filteredB[which(plldt_filteredB > 0)]
    plldt_filteredA <- plldt_filtered[which(plldt_filtered <= length_A)]
    
    afm_pae_interAB <- afm_pae[[name]][seq(length_A+1, nrow(afm_pae[[name]])), seq(1,length_A)]
    afm_pae_interAB <- afm_pae_interAB[plldt_filteredB, plldt_filteredA]
    afm_pae_interBA <- afm_pae[[name]][seq(1,length_A),seq(length_A+1,nrow(afm_pae[[name]]))]
    afm_pae_interBA <- afm_pae_interBA[plldt_filteredA,plldt_filteredB]
    
    heatmap_interAB <- ComplexHeatmap::Heatmap(afm_pae_interAB, cluster_rows = TRUE, cluster_columns = TRUE,
                                               border = TRUE, use_raster = TRUE, row_km = 2, column_km = 2,
                                               heatmap_legend_param = list(title = "PAE (Å)"),
                                               row_title = NULL, column_title = NULL, circlize::colorRamp2(c(0, 15, 30), c("blue", "white", "red")))
    heatmap_interBA <- ComplexHeatmap::Heatmap(afm_pae_interBA, cluster_rows = TRUE, cluster_columns = TRUE,
                                               border = TRUE, use_raster = TRUE, row_km = 2, column_km = 2,
                                               heatmap_legend_param = list(title = "PAE (Å)"),
                                               row_title = NULL, column_title = NULL, circlize::colorRamp2(c(0, 15, 30), c("blue", "white", "red")))
    
    afm_pae_df <- pae.interface(import_afm, plldt = plldt)
    
    barplot <- afm_pae_df %>%
      dplyr::filter(complex == str_replace(afm_complex, "-", "_") | complex == afm_complex) %>%
      tidyr::pivot_longer(cols = contains("cluster")) %>%
      dplyr::select(-file) %>%
      tidyr::separate(col = "name", into = c("chain", "cluster", "data"), sep = "\\.") %>%
      tidyr::unite(col = "name", chain, cluster) %>%
      dplyr::filter(data == "mean", rank == afm_rank) %>%
      unique() %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(name, value, decreasing = TRUE), y = value, fill = value)) +
      ggplot2::geom_bar(stat = "identity", color = "black") +
      ggplot2::scale_fill_gradientn(colors = c("blue", "white", "red"), limit = c(0, NA)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none",
                     text = ggplot2::element_text(size = 8, family = "Avenir"),
                     axis.text = ggplot2::element_text(size = 8),
                     axis.text.x = ggplot2::element_blank(),
                     axis.title = ggplot2::element_text(size = 10, family = "Avenir Medium"),
                     strip.text = ggplot2::element_text(size = 10),
                     legend.text = ggplot2::element_text(size = 8),
                     strip.background = ggplot2::element_rect(color = "white"),
                     panel.border = ggplot2::element_blank(),
                     axis.line = ggplot2::element_line(size = 0.5)) +
      ggplot2::labs(x = "cluster", y = "mean PAE (Å)", fill = "PAE")
    
    if(length_A > length_B) {
      plot <- cowplot::plot_grid(
        cowplot::plot_grid(grid.grabExpr(draw(heatmap_interBA)), NULL, ncol = 1, rel_heights = c(1,0)),
        cowplot::plot_grid(grid.grabExpr(draw(heatmap_interAB)), NULL,
                           cowplot::plot_grid(barplot, NULL, ncol = 2),
                           ncol = 1, rel_heights = c(ifelse(nrow(afm_pae_interAB)/nrow(afm_pae_interBA) < 0.5, 0.5, nrow(afm_pae_interAB)/nrow(afm_pae_interBA)), 0.5, 0.5)),
        ncol = 2, rel_widths = c(1, ifelse(ncol(afm_pae_interAB)/ncol(afm_pae_interBA) > 5, 5, ncol(afm_pae_interAB)/ncol(afm_pae_interBA))))
    }
    if(length_B > length_A) {
      plot <- cowplot::plot_grid(
        cowplot::plot_grid(grid.grabExpr(draw(heatmap_interAB)), NULL, ncol = 1, rel_heights = c(1,0)),
        cowplot::plot_grid(grid.grabExpr(draw(heatmap_interBA)), NULL,
                           cowplot::plot_grid(barplot, NULL, ncol = 2),
                           ncol = 1, rel_heights = c(ifelse(nrow(afm_pae_interBA)/nrow(afm_pae_interAB) < 0.5, 0.5, nrow(afm_pae_interBA)/nrow(afm_pae_interAB)), 0.5, 0.5)),
        ncol = 2, rel_widths = c(1, ifelse(ncol(afm_pae_interBA)/ncol(afm_pae_interAB) > 5, 5, ncol(afm_pae_interBA)/ncol(afm_pae_interAB))))
    }
    
  }
  
  return(plot)
}

