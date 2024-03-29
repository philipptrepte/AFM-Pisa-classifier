#' Heatmap for PAE, deltaG, interfaceArea and surfaceArea for predicted AFM models
#'
#' @param import_afm
#' @param import_pisa
#'
#' @return ggplot geom_tile
#' @export
#'
#' @examples
afm.pisa.heatmap <- function(import_afm, import_pisa, plldt = 0) {
  interface <- import_pisa
  afm_pae <- pae.interface(import_afm, plldt = plldt)

  plot <- cowplot::plot_grid(
    interface %>%
      right_join(afm_pae %>% dplyr::select(complex, A_protein, B_protein, model, rank, pae), by=c("complex", "model", "rank")) %>%
      dplyr::mutate(interaction = paste(A_protein,"+", B_protein),
                    pae = ifelse(pae > 30, 30, pae)) %>%
      ggplot(aes(x = rank, y = interaction, fill = pae)) +
      geom_tile() +
      scale_fill_gradientn(colours = c("blue", "white", "red"), breaks = c(0, 10, 20, 30), limits = c(0,30)) +
      theme(text = element_text(size = 8),
            axis.text.y = element_text(size = 6),
            axis.text.x = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 8),
            panel.background = element_rect(fill="white"),
            legend.position = "top") +
      labs(fill = "PAE(Å)"),

    interface %>%
      right_join(afm_pae %>% dplyr::select(complex, A_protein, B_protein, model, rank, pae), by=c("complex", "model", "rank")) %>%
      dplyr::mutate(interaction = paste(A_protein,"+", B_protein),
                    deltaG = ifelse(deltaG < -50, -50, deltaG)) %>%
      ggplot(aes(x = rank, y = interaction, fill = deltaG)) +
      geom_tile() +
      scale_fill_gradientn(colours = viridisLite::viridis(option = "C", n = 100), labels = c("0", "-10", "-20", "-30", "-40", "< -50"),
                           breaks = c(0, -10, -20, -30, -40, -50), limits = c(NA,0)) +
      theme(text = element_text(size = 8),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 8),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill="white"),
            legend.position = "top") +
      labs(fill = "∆G"),

    interface %>%
      right_join(afm_pae %>% dplyr::select(complex, A_protein, B_protein, model, rank, pae), by=c("complex", "model", "rank")) %>%
      dplyr::mutate(interaction = paste(A_protein,"+", B_protein),
                    interfaceArea = ifelse(interfaceArea >3000, 3000, interfaceArea),
                    interfaceArea = interfaceArea/1000) %>%
      ggplot(aes(x = rank, y = interaction, fill = interfaceArea)) +
      geom_tile() +
      scale_fill_gradientn(colours = viridisLite::viridis(n = 100), labels = c("0", "1", "2", ">3"), limits = c(0,3),
                           breaks = c(0, 1, 2, 3)) +
      theme(text = element_text(size = 8),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 8),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill="white"),
            legend.position = "top") +
      labs(fill = "iA"),

    interface %>%
      right_join(afm_pae %>% dplyr::select(complex, A_protein, B_protein, model, rank, pae), by=c("complex", "model", "rank")) %>%
      dplyr::mutate(interaction = paste(A_protein,"+", B_protein)) %>%
      ggplot(aes(x = rank, y = interaction, fill = log2(surfaceArea))) +
      geom_tile() +
      scale_fill_gradientn(colours = viridisLite::viridis(n = 100, option = "E"), limits = c(0,NA)) +
      theme(text = element_text(size = 8),
            axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 8),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.background = element_rect(fill="white"),
            legend.position = "top")  +
      labs(fill = "sA"),

    align = "hv", axis = "bt", ncol = 4)

  return(plot)
}
