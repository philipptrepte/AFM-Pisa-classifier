#' Boxplot of PISA results
#'
#' @param import_pisa
#'
#' @return ggplot2 boxplot
#' @export
#'
#' @examples
pisa.boxplot <- function(import_pisa) {
  import_pisa %>%
    dplyr::mutate(complex = str_replace_all(complex, "_", "-")) %>%
    dplyr::mutate("-deltaG" = -deltaG) %>%
    tidyr::pivot_longer(c(interfaceArea, "-deltaG", surfaceArea), names_to = "data", values_to = "score") %>%
    ggplot2::ggplot(ggplot2::aes(x = complex, y = score)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 1, fill = "#8ED2C6", color = "black") +
    ggplot2::geom_jitter(shape = 21, alpha = 1, fill = "white", size = 0.75) +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~ factor(data, levels = c("-deltaG", "interfaceArea", "surfaceArea")), scales = "free", ncol = 3) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
          axis.text = ggplot2::element_text(size = 8, color = "black"),
          axis.title = ggplot2::element_text(size = 10, family = "Avenir Medium"),
          text = ggplot2::element_text(family = "Avenir"),
          strip.text = ggplot2::element_text(size = 10, family = "Avenir Medium"),
          strip.background = ggplot2::element_blank(),
          legend.position = "none",
          axis.line = ggplot2::element_line(color = "black"),
          axis.ticks = ggplot2::element_line(color = "black"))
}
