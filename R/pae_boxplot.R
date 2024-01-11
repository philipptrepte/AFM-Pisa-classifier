#' Boxplot to visualize PAE for all models of a AFM predicted complex
#'
#' @param pae_interface: object from the pae.interface() function
#'
#' @return boxplot
#' @export
#'
#' @examples
pae.boxplot <- function(pae_interface) {
  pae_interface %>%
    ggplot2::ggplot(ggplot2::aes(x = interaction, y = pae)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 1, fill = "#8ED2C6", color = "black", stroke = 0.5) +
    ggplot2::geom_jitter(shape = 21, alpha = 1, fill = "white", size = 0.75) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = ggplot2::element_text(size = 8, color = "black"),
          axis.title = ggplot2::element_text(size = 10),
          text = ggplot2::element_text(family = "Avenir"),
          legend.position = "none",
          axis.line = ggplot2::element_line(color = "black"), axis.ticks = ggplot2::element_line(color = "black")) +
    ggplot2::labs(y = "PAE (Å)")
}
