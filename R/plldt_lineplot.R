#' Plot pLLDT as lineplot
#'
#' @param import_afm: list in the format from import.afm()
#' @param afm_complex: AFM predicted complex in the format "NSP10-NSP16"
#' @param afm_rank: "all" or integer: which afm_rank from AFM predictions to plot
#'
#' @return ggplot
#' @export
#'
#' @examples
plldt.lineplot <- function(import_afm, afm_complex = "NSP10-NSP16", afm_rank = "all") {
  set.seed = 555
  assertthat::assert_that(any(str_detect(names(import_afm$plldt), afm_complex)))
  assertthat::assert_that(afm_rank == "all" | is.double(afm_rank))
  proteinInfo <- import_afm$protein %>%
    dplyr::filter(A_protein == stringr::str_split(afm_complex, "-", simplify = TRUE)[1,1] & B_protein == stringr::str_split(afm_complex, "-", simplify = TRUE)[1,2])

  if(is.double(afm_rank)) {
    afm_rank <- as.integer(afm_rank)
    name <- names(import_afm$plldt)[stringr::str_detect(names(import_afm$plldt), afm_complex) & stringr::str_detect(names(import_afm$plldt), paste0("rank_", afm_rank))]
    plot <- data.frame(plldt = import_afm$plldt[[name]],
                       rank = afm_rank) %>% dplyr::mutate(id = dplyr::row_number()) %>%
      ggplot2::ggplot(ggplot2::aes(x = id, y = plldt, color = factor(rank))) +
      ggplot2::geom_line(size = 0.25) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~ stringr::str_replace(afm_complex, "-|_", "+")) +
      viridis::scale_color_viridis(option = "D", direction = 1, discrete = TRUE) +
      ggplot2::theme(legend.position = "right",
            text = ggplot2::element_text(size = 8, family = "Avenir"),
            axis.text = ggplot2::element_text(size = 8),
            axis.title = ggplot2::element_text(size = 10, family = "Avenir Medium"),
            strip.text = ggplot2::element_text(size = 10, family = "Avenir Medium"),
            legend.text = ggplot2::element_text(size = 8),
            strip.background = ggplot2::element_rect(color = "white"),
            panel.border = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(size = 0.5)) +
      ggplot2::geom_vline(xintercept = proteinInfo$A_length) +
      ggplot2::geom_hline(yintercept = 50) +
      ggplot2::labs(x = "amino acid", y = "pLLDT", fill = "pLLDT", col = "pLLDT") +
      ggplot2::lims(y = c(0,100))
  }

  if(afm_rank == "all") {
    ranks <- names(import_afm$plldt)[stringr::str_detect(names(import_afm$plldt), afm_complex)]
    df <- data.frame()
    for(i in 1:length(ranks)) {
      name <- names(import_afm$plldt)[stringr::str_detect(names(import_afm$plldt), afm_complex) & stringr::str_detect(names(import_afm$plldt), paste0("rank_", i))]
      df <- rbind(df, data.frame(plldt = import_afm$plldt[[name]],
                       rank = i) %>% dplyr::mutate(id = dplyr::row_number()))
      plot <- df %>%
        ggplot2::ggplot(ggplot2::aes(x = id, y = plldt, color = factor(rank))) +
        ggplot2::geom_line(size = 0.25) +
        ggplot2::theme_bw() +
        ggplot2::facet_wrap(~ stringr::str_replace(afm_complex, "-|_", "+")) +
        viridis::scale_color_viridis(option = "D", direction = -1, discrete = TRUE) +
        ggplot2::theme(legend.position = "right",
                       text = ggplot2::element_text(size = 8, family = "Avenir"),
                       axis.text = ggplot2::element_text(size = 8),
                       axis.title = ggplot2::element_text(size = 10, family = "Avenir Medium"),
                       strip.text = ggplot2::element_text(size = 10, family = "Avenir Medium"),
                       legend.text = ggplot2::element_text(size = 8),
                       strip.background = ggplot2::element_rect(color = "white"),
                       panel.border = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(size = 0.5)) +
        ggplot2::geom_vline(xintercept = proteinInfo$A_length) +
        ggplot2::geom_hline(yintercept = 50) +
        ggplot2::labs(x = "amino acid", y = "pLLDT", fill = "pLLDT", col = "rank") +
        ggplot2::lims(y = c(0,100))
    }

  }
  return(plot)
}
