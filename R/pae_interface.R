#' Perform kmeans clustering to identify amino acid cluster with the lowest PAE
#'
#' @param afm: list; result from import.afm() function
#' @param plldt: pLLDT cutoff
#'
#' @return data frame
#' @export
#'
#' @examples
#' import_afm <- import.afm()
#' pae_interface <- pae.interface()
pae.interface <- function(import_afm, plldt = 50) {
  set.seed(555)
  afm_pae <- import_afm$pae
  afm_plldt <- import_afm$plldt

  for(i in 1:length(afm_pae)){
    if(i == 1){
      cat(paste0("pae mean of stats::kmeans-clustering: ", i, " of ", length(afm_pae), "\n"))
      afm_pae_df <- data.frame(matrix(nrow = 0, ncol = 33))
    }
    if(i %% 10==0) {
      # Print on the screen some message
      cat(paste0("pae mean of stats::kmeans-clustering: ", i, " of ", length(afm_pae), "\n"))
    }

    A = unlist(str_split(unlist(str_extract_all(names(afm_pae)[i], pattern = "[^/]*$"))[1], pattern = "_"))[1:2][1]
    B = unlist(str_split(unlist(str_extract_all(names(afm_pae)[i], pattern = "[^/]*$"))[1], pattern = "_"))[1:2][2]
    proteinInfo <- import_afm$protein %>%
      dplyr::filter(A_protein == A & B_protein == B)
    plldt_filtered <- as.numeric(which(afm_plldt[[i]] > plldt))
    plldt_filteredB <- plldt_filtered-proteinInfo$A_length
    plldt_filteredB <- plldt_filteredB[which(plldt_filteredB > 0)]
    plldt_filteredA <- plldt_filtered[which(plldt_filtered<=proteinInfo$A_length)]
    if(length(plldt_filteredB) < 10) {
      qualityA = FALSE
    } else{qualityA = TRUE}
    if(length(plldt_filteredA) < 1) {
      qualityB = FALSE
    } else{qualityB = TRUE}
    afm_pae_interAB <- afm_pae[[i]][seq(proteinInfo$A_length+1,nrow(afm_pae[[i]])),seq(1,proteinInfo$A_length)]
    afm_pae_interAB <- afm_pae_interAB[plldt_filteredB,plldt_filteredA]
    afm_pae_interBA <- afm_pae[[i]][seq(1,proteinInfo$A_length),seq(proteinInfo$A_length+1,nrow(afm_pae[[i]]))]
    afm_pae_interBA <- afm_pae_interBA[plldt_filteredA,plldt_filteredB]

    if(qualityA == FALSE | qualityB == FALSE) {
      tmp <- data.frame(referenceset.af.proteinInfo,
                        file = names(afm_pae)[i],
                        complex = str_extract(names(afm_pae)[i], pattern = "[:alnum:]*_[:alnum:]*(?=_)"),
                        model = as.numeric(str_extract(names(afm_pae)[i], pattern = "(?<=_model_)[:digit:]*")),
                        rank = as.numeric(str_extract(names(afm_pae)[i], pattern = "(?<=_rank_)[:digit:]*")),
                        interAB.cluster11.median = NA,
                        interAB.cluster11.mean = NA,
                        interAB.cluster11.size = NA,

                        interAB.cluster22.median = NA,
                        interAB.cluster22.mean = NA,
                        interAB.cluster22.size = NA,

                        interAB.cluster21.median = NA,
                        interAB.cluster21.mean = NA,
                        interAB.cluster21.size = NA,

                        interAB.cluster12.median = NA,
                        interAB.cluster12.mean = NA,
                        interAB.cluster12.size = NA,

                        interBA.cluster11.median = NA,
                        interBA.cluster11.mean = NA,
                        interBA.cluster11.size = NA,

                        interBA.cluster22.median = NA,
                        interBA.cluster22.mean = NA,
                        interBA.cluster22.size = NA,

                        interBA.cluster21.median = NA,
                        interBA.cluster21.mean = NA,
                        interBA.cluster21.size = NA,

                        interBA.cluster12.median = NA,
                        interBA.cluster12.mean = NA,
                        interBA.cluster12.size = NA
      )
    } else {
      if(class(try(stats::kmeans(afm_pae_interAB, 2), silent = TRUE)) == "try-error"){
        afm_pae_interAB_cluster_row <- list(`1` = seq(1, nrow(afm_pae_interAB)),
                                                        `2` = seq(1, nrow(afm_pae_interAB)))
      } else{afm_pae_interAB_cluster_row <- with(data.frame(cluster = stats::kmeans(afm_pae_interAB, 2)$cluster) %>% dplyr::mutate(id = dplyr::row_number()), split(id, cluster))}

      if(class(try(stats::kmeans(t(afm_pae_interAB), 2), silent = TRUE)) == "try-error"){
        afm_pae_interAB.cluster.col <- list(`1` = seq(1, ncol(afm_pae_interAB)),
                                                        `2` = seq(1, ncol(afm_pae_interAB)))
      } else{afm_pae_interAB.cluster.col <- with(data.frame(cluster = stats::kmeans(t(afm_pae_interAB), 2)$cluster) %>% dplyr::mutate(id = dplyr::row_number()), split(id, cluster))}

      if(class(try(stats::kmeans(afm_pae_interBA, 2), silent = TRUE)) == "try-error"){
        afm_pae_interBA.cluster.row <- list(`1` = seq(1, nrow(afm_pae_interBA)),
                                                        `2` = seq(1, nrow(afm_pae_interBA)))
      } else{afm_pae_interBA.cluster.row <- with(data.frame(cluster = stats::kmeans(afm_pae_interBA, 2)$cluster) %>% dplyr::mutate(id = dplyr::row_number()), split(id, cluster))}

      if(class(try(stats::kmeans(t(afm_pae_interBA), 2), silent = TRUE)) == "try-error"){
        afm_pae_interBA.cluster.col <- list(`1` = seq(1, ncol(afm_pae_interBA)),
                                                        `2` = seq(1, ncol(afm_pae_interBA)))
      } else{afm_pae_interBA.cluster.col <- with(data.frame(cluster = stats::kmeans(t(afm_pae_interBA), 2)$cluster) %>% dplyr::mutate(id = dplyr::row_number()), split(id, cluster))}

      tmp <- data.frame(proteinInfo,
                        file = names(afm_pae)[i],
                        complex = str_extract(names(afm_pae)[i], pattern = "[:alnum:]*_[:alnum:]*(?=_)"),
                        model = as.numeric(str_extract(names(afm_pae)[i], pattern = "(?<=_model_)[:digit:]*")),
                        rank = as.numeric(str_extract(names(afm_pae)[i], pattern = "(?<=_rank_)[:digit:]*")),
                        interAB.cluster11.median = median(afm_pae_interAB[afm_pae_interAB_cluster_row$`1`, afm_pae_interAB.cluster.col$`1`], na.rm = TRUE),
                        interAB.cluster11.mean = mean(afm_pae_interAB[afm_pae_interAB_cluster_row$`1`, afm_pae_interAB.cluster.col$`1`], na.rm = TRUE),
                        interAB.cluster11.size = length(afm_pae_interAB_cluster_row$`1`) * length(afm_pae_interAB.cluster.col$`1`),

                        interAB.cluster22.median = median(afm_pae_interAB[afm_pae_interAB_cluster_row$`2`, afm_pae_interAB.cluster.col$`2`], na.rm = TRUE),
                        interAB.cluster22.mean = mean(afm_pae_interAB[afm_pae_interAB_cluster_row$`2`, afm_pae_interAB.cluster.col$`2`], na.rm = TRUE),
                        interAB.cluster22.size = length(afm_pae_interAB_cluster_row$`2`) * length(afm_pae_interAB.cluster.col$`2`),

                        interAB.cluster21.median = median(afm_pae_interAB[afm_pae_interAB_cluster_row$`2`, afm_pae_interAB.cluster.col$`1`], na.rm = TRUE),
                        interAB.cluster21.mean = mean(afm_pae_interAB[afm_pae_interAB_cluster_row$`2`, afm_pae_interAB.cluster.col$`1`], na.rm = TRUE),
                        interAB.cluster21.size = length(afm_pae_interAB_cluster_row$`2`) * length(afm_pae_interAB.cluster.col$`1`),

                        interAB.cluster12.median = median(afm_pae_interAB[afm_pae_interAB_cluster_row$`1`, afm_pae_interAB.cluster.col$`2`], na.rm = TRUE),
                        interAB.cluster12.mean = mean(afm_pae_interAB[afm_pae_interAB_cluster_row$`1`, afm_pae_interAB.cluster.col$`2`], na.rm = TRUE),
                        interAB.cluster12.size = length(afm_pae_interAB_cluster_row$`1`) * length(afm_pae_interAB.cluster.col$`2`),

                        interBA.cluster11.median = median(afm_pae_interBA[afm_pae_interBA.cluster.row$`1`, afm_pae_interBA.cluster.col$`1`], na.rm = TRUE),
                        interBA.cluster11.mean = mean(afm_pae_interBA[afm_pae_interBA.cluster.row$`1`, afm_pae_interBA.cluster.col$`1`], na.rm = TRUE),
                        interBA.cluster11.size = length(afm_pae_interBA.cluster.row$`1`) * length(afm_pae_interBA.cluster.col$`1`),

                        interBA.cluster22.median = median(afm_pae_interBA[afm_pae_interBA.cluster.row$`2`, afm_pae_interBA.cluster.col$`2`], na.rm = TRUE),
                        interBA.cluster22.mean = mean(afm_pae_interBA[afm_pae_interBA.cluster.row$`2`, afm_pae_interBA.cluster.col$`2`], na.rm = TRUE),
                        interBA.cluster22.size = length(afm_pae_interBA.cluster.row$`2`) * length(afm_pae_interBA.cluster.col$`2`),

                        interBA.cluster21.median = median(afm_pae_interBA[afm_pae_interBA.cluster.row$`2`, afm_pae_interBA.cluster.col$`1`], na.rm = TRUE),
                        interBA.cluster21.mean = mean(afm_pae_interBA[afm_pae_interBA.cluster.row$`2`, afm_pae_interBA.cluster.col$`1`], na.rm = TRUE),
                        interBA.cluster21.size = length(afm_pae_interBA.cluster.row$`2`) * length(afm_pae_interBA.cluster.col$`1`),

                        interBA.cluster12.median = median(afm_pae_interBA[afm_pae_interBA.cluster.row$`1`, afm_pae_interBA.cluster.col$`2`], na.rm = TRUE),
                        interBA.cluster12.mean = mean(afm_pae_interBA[afm_pae_interBA.cluster.row$`1`, afm_pae_interBA.cluster.col$`2`], na.rm = TRUE),
                        interBA.cluster12.size = length(afm_pae_interBA.cluster.row$`1`) * length(afm_pae_interBA.cluster.col$`2`)
      )
    }

    afm_pae_df <- rbind(afm_pae_df, tmp)

    rm(afm_pae_interAB_cluster_row, afm_pae_interAB.cluster.col, afm_pae_interBA.cluster.row, afm_pae_interBA.cluster.col)

    if(i == length(afm_pae)){
      afm_pae_df <- afm_pae_df %>%
        dplyr::rowwise() %>%
        dplyr::mutate(pae = min(dplyr::c_across(cols = contains("mean"))),
                      interaction = base::paste0(A_protein,"+",B_protein)) %>%
        dplyr::ungroup()
      cat(paste0("Done.", "\n"))
    }
  }
  return(afm_pae_df)
}
