#' Import PAE, pLLDT and PTM from AFM predictions .json files
#'
#' @param dir: specify directory where '*scores.json' files are stored
#'
#' @return list containing elements
#' @export
#'
#' @examples
#' afm_json <- import.afm()
import.afm <- function(dir = 'data/') {
  config_files <- list.files(path = dir, pattern = "config.json", recursive = TRUE)
  for(i in 1:length(config_files)) {
    if(i == 1) {
      config_json <- list()
      cat(paste0("read AFM `config.json` files: ", i, " of ", length(config_files), "\n"))
    }
    if(i %% 50==0) {
      # Print on the screen some message
      cat(paste0("read AFM `config.json` files: ", i, " of ", length(config_files), "\n"))
    }
    config_json[[i]] <- rjson::fromJSON(file = paste0(dir,"/",config_files[i]))
  }
  num_models <- c()
  for(i in 1:length(config_json)) {
    num_models <- c(num_models, config_json[[1]]$num_models)
  }

  json_files <- list.files(path = dir, pattern = "scores.json", recursive = TRUE)

  for(i in 1:length(json_files)) {
    if(i == 1) {
      afm_json <- list()
      cat(paste0("read AFM `*scores.json` files: ", i, " of ", length(json_files), "\n"))
    }
    if(i %% 50==0) {
      # Print on the screen some message
      cat(paste0("read AFM `*scores.json` files: ", i, " of ", length(json_files), "\n"))
    }
    afm_json[[i]] <- rjson::fromJSON(file = paste0(dir,"/",json_files[i]))
  }

  for(i in 1:length(afm_json)) {
    if(i == 1) {
      afm_ptm <- data.frame(matrix(nrow = 0, ncol = 3))
      colnames(afm_ptm) <- c("file", "ptm", "max_pae")
      afm_plldt <- list()
      afm_pae <- list()
      cat(paste0("extract ptm, plldt and pae values: ", i, " of ", length(json_files), "\n"))
    }
    if(i %% 50==0) {
      # Print on the screen some message
      cat(paste0("extract ptm, plldt and pae values: ", i, " of ", length(json_files), "\n"))
    }
    tmp <- data.frame(file = json_files[i],
                      ptm = afm_json[[i]]$ptm,
                      max_pae = afm_json[[i]]$max_pae)
    afm_ptm <- rbind(afm_ptm, tmp)
    afm_plldt[[i]] <- afm_json[[i]]$plddt
    afm_pae[[i]] <- matrix(unlist(afm_json[[i]]$pae),
                                       ncol = length(afm_plldt[[i]]),
                                       nrow = length(afm_plldt[[i]]))

    if(i == length(afm_json)) {
      names(afm_plldt) <- json_files
      names(afm_pae) <- json_files
      cat(paste0("Done.", "\n"))
    }
  }

  for(i in 1:length(afm_json)) {
    if(i == 1) {
      afm_ptm <- data.frame(matrix(nrow = 0, ncol = 3))
      colnames(afm_ptm) <- c("file", "ptm", "max_pae")
      afm_plldt <- list()
      afm_pae <- list()
      cat(paste0("extract ptm, plldt and pae values: ", i, " of ", length(json_files), "\n"))
    }
    if(i %% 50==0) {
      # Print on the screen some message
      cat(paste0("extract ptm, plldt and pae values: ", i, " of ", length(json_files), "\n"))
    }
    tmp <- data.frame(file = json_files[i],
                      ptm = afm_json[[i]]$ptm,
                      max_pae = afm_json[[i]]$max_pae)
    afm_ptm <- rbind(afm_ptm, tmp)
    afm_plldt[[i]] <- afm_json[[i]]$plddt
    afm_pae[[i]] <- matrix(unlist(afm_json[[i]]$pae),
                           ncol = length(afm_plldt[[i]]),
                           nrow = length(afm_plldt[[i]]))

    if(i == length(afm_json)) {
      names(afm_plldt) <- json_files
      names(afm_pae) <- json_files
      cat(paste0("Done.", "\n"))
    }
  }

  a3m_files <- list.files(path = dir, pattern = ".a3m", recursive = TRUE, full.names = TRUE)
  for(i in 1:length(a3m_files)) {
    if(i == 1) {
      a3m_json <- list()
      protein_info <- data.frame()
      cat(paste0("read `*.a3m` files: ", i, " of ", length(a3m_files), "\n"))
      protein_complex <- data.frame(file = unlist(str_extract_all(json_files[i], pattern = "[^/]*$"))[1])
    }
    if(i %% 50==0) {
      # Print on the screen some message
    }
    protein_info <- rbind(protein_info, readLines(a3m_files[i], n = 1))
    protein_complex <- rbind(protein_complex, unlist(str_extract_all(json_files[i+num_models[i-1]], pattern = "[^/]*$"))[1])

    if(i == length(a3m_files)) {
      colnames(protein_info) <- "protein_info"
      protein_info <- protein_info %>%
        tidyr::separate(protein_info, into = c("first_values", "second_values"), sep = "\t") %>%
        tidyr::separate(first_values, into = c("A_length", "B_length"), sep = ",") %>%
        dplyr::select(-second_values) %>%
        cbind(protein_complex) %>%
        dplyr::mutate(A_length = as.numeric(str_remove_all(A_length, "#")),
                      B_length = as.numeric(B_length)) %>%
        tidyr::separate(col = "file", into = c("A_protein", "B_protein"), sep = "_", extra = "drop", remove = TRUE)
      cat(paste0("read `*.a3m` files: ", i, " of ", length(a3m_files), "\n"))
    }
  }

  return(list(ptm = afm_ptm,
              plldt = afm_plldt,
              pae = afm_pae,
              protein = protein_info,
              num_models = num_models))
}
