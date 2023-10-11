#' Import interface summary from PISA .xml files
#'
#' @param dir: directory where .xml files are saved
#'
#' @return data.frame
#' @export
#'
#' @examples
#' import_pisa <- import.pisa(dir="data")
import.pisa <- function(dir = "data") {
  files <- list.files(path = dir, recursive = TRUE, full.names = TRUE, pattern = "(interface)")
  for(i in 1:length(files)) {
    if(i == 1){
      pisa_info <- list()
    }
    if(i %% 100==0) {
      # Print on the screen some message
      cat(paste0("read xml files:", i, " of ", length(files), "\n"))
    }
    pisa_info[[i]] <- xml2::as_list(xml2::read_xml(x = files[i]))
    if(i == max(length(files))){
      names(pisa_info) <- stringr::str_remove(files, pattern = ".*/xml_files")
    }
  }

  for(i in seq(from = 1, to = length(files), by = 2)) {
    if(i == 1) {
      interface_info <- data.frame(matrix(nrow = 0, ncol = 15))
    }
    interface_info <- rbind(interface_info,
                            data.frame(complex = str_extract(names(pisa_info)[i], pattern = "[:alnum:]*_[:alnum:]*(?=_)"),
                                       model = as.numeric(str_extract(names(pisa_info)[i], pattern = "(?<=_model_)[:digit:]*")),
                                       rank = as.numeric(str_extract(names(pisa_info)[i], pattern = "(?<=_rank_)[:digit:]*")),
                                       surfaceAreaA = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$TOTALSURFACEAREA1),
                                       surfaceAreaB = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$TOTALSURFACEAREA2),
                                       interfaceAreaA = as.numeric(pisa_info[[i]]$INTERFACESUMMARY$STRUCTURE1$SOLVENTAREA1$INTERFACEAREA),
                                       interfaceAreaB = as.numeric(pisa_info[[i]]$INTERFACESUMMARY$STRUCTURE2$SOLVENTAREA2$INTERFACEAREA),
                                       interfaceResiduesA = as.numeric(pisa_info[[i]]$INTERFACESUMMARY$STRUCTURE1$NUMBEROFREDSIDUES1$NINTERFACEREDSIDUES),
                                       interfaceResiduesB = as.numeric(pisa_info[[i]]$INTERFACESUMMARY$STRUCTURE2$NUMBEROFREDSIDUES2$NINTERFACEREDSIDUES),
                                       deltaG = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACEAREA),
                                       pvalue = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACEDELTAGPVALUE),
                                       hbonds = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACENHBONDS),
                                       saltbridges = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACENSALTBRIDGES),
                                       disulfide = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACENDISULFIDEBONDS),
                                       css = as.numeric(pisa_info[[i+1]]$INTERFACELIST$INTERFACEROW$INTERFACECSS)
                            ))
    if(i == max(length(files))-1) {
      interface_info <- interface_info %>%
        dplyr::rowwise() %>%
        dplyr::mutate(interfaceArea = rowMeans(dplyr::across(contains("interfaceArea")), na.rm = TRUE),
                      surfaceArea = sum(x = dplyr::across(contains("surfaceArea")), na.rm = TRUE)) %>%
        dplyr::ungroup()
    }
  }
  return(interface_info)
}
