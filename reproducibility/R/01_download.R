#' Download dataset from Zenodo
#' @param url Zenodo download URL
#' @param dest Local destination path
#' @return The path to the downloaded file
download_zenodo_h5 <- function(url, dest) {
    if (!dir.exists(dirname(dest))) {
        dir.create(dirname(dest), recursive = TRUE)
    }
    
    # Download file using httr to handle redirects cleanly
    response <- httr::GET(url, httr::write_disk(dest, overwrite = TRUE), httr::progress())
    httr::stop_for_status(response)
    
    return(dest)
}
