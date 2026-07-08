#' Get Configuration Paths
#'
#' Retrieves path configurations from the `config.yml` file. It searches for the
#' configuration file within the package context and falls back to a local `inst/config.yml`
#' if not found.
#'
#' @return A list of paths specified in the configuration file.
#' @export
get_config_paths <- function() {
    cfg_file <- app_sys("config.yml")
    if (cfg_file == "") cfg_file <- "inst/config.yml"
    cfg <- config::get(file = cfg_file)
    return(cfg$paths)
}
