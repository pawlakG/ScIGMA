get_config_paths <- function() {
  cfg_file <- app_sys("config.yml")
  if (cfg_file == "") cfg_file <- "inst/config.yml"
  cfg <- config::get(file = cfg_file)
  return(cfg$paths)
}
