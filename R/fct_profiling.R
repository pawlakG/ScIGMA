#' Profiling helper for ScIGMA
#'
#' @description
#' Wraps an expression to measure time and memory RAM impact,
#' and prints a report.
#'
#' @param step_name Character. The name of the profiled step.
#' @param expr Expression. The R code to profile.
#' @param filepath Character (optional). Path to the loaded file to compute RAM/File size ratio.
#' @return The result of the expression.
#' @noRd
ScIGMA_profile <- function(step_name, expr, filepath = NULL) {
    # File size in MB
    file_size_mb <- NA
    if (!is.null(filepath) && file.exists(filepath)) {
        file_size_mb <- file.size(filepath) / (1024^2)
    }

    # RAM before
    gc_info_before <- gc(reset = TRUE)
    gc_before <- sum(gc_info_before[, 2])

    # Time before
    t0 <- Sys.time()

    # Execute expression in calling environment
    res <- withVisible(force(expr))

    # Time after
    t1 <- Sys.time()
    time_taken <- as.numeric(difftime(t1, t0, units = "secs"))

    # RAM after
    gc_info_after <- gc()
    gc_after <- sum(gc_info_after[, 2])
    ram_diff <- gc_after - gc_before

    # Ratio
    ratio_str <- "N/A"
    if (!is.na(file_size_mb) && file_size_mb > 0) {
        ratio_str <- sprintf("%.2f", ram_diff / file_size_mb)
    }

    # Machine specs
    sys_info <- Sys.info()
    os <- sys_info["sysname"]
    cpu <- "Unknown"
    tot_ram <- "Unknown"
    if (os == "Darwin") {
        cpu <- tryCatch(system("sysctl -n machdep.cpu.brand_string", intern = TRUE), error = function(e) "Unknown")
        tot_ram_bytes <- tryCatch(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)), error = function(e) NA)
        if (!is.na(tot_ram_bytes)) tot_ram <- sprintf("%.2f GB", tot_ram_bytes / (1024^3))
    } else if (os == "Linux") {
        cpu <- tryCatch(system("lscpu | grep 'Model name' | sed -r 's/Model name:\\s+//g'", intern = TRUE), error = function(e) "Unknown")
        tot_ram_kb <- tryCatch(as.numeric(system("awk '/MemTotal/ {print $2}' /proc/meminfo", intern = TRUE)), error = function(e) NA)
        if (!is.na(tot_ram_kb)) tot_ram <- sprintf("%.2f GB", tot_ram_kb / (1024^2))
    } else if (os == "Windows") {
        cpu <- tryCatch(system("wmic cpu get name", intern = TRUE)[2], error = function(e) "Unknown")
        cpu <- trimws(cpu)
        tot_ram_bytes <- tryCatch(as.numeric(system("wmic computersystem get TotalPhysicalMemory", intern = TRUE)[2]), error = function(e) NA)
        if (!is.na(tot_ram_bytes)) tot_ram <- sprintf("%.2f GB", tot_ram_bytes / (1024^3))
    }

    report <- sprintf(
        "=== ScIGMA Performance Report ===\nStep: %s\nTime: %.2f seconds\nRAM Before: %.2f MB\nRAM After: %.2f MB\nRAM Impact (Delta): %.2f MB\nFile Size: %s MB\nRAM/File Size Ratio: %s\nOS: %s\nCPU: %s\nTotal RAM: %s\n=================================\n",
        step_name, time_taken, gc_before, gc_after, ram_diff,
        ifelse(is.na(file_size_mb), "N/A", sprintf("%.2f", file_size_mb)),
        ratio_str,
        os, cpu, tot_ram
    )

    message(report)
    if (res$visible) {
        return(res$value)
    } else {
        return(invisible(res$value))
    }
}
