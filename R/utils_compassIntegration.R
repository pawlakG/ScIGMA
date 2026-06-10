#' Run Compiled COMPASS Binary
#'
#' @description
#' Locates the compiled C++ binary within the package installation directory
#' and executes it as a system process.
#'
#' @param input_file Character. Path to the preprocessed AML data.
#' @param output_prefix Character. Prefix for the output results.
#'
#' @return Integer exit code from the system process (0 indicates success).
#' @noRd
execute_compass_binary <- function(input_file, output_prefix) {
    # system.file will look into the 'inst' folder of the installed package
    binary_path <- system.file("bin", "COMPASS", package = "your_package")

    if (binary_path == "") {
        stop("COMPASS binary not found. Compilation failed during installation.")
    }

    cli_arguments <- c(
        "-i", input_file,
        "-o", output_prefix,
        "--nchains", "4",
        "--chainlength", "5000",
        "--CNA", "1"
    )

    # Execute the binary
    exit_status <- system2(command = binary_path, args = cli_arguments)

    return(exit_status)
}
