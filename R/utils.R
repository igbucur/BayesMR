#' A more portable function for invoking a system command
#'
#' @param command String containing the system command to be invoked.
#' @param Windows_shell Full path to the Windows shell to be used for invoking 
#' the system command. MSYS2+MinGW has been tested and is recommended.
#' @param ... Further shell parameters.
#'
#' @return See \link[base]{system} for Linux and \link[base]{shell} for Windows.
#' @export
#'
#' @examples run_system_command("echo Hello!")
run_system_command <- function(command, Windows_shell = 'C:/msys64/msys2_shell.cmd -defterm -here -no-start -mingw64', ...) {
  
  if(.Platform$OS.type == "windows") {
    shell(shQuote(command), Windows_shell, flag = "-c", ...)
  } else {
    system(command)
  }
}