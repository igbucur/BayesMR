% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{run_system_command}
\alias{run_system_command}
\title{A more portable function for invoking a system command}
\usage{
run_system_command(
  command,
  Windows_shell = "C:/msys64/msys2_shell.cmd -defterm -here -no-start -mingw64",
  ...
)
}
\arguments{
\item{command}{String containing the system command to be invoked.}

\item{Windows_shell}{Full path to the Windows shell to be used for invoking
the system command. MSYS2+MinGW has been tested and is recommended.}

\item{...}{Further shell parameters.}
}
\value{
See \link[base]{system} for Linux and \link[base]{shell} for Windows.
}
\description{
A more portable function for invoking a system command
}
\examples{
run_system_command("echo Hello!")
}
