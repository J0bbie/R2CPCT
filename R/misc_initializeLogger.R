#' @title Initialize a parallel-ready logger.
#'
#' @param name (character): Name of the logger (can be empty).
#' @param file (character): Path to file where logging will be written to (handy when performing parallel).
#'
#' @examples
#' \donttest{
#'
#'  initializeLogger()
#'
#' }
#' @return (NULL) Returns instance of a ParallelLogger.
#' @export
initializeLogger <- function(name = NULL, file = NULL){

    # Input validation --------------------------------------------------------

    checkmate::assertCharacter(name, null.ok = TRUE)
    if(!is.null(file)) checkmate::checkFile(file, access = 'rw')


    # Initialize Logger -------------------------------------------------------

    ParallelLogger::clearLoggers()

    ParallelLogger::registerLogger(
        ParallelLogger::createLogger(
            name = name,
            threshold = "INFO",
            appenders = list(ParallelLogger::createConsoleAppender(layout = ParallelLogger::layoutTimestamp)))
    )

    # (Also) create a file logger.
    if(!is.null(file)){
        ParallelLogger::registerLogger(
            ParallelLogger::createLogger(
                name = name,
                threshold = "INFO",
                appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutTimestamp, file = file, overwrite = TRUE)))
        )
    }
}
