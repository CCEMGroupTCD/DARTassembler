import logging
import pathlib as pl


def setup_logger(name: str = "DARTassembler", log_file: pl.Path = pl.Path("assembly.log"), level=logging.DEBUG) -> logging.Logger:
    """
    Sets up a logger for the DART assembly process.
    :param name: The name of the logger, defaults to "DART"
    :param log_file: The name of the log file, defaults to "assembly.log"
    :param level: The level of logging, defaults to INFO
    :return: The configured logger
    """

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid duplicate handlers
    if not logger.handlers:
        formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(name)s: %(message)s')

        # Console output
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        # File output
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


# Create a default global logger instance
default_logger = setup_logger()

# Custom exception that auto-logs
class LoggedValueError(ValueError):
    def __init__(self, message, logger=default_logger):
        logger.error(message)
        super().__init__(message)