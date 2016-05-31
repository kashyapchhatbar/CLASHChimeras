import logging
import coloredlogs


def debug_logger(name):
    logger = logging.getLogger(name)
    coloredlogs.install(level=logging.DEBUG, show_name=False,
                        show_hostname=True)
    return logger


def info_logger(name):
    logger = logging.getLogger(name)
    coloredlogs.install(level=logging.INFO, show_name=False,
                        show_hostname=True)
    return logger


def warning_logger(name):
    logger = logging.getLogger(name)
    coloredlogs.install(level=logging.WARNING, show_name=False,
                        show_hostname=True)
    return logger


def error_logger(name):
    logger = logging.getLogger(name)
    coloredlogs.install(level=logging.ERROR, show_name=False,
                        show_hostname=True)
    return logger
