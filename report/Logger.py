import logging
import multiprocess


def configure_logger(filename):
    """
    Configure logger file.
    :param filename: str - filename where to log
    :return: None
    """
    format_str = '%(asctime)s %(levelname)10s %(pid)10d: %(message)s'
    logging.basicConfig(filename=filename, level=logging.DEBUG, format=format_str, filemode="w")
    logging.getLogger().setLevel(logging.INFO)


def log_str(to_log, stdout_too=True, priority=logging.INFO, flush=False):
    """
    Write a string to log.
    :param to_log: str - string to write to log
    :param stdout_too: bool - write the string to stdout too
    :param priority: int - logging priority
    :param flush: bool - flush this logger?
    :return: None
    """
    if stdout_too:
        print(to_log)

    # extract pid
    pid = multiprocess.current_process().pid

    # delete \r and make it one-liners:
    to_log = to_log.replace('\r', '')
    to_log = to_log.replace('\n', '    ')

    # finally log it!
    logging.log(priority, to_log, extra={"pid": pid})

    # flush it?
    if flush:
        logging.getLogger().handlers[0].flush()
