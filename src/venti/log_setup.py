# SPDX-FileCopyrightText: 2025-2026 opera-adt
# SPDX-License-Identifier: BSD-3-Clause
# Part of Venti, https://github.com/opera-adt/Venti (BSD-3-Clause, see LICENSE).
"""Logging setup shared by every way of running Venti.

Every Venti module logs to a child of the ``venti`` logger. `configure_logging`
attaches handlers to that logger only, so it never changes the logging of
the program importing Venti, and it is safe to call repeatedly.
"""

from __future__ import annotations

import logging
from pathlib import Path

LOGGER_NAME = "venti"
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Marks handlers added here, so repeated calls can recognise them.
_HANDLER_TAG = "_venti_handler"


def configure_logging(
    level: int | str | None = None,
    log_file: str | Path | None = None,
) -> logging.Logger:
    """Send Venti's log messages to the console and, optionally, a file.

    Called by the CLI, `run_workflow`, the workflow classes and the staging
    functions; call it yourself only when using lower-level functions
    directly.

    Parameters
    ----------
    level : int or str, optional
        Level for all Venti messages (e.g. ``"DEBUG"``). If omitted, a level
        already set on the ``venti`` logger is kept, otherwise ``INFO``.
    log_file : str or Path, optional
        Also append messages to this file (parent directories are created).
        Calling again with the same file adds no second handler.

    Returns
    -------
    logging.Logger
        The ``venti`` logger.

    Raises
    ------
    ValueError
        If `level` is not a valid logging level name.

    Notes
    -----
    The console handler stays silent whenever the root logger has a handler
    (for example after `logging.basicConfig`, even if called later), so
    messages are never printed twice.

    Examples
    --------
    >>> from venti.log_setup import configure_logging
    >>> configure_logging(log_file=None).name
    'venti'

    """
    logger = logging.getLogger(LOGGER_NAME)

    if level is not None:
        logger.setLevel(level.upper() if isinstance(level, str) else level)
    elif logger.level == logging.NOTSET:
        logger.setLevel(logging.INFO)

    formatter = logging.Formatter(LOG_FORMAT)
    if not any(getattr(h, _HANDLER_TAG, None) == "console" for h in logger.handlers):
        console = logging.StreamHandler()
        console.setFormatter(formatter)
        # Checked per message, not once: code imported later (e.g. the
        # staging scripts) may call logging.basicConfig, and records then
        # reach the root handler by propagation instead.
        console.addFilter(lambda _record: not logging.getLogger().handlers)
        setattr(console, _HANDLER_TAG, "console")
        logger.addHandler(console)

    if log_file is not None:
        path = Path(log_file).resolve()
        if not any(
            getattr(h, _HANDLER_TAG, None) == str(path) for h in logger.handlers
        ):
            path.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(path)
            file_handler.setFormatter(formatter)
            setattr(file_handler, _HANDLER_TAG, str(path))
            logger.addHandler(file_handler)
            logger.info("Writing log to %s", path)

    return logger
