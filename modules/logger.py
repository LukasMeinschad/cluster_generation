import os 
import sys
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from contextlib import contextmanager


class Logger:
    """  
    Logger class for output file logging.
    """
    def __init__(
        self,
        name: str = "cluster_gen",
        log_file: Optional[str] = "cluster_gen.log",
        level: int = logging.INFO,
        file_mode: str = "a"
        ):
        """   
        Args:
            name (str): Name of the logger.
            log_file (str, optional): Path to the log file. If None, logging to file is disabled. Defaults to "cluster_gen.log".
            level (int): Logging level (e.g., logging.INFO, logging.DEBUG). Defaults to logging.INFO.
            file_mode (str): Mode for opening the log file ('a' for append, 'w' for write). Defaults to 'a'.
        """
        self._logger = logging.getLogger(name) 
        self._logger.setLevel(level)
        self._logger.handlers.clear() # Clear existing handlers

        formatter = logging.Formatter(
            fmt="[%(asctime)s] [%(levelname)-8s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )

        # File Handler 
        if log_file:
            path = Path(log_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(path, mode=file_mode, encoding="utf-8")
            fh.setFormatter(formatter)
            self._logger.addHandler(fh)

        
    # ========================= Core Logging Methods =============================
    def debug(self, msg: str) -> None:
        """   
        Log a Debug Message
        """
        self._logger.debug(msg) 
    def info(self, msg: str) -> None:
        """  
        Log an info message
        """
        self._logger.info(msg)
    def warning(self, msg: str) -> None:
        """  
        Log a warning message
        """
        self._logger.warning(msg)
    def error(self, msg: str) -> None:
        """  
        Log an error message
        """
        self._logger.error(msg)
    def critical(self, msg: str) -> None:
        """  
        Log a critical message
        """
        self._logger.critical(msg)

    # ======================== Sections & Separators ==============================

    def separator(self, char: str = "-", width: int = 72) -> None:
        """   
        Write a visual separator line
        """
        self._logger.info(char * width)
    def header(self, title: str, width: int = 72) -> None:
        """   
        Writes a prominent section header
        """
        line = "=" * width
        self._logger.info(line)
        self._logger.info(f"{title.center(width)}")
        self._logger.info(line)
    def section(self, title: str) -> None:
        """   
        Writes a subsection header-
        """
        self.separator()
        self._logger.info(f"{title}")
        self.separator()
    def write_program_header(self) -> None:
        """  
        Writes the Standardized Big Header for the Program
        """
        self.header("CLUSTER GEN", width=80)
        # Metadata
        self._logger.info(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self._logger.info(f"Python Version: {sys.version.split()[0]}")
        self._logger.info(f"Working Directory: {os.getcwd()}")
        self._logger.info(f"Log File: {self._logger.handlers[0].baseFilename if self._logger.handlers else 'None'}")
        self._logger.info("Made by: Lukas Meinschad")






