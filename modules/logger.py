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

    # ======================== General Logging Functionalities ==============================
    def write_xyz_trajectory(
            self,
            molecules: List[Any],
            filepath: str,
            energies: Optional[List[float]] = None,
            comments: Optional[List[str]] = None,
            append: bool = False
        ) -> None:
        """
        Write a list of Molecuel objects to an XYZ file

        Args:
            molecules (List[Any]): List of Molecule objects to write.
            filepath (str): Path to the output XYZ file.
            energies (List[float], optional): List of energies corresponding to each molecule. Defaults to None.
            comments (List[str], optional): List of comments for each molecule. Defaults to None.
            append (bool): Whether to append to the file if it exists. Defaults to False.
        """
        if not molecules:
            self.warning("No molecules provided to write to XYZ file.")
            return
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        mode = "a" if append else "w"
        n_written = 0
        with open(path, mode, encoding="utf-8") as f:
            for idx, mol in enumerate(molecules):
                labels = mol.atom_labels
                labels = [self.remove_digits(label) for label in labels]
                coords = mol.coordinates
                n_atoms = len(labels)

                comment_parts = []
                if comments and idx < len(comments):
                    comment_parts.append(comments[idx])
                else:
                    comment_parts.append(f"Structure {idx+1}")
                if energies and idx < len(energies):
                    comment_parts.append(f"Energy: {energies[idx]:.6f} Hartree")
                comment_line = " | ".join(comment_parts)
                f.write(f"{n_atoms}\n)")
                f.write(f"{comment_line}\n")
                for i in range(n_atoms):
                    x,y,z = coords[i]
                    f.write(f"{labels[i]} {x:.6f} {y:.6f} {z:.6f}\n")
                n_written += 1
        self.info(f"Wrote {n_written} structures to {filepath}")

    def write_single_xyz(
            self,
            molecule: Any,
            filepath: str,
            energy: Optional[float] = None,
            comment: Optional[str] = None,
            append: bool = False 
        ) -> None:
        """  
        Write a single molecule to an XYZ file

        Args:
            molecule (Any): Molecule object to write.
            filepath (str): Path to the output XYZ file.
            energy (float, optional): Energy of the molecule. Defaults to None.
            comment (str, optional): Comment for the molecule. Defaults to None.
            append (bool): Whether to append to the file if it exists. Defaults to False.
        """
        energies = [energy] if energy is not None else None
        comments = [comment] if comment is not None else None
        self.write_xyz_trajectory(
            molecules=[molecule],
            filepath=filepath,
            energies=energies,
            comments=comments,
            append=append
        )

    def log_matrix(
            self,
            matrix: "np.ndarray",
            row_labels: Optional[List[str]] = None,
            col_labels: Optional[List[str]] = None,
            title: str = "Matrix",
            fmt: str = ".4f",
            max_col_width: int = 14
        ) -> None:
        """  
        Log a 2D numpy array as a formatted table

        Args:
            matrix (np.ndarray): 2D array to log.
            row_labels (List[str], optional): Labels for rows. Defaults to None.
            col_labels (List[str], optional): Labels for columns. Defaults to None.
            title (str): Title for the matrix. Defaults to "Matrix".
            fmt (str): Format string for values (e.g., ".4f"). Defaults to ".4f".
            max_col_width (int): Maximum width for each column in characters. Defaults to 14.
        """
        import numpy as np
        if matrix.ndim != 2:
            self.warning(f"log_matrix expects a 2D array, but got {matrix.ndim}D array.")
            return
        n_rows, n_cols = matrix.shape 
        if row_labels is None:
            row_labels = [f"Row {i}" for i in range(n_rows)]
        if col_labels is None:
            col_labels = [f"Col {j}" for j in range(n_cols)]
        # Truncate labels to max_col_width
        col_labels = [l[:max_col_width] for l in col_labels]
        row_label_width = max(len(l) for l in row_labels) + 2
        self.section(title)
        # Header Line
        header = " " * row_label_width + " ".join(f"{c:{max_col_width}}" for c in col_labels)
        self._logger.info(header)
        self.separator(char="-", width=len(header))

        # Data Lines
        for i, row in enumerate(matrix):
            values = " ".join(f"{v:{max_col_width}{fmt}}" for v in row)
            self._logger.info(f"{row_labels[i]:<{row_label_width}}{values}")

        self.separator(char="=", width=len(header))


    @staticmethod
    def remove_digits(s: str) -> str:
        """  
        Utility function to remove digits from a string (used for cleaner logging)
        """
        return ''.join(filter(lambda x: not x.isdigit(), s))





