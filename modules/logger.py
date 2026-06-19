import os
import sys
import logging
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
from contextlib import contextmanager

import numpy as np

# Silence and supress warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pyfiglet")
warnings.filterwarnings("ignore", category=UserWarning, module="pkg_resources")

import pyfiglet
import pkg_resources


class Logger:
    """  
    Logger class for output file logging.
    """
    # Class wide variables for formatting
    INDENT_UNIT: int = 4 # single source of truth for indentation :)


    def __init__(
        self,
        name: str = "cluster_gen",
        log_file: Optional[str] = "cluster_gen.log",
        level: int = logging.INFO,
        file_mode: str = "a",
        include_timestamp: bool = False,
        include_level: bool = False,
        ):
        """   
        Args:
            name (str): Name of the logger.
            log_file (str, optional): Path to the log file. If None, messages are printed to the console instead. Defaults to "cluster_gen.log".
            level (int): Logging level (e.g., logging.INFO, logging.DEBUG). Defaults to logging.INFO.
            file_mode (str): Mode for opening the log file ('a' for append, 'w' for write). Defaults to 'a'.
        """
        self._logger = logging.getLogger(name) 
        self._logger.setLevel(level)
        self._logger.handlers.clear() # Clear existing handlers
        self._logger.propagate = False

        if include_timestamp and include_level:
            fmt = "[%(asctime)s] [%(levelname)-8s] %(message)s"
        elif include_timestamp:
            fmt = "[%(asctime)s] %(message)s"
        elif include_level:
            fmt = "[%(levelname)-8s] %(message)s"
        else:
            fmt = "%(message)s"

        formatter = logging.Formatter(fmt=fmt, datefmt="%Y-%m-%d %H:%M:%S")

        # File Handler
        if log_file:
            path = Path(log_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            fh = logging.FileHandler(path, mode=file_mode, encoding="utf-8")
            fh.setFormatter(formatter)
            self._logger.addHandler(fh)
        else:
            # No log file given - fall back to printing everything to the console
            ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(formatter)
            self._logger.addHandler(ch)

        
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

    # ======================== Utilities for Step Logging summaries and so on ==============================
    def step(self, current: int, total: int, title: str) -> None:
        """  
        Makes steps like [2/5] Fragmentation of Molecules --> always un-intendend
        """
        self._logger.info(f"[{current}/{total}] {title}")

    def substep(self, msg: str, level: int = 1) -> None:
        """    
        Detail line nested under the most recent step()/section()
        level is the nesting depth: intendation = INDENT_UNIT * level spaces
        """
        self._logger.info(f"{' ' * (self.INDENT_UNIT * level)}{msg}")

    def parameters(
        self,
        title: str,
        items: Union[Dict[str, Any], List[Tuple[str, Any]]],
        indent: int = 1,
    ) -> None:
        """
        Aligned 'label : value' block, "items" is a dict (insertion order preserved)
        or a list of (label, value) pairs when duplicate labels or conditional
        ordering is needed.
        """
        pairs = list(items.items()) if isinstance(items, dict) else list(items)
        if not pairs:
            return
        prefix = " " * (self.INDENT_UNIT * indent)
        label_width = max(len(str(label)) for label, _ in pairs)
        self._logger.info(f"{title}:")
        for label, value in pairs:
            self._logger.info(f"{prefix}{str(label):<{label_width}} : {value}")

    def worker_summary(
        self,
        title: str,
        per_worker: Dict[int, Any],
        indent: int = 1,
    ) -> None:
        """
        Log aggregate min/mean/max across workers instead of a per-worker dump.

        Two shapes are supported:
          - per_worker: Dict[int, float]            -> one summary line
          - per_worker: Dict[int, Dict[str, float]]  -> one line per inner key
            (e.g. one line per operator, averaged across workers)
        """
        if not per_worker:
            return
        n_workers = len(per_worker)
        prefix = " " * (self.INDENT_UNIT * indent)
        first = next(iter(per_worker.values()))

        if isinstance(first, dict):
            keys = list(first.keys())
            label_width = max(len(str(k)) for k in keys)
            self._logger.info(f"{title} ({n_workers} workers):")
            for key in keys:
                vals = [w.get(key, 0.0) for w in per_worker.values()]
                self._logger.info(
                    f"{prefix}{str(key):<{label_width}} : "
                    f"mean={np.mean(vals):5.1f}%  min={np.min(vals):5.1f}%  max={np.max(vals):5.1f}%"
                )
        else:
            vals = list(per_worker.values())
            self._logger.info(
                f"{title} ({n_workers} workers): "
                f"mean={np.mean(vals):.2f}  min={np.min(vals):.2f}  max={np.max(vals):.2f}"
            )

    # ======================== Sections & Separators ==============================
    def program_header(self) -> None:
        """ 
        Writes a Custom Header with the Program Name and Metadata
        """
        self._logger.info("=" * 80)
        ascii_title = pyfiglet.figlet_format("ML-BHOP", font="slant")
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        explanation = "A Machine Learning assisted Basin Hopping Optimization for Molecular Cluster Generation"
        programmer1 = "Made by: Lukas Meinschad"
        assistance = "Assisted by Jonas Schlagin, Klaus R. Liedl at University of Innsbruck, Liedl Lab"
        self._logger.info(ascii_title)
        self._logger.info(f"{explanation}")
        self._logger.info(f"{timestamp}")
        self._logger.info(f"{programmer1}")
        self._logger.info(f"{assistance}")
        self._logger.info("=" * 80)

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
            append: bool = False,
            sort_by_energy: bool = False
        ) -> None:
        """
        Write a list of Molecuel objects to an XYZ file

        Args:
            molecules (List[Any]): List of Molecule objects to write.
            filepath (str): Path to the output XYZ file.
            energies (List[float], optional): List of energies corresponding to each molecule. Defaults to None.
            comments (List[str], optional): List of comments for each molecule. Defaults to None.
            append (bool): Whether to append to the file if it exists. Defaults to False.
            sort_by_energy (bool): Whether to sort the molecules by energy before writing. Defaults to False.
        """
        if not molecules:
            self.warning("No molecules provided to write to XYZ file.")
            return
        
        # Optional energy sorting
        if sort_by_energy and energies:
            if len(energies) != len(molecules):
                self.warning("Length of energies does not match number of molecules. Skipping energy sorting.")
            
            combined = []
            for idx, mol in enumerate(molecules):
                eng = energies[idx] if energies and idx < len(energies) else None
                comm = comments[idx] if comments and idx < len(comments) else None
                combined.append((mol, eng, comm))
            
            combined.sort(key=lambda x: x[1] if x[1] is not None else float('inf'))
            # Unpack and sort back
            molecules = [item[0] for item in combined]
            energies = [item[1] for item in combined] if energies else None
            comments = [item[2] for item in combined] if comments else None

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

                if energies and idx < len(energies) and energies[idx] != float('inf'):
                    comment_parts.append(f"Energy: {energies[idx]:.6f} Hartree")

                comment_line = " | ".join(comment_parts)
                f.write(f"{n_atoms}\n")
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


if __name__ == "__main__":
    # Example usage of the Logger class
    logger = Logger(level=logging.INFO, include_timestamp=True, include_level=True)
    
    logger.program_header()
    logger.step(1, 3, "Starting the process")
    logger.substep("Loading data", level=1)
    logger.parameters("Configuration Parameters", {"param1": 10, "param2": 20}, indent=2)
    logger.worker_summary("Worker Performance", {0: {"accuracy": 95.5, "time": 120}, 1: {"accuracy": 96.0, "time": 110}}, indent=1)
    logger.write_program_header()





