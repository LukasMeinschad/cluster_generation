import os
from datetime import datetime
from typing import List, Tuple, Optional, Any, TextIO
from enum import Enum
from contextlib import contextmanager
from pathlib import Path
import json


class LogLevel(Enum):
    """Log levels for message severity"""
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    CRITICAL = "CRITICAL"


class Formatter:
    """Handles formatting of log messages"""
    
    @staticmethod
    def header_line(width: int = 50, char: str = "=") -> str:
        """Create header line"""
        return char * width
    
    @staticmethod
    def separator_line(width: int = 50, char: str = "-") -> str:
        """Create separator line"""
        return char * width
    
    @staticmethod
    def timestamp(fmt: str = "%Y-%m-%d %H:%M:%S") -> str:
        """Get formatted timestamp"""
        return datetime.now().strftime(fmt)
    
    @staticmethod
    def section_header(title: str, width: int = 50) -> str:
        """Create formatted section header"""
        lines = [
            Formatter.separator_line(width),
            f"{title} - {Formatter.timestamp()}",
            Formatter.separator_line(width)
        ]
        return "\n".join(lines)
    
    @staticmethod
    def molecule_coords(atom_labels: List[str], coordinates: List[Tuple[float, float, float]]) -> str:
        """Format molecule coordinates"""
        lines = []
        for label, coord in zip(atom_labels, coordinates):
            lines.append(f"{label:4s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}")
        return "\n".join(lines)
    
    @staticmethod
    def bond_info(bonds: List[Tuple], bond_type: str) -> str:
        """Format bond information"""
        lines = [f"{bond_type} Bonds ({len(bonds)}):"]
        for bond in bonds:
            atom1 = bond.atom1
            atom2 = bond.atom2
            strength = bond.strength
            lines.append(f"  {atom1:4s} - {atom2:4s} (Strength: {strength:.3f})") 
        return "\n".join(lines)
    
    @staticmethod
    def remove_digits(text: str) -> str:
        """Remove all digits from a string"""
        return ''.join(char for char in text if not char.isdigit())


class FileWriter:
    """Handles file writing operations"""
    
    def __init__(self, filepath: Path, mode: str = "a"):
        """
        Initialize FileWriter
        
        Args:
            filepath: Path to output file
            mode: File mode ('a' for append, 'w' for write)
        """
        self.filepath = Path(filepath)
        self.mode = mode
    
    @contextmanager
    def open_file(self):
        """Context manager for safe file operations"""
        f = None
        try:
            self.filepath.parent.mkdir(parents=True, exist_ok=True)
            f = open(self.filepath, self.mode, encoding='utf-8')
            yield f
        except IOError as e:
            raise IOError(f"Failed to open file {self.filepath}: {e}") from e
        finally:
            if f:
                f.close()
    
    def write(self, content: str) -> None:
        """Write content to file"""
        with self.open_file() as f:
            f.write(content)
            if not content.endswith("\n"):
                f.write("\n")
    
    def write_lines(self, lines: List[str]) -> None:
        """Write multiple lines to file"""
        with self.open_file() as f:
            for line in lines:
                f.write(line)
                if not line.endswith("\n"):
                    f.write("\n")


class Logger:
    """
    Enhanced Logger class for cluster generation with support for
    different log levels, formatted output, and multiple output types.
    """
    
    def __init__(
        self, 
        out_file: str = "cluster_gen.out",
        mode: str = "a",
        min_level: LogLevel = LogLevel.INFO,
        line_width: int = 80
    ):
        """
        Initialize Logger
        
        Args:
            out_file: Name of the output file
            mode: File mode ('a' for append, 'w' for write)
            min_level: Minimum log level to write
            line_width: Width of separator lines
        """
        self.writer = FileWriter(out_file, mode)
        self.formatter = Formatter()
        self.min_level = min_level
        self.line_width = line_width
        self._log_count = 0
    
    # ==================== Core Logging Methods ====================
    
    def log(self, message: str, level: LogLevel = LogLevel.INFO) -> None:
        """
        Write a log message with specified level
        
        Args:
            message: Message to log
            level: Log level
        """
        if self._should_log(level):
            formatted = self._format_log_message(message, level)
            self.writer.write(formatted)
            self._log_count += 1
    
    def debug(self, message: str) -> None:
        """Write debug message"""
        self.log(message, LogLevel.DEBUG)
    
    def info(self, message: str) -> None:
        """Write info message"""
        self.log(message, LogLevel.INFO)
    
    def warning(self, message: str) -> None:
        """Write warning message"""
        self.log(message, LogLevel.WARNING)
    
    def error(self, message: str) -> None:
        """Write error message"""
        self.log(message, LogLevel.ERROR)
    
    def critical(self, message: str) -> None:
        """Write critical message"""
        self.log(message, LogLevel.CRITICAL)
    
    # ==================== Header and Sections ====================
    
    def write_header(self, title: str = "Cluster Generation Log") -> None:
        """Write general header to output file"""
        header = [
            self.formatter.header_line(self.line_width, "="),
            f"{title} - {self.formatter.timestamp()}",
            self.formatter.header_line(self.line_width, "="),
            ""
        ]
        self.writer.write_lines(header)
    
    def write_section(self, title: str, content: str) -> None:
        """
        Write a formatted section with title and content
        
        Args:
            title: Section title
            content: Section content
        """
        section = [
            self.formatter.separator_line(self.line_width),
            f"{title} - {self.formatter.timestamp()}",
            self.formatter.separator_line(self.line_width),
            content,
            ""
        ]
        self.writer.write_lines(section)
    
    # ==================== Molecule Logging ====================
    
    def write_molecule_info(self, molecule: Any) -> None:
        """
        Write comprehensive molecule information
        
        Args:
            molecule: Molecule object with atom_labels, coordinates, etc.
        """
        lines = [
            f"Molecule Name: {getattr(molecule, 'name', 'Unknown')}",
            f"Number of Atoms: {len(molecule.atom_labels)}",
            f"Charge: {getattr(molecule, 'charge', 0)}",
            f"Spin Multiplicity: {getattr(molecule, 'spin_mult', 1)}",
            "",
            "Atom Labels and Coordinates:",
            self.formatter.molecule_coords(molecule.atom_labels, molecule.coordinates),
            ""
        ]
        
        # Add bond information if available
        try:
            covalent_bonds, hydrogen_bonds = molecule.get_bonds()

            lines.extend([
                self.formatter.bond_info(covalent_bonds, "Covalent"),
                "",
                self.formatter.bond_info(hydrogen_bonds, "Hydrogen"),
            ])
        except Exception as e:
            self.debug(f"Could not retrieve bond information: {e}")
        
        self.write_section("Molecule Information", "\n".join(lines))

    def write_submolecule_info(self, submolecules: List[Any]) -> None: 
        """   
        Writes information about the submolecules, i.e the fragments obtained by connectivity analysis.

        We log the Number of Submolecules and for each submolecule the index, number of atoms, atom labels and coordinates.
        """
        lines = [f"Number of Submolecules: {len(submolecules)}"]
        for idx, submol in enumerate(submolecules, 1):
            lines.extend([
                f"Submolecule {idx}:",
                f"  Number of Atoms: {len(submol.atom_labels)}",
                "  Atom Labels and Coordinates:",
                *[f"    {line}" for line in 
                  self.formatter.molecule_coords(
                      submol.atom_labels, 
                      submol.coordinates
                  ).split('\n')],
                ""
            ])

        self.write_section(
            "Submolecule Information",
            "\n".join(lines)
        )

    def write_hbond_configurations(self, configurations: List[Any]) -> None:
        """
        Write hydrogen bond configurations
        
        Args:
            configurations: List of HBondConfiguration or SubMolecule objects
        """
        lines = []
        for idx, config in enumerate(configurations, 1):
            lines.extend([
                f"Configuration {idx}:",
                f"  Type: {type(config).__name__}",
            ])
            
            # Add specific information based on type
            if hasattr(config, 'donor_idx'):
                lines.append(f"  Donor Index: {config.donor_idx}")
            if hasattr(config, 'hydrogen_idx'):
                lines.append(f"  Hydrogen Index: {config.hydrogen_idx}")
            if hasattr(config, 'acceptor_idx'):
                lines.append(f"  Acceptor Index: {config.acceptor_idx}")
            if hasattr(config, 'angle'):
                lines.append(f"  Angle: {config.angle:.2f}°")
            
            # Add atom coordinates
            if hasattr(config, 'atom_labels') and hasattr(config, 'coordinates'):
                lines.extend([
                    "  Atom Labels and Coordinates:",
                    *[f"    {line}" for line in 
                      self.formatter.molecule_coords(
                          config.atom_labels, 
                          config.coordinates
                      ).split('\n')]
                ])
            
            lines.append("")
        
        self.write_section(
            f"Hydrogen Bond Configurations ({len(configurations)})",
            "\n".join(lines)
        )

    # =================== Transformation Logging ====================

    def write_reference_frame_info(self, ref_frame: Any) -> None:
        """  
        Writes information about the reference frame used for centering the molecule.
        """
        lines = [
            "Reference Frame Information:",
            f"Origin: {ref_frame.origin}",
            f"Eigenvalues of Inertia Tensor: {ref_frame.eigenvalues}",
            f"X-axis: {ref_frame.x_axis}",
            f"Y-axis: {ref_frame.y_axis}",
            f"Z-axis: {ref_frame.z_axis}",
        ]
        self.write_section(
            "Reference Frame Details",
            "\n".join(lines)
        )
    
    # ==================== Calculation Results ====================
    
    def write_scf_batch_result(
        self, 
        results: List[Tuple[Any, float]]
    ) -> None:
        """
        Write SCF batch calculation results
        
        Args:
            results: List of (molecule, energy) tuples
        """
        lines = []
        for idx, (mol, energy) in enumerate(results, 1):
            mol_name = getattr(mol, 'name', f'Molecule {idx}')
            lines.extend([
                f"Configuration {idx}: {mol_name}",
                f"  SCF Energy: {energy:.8f} Hartree",
                f"  SCF Energy: {energy * 27.211386245988:.6f} eV",
                f"  Number of Atoms: {len(mol.atom_labels)}",
                ""
            ])
        
        self.write_section(
            f"SCF Batch Results ({len(results)} configurations)",
            "\n".join(lines)
        )
    
    def write_optimization_results(
        self,
        results: List[Tuple[Any, float]],
        include_geometries: bool = True
    ) -> None:
        """
        Write optimization results with optional geometry details
        
        Args:
            results: List of (optimized_molecule, energy) tuples
            include_geometries: Whether to include full geometries
        """
        lines = []
        for idx, (mol, energy) in enumerate(results, 1):
            mol_name = getattr(mol, 'name', f'Optimized {idx}')
            lines.extend([
                f"Optimization {idx}: {mol_name}",
                f"  Final Energy: {energy:.8f} Hartree",
                f"  Number of Atoms: {len(mol.atom_labels)}",
            ])
            
            if include_geometries:
                lines.extend([
                    "  Optimized Geometry:",
                    *[f"    {line}" for line in 
                      self.formatter.molecule_coords(
                          mol.atom_labels,
                          mol.coordinates
                      ).split('\n')]
                ])
            
            lines.append("")
        
        self.write_section(
            f"Optimization Results ({len(results)} structures)",
            "\n".join(lines)
        )
    
    # ==================== File Export Methods ====================
    
    @staticmethod
    def write_trajectory_sampling(
        sampled_mols: List[Any],
        filename: str = "sampled_configurations_trj.xyz"
    ) -> None:
        """
        Write XYZ trajectory file of sampled molecules
        
        Args:
            sampled_mols: List of molecule objects
            filename: Output filename
        """
        filepath = Path(filename)
        
        try:
            with open(filepath, "w", encoding='utf-8') as trj_file:
                for i, mol in enumerate(sampled_mols, 1):
                    trj_file.write(f"{len(mol.atom_labels)}\n")
                    trj_file.write(f"Sampled Configuration {i}\n")
                    
                    for label, coord in zip(mol.atom_labels, mol.coordinates):
                        element = Formatter.remove_digits(label)
                        trj_file.write(
                            f"{element:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n"
                        )
        except IOError as e:
            raise IOError(f"Failed to write trajectory to {filepath}: {e}") from e
    
    @staticmethod
    def write_optimization_xyz_batch(
        opt_results: List[Tuple[Any, float]],
        output_dir: str = "optimized_geometries"
    ) -> None:
        """
        Write optimized geometries to individual XYZ files
        
        Args:
            opt_results: List of (optimized_molecule, energy) tuples
            output_dir: Output directory name
        """
        dirpath = Path(output_dir)
        dirpath.mkdir(parents=True, exist_ok=True)
        
        for idx, (opt_mol, opt_energy) in enumerate(opt_results, 1):
            file_path = dirpath / f"optimized_molecule_{idx:03d}.xyz"
            
            try:
                with open(file_path, "w", encoding='utf-8') as xyz_file:
                    xyz_file.write(f"{len(opt_mol.atom_labels)}\n")
                    xyz_file.write(
                        f"Optimized Molecule {idx} | "
                        f"Energy: {opt_energy:.8f} Hartree\n"
                    )
                    
                    for label, coord in zip(opt_mol.atom_labels, opt_mol.coordinates):
                        element = Formatter.remove_digits(label)
                        xyz_file.write(
                            f"{element:2s} {coord[0]:12.6f} {coord[1]:12.6f} {coord[2]:12.6f}\n"
                        )
            except IOError as e:
                raise IOError(f"Failed to write to {file_path}: {e}") from e
    
    @staticmethod
    def write_json_summary(
        results: List[Tuple[Any, float]],
        filename: str = "calculation_summary.json"
    ) -> None:
        """
        Write calculation summary as JSON
        
        Args:
            results: List of (molecule, energy) tuples
            filename: Output JSON filename
        """
        filepath = Path(filename)
        
        summary = {
            "timestamp": datetime.now().isoformat(),
            "num_structures": len(results),
            "structures": []
        }
        
        for idx, (mol, energy) in enumerate(results, 1):
            struct_data = {
                "id": idx,
                "name": getattr(mol, 'name', f'Structure_{idx}'),
                "energy_hartree": energy,
                "energy_ev": energy * 27.211386245988,
                "num_atoms": len(mol.atom_labels),
                "charge": getattr(mol, 'charge', 0),
                "spin_mult": getattr(mol, 'spin_mult', 1),
            }
            summary["structures"].append(struct_data)
        
        try:
            with open(filepath, "w", encoding='utf-8') as json_file:
                json.dump(summary, json_file, indent=2)
        except IOError as e:
            raise IOError(f"Failed to write JSON to {filepath}: {e}") from e
    
    # ==================== Private Helper Methods ====================
    
    def _should_log(self, level: LogLevel) -> bool:
        """Check if message should be logged based on level"""
        level_order = {
            LogLevel.DEBUG: 0,
            LogLevel.INFO: 1,
            LogLevel.WARNING: 2,
            LogLevel.ERROR: 3,
            LogLevel.CRITICAL: 4
        }
        return level_order[level] >= level_order[self.min_level]
    
    def _format_log_message(self, message: str, level: LogLevel) -> str:
        """Format log message with timestamp and level"""
        timestamp = self.formatter.timestamp()
        return f"[{timestamp}] [{level.value}] {message}\n"
    
    # ==================== Context Manager Support ====================
    
    def __enter__(self):
        """Enter context manager"""
        self.write_header()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Exit context manager"""
        if exc_type:
            self.error(f"Exception occurred: {exc_val}")
        
        footer = [
            "",
            self.formatter.header_line(self.line_width, "="),
            f"Log completed - {self.formatter.timestamp()}",
            f"Total log entries: {self._log_count}",
            self.formatter.header_line(self.line_width, "=")
        ]
        self.writer.write_lines(footer)
        
        return False  # Don't suppress exceptions


# ==================== Convenience Functions ====================

def create_logger(
    out_file: str = "cluster_gen.out",
    mode: str = "a",
    level: LogLevel = LogLevel.INFO
) -> Logger:
    """
    Factory function to create a Logger instance
    
    Args:
        out_file: Output file path
        mode: File mode
        level: Minimum log level
    
    Returns:
        Configured Logger instance
    """
    return Logger(out_file=out_file, mode=mode, min_level=level)





