""" 
Input Reader Module for Cluster Generation Algorithm
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

def _coerce_value(raw_value: str) -> Any:
    """  
    A
    """
    raw_value = raw_value.strip()

    if raw_value == "None":
        return None
    if raw_value == "True":
        return True
    if raw_value == "False":
        return False
    if (
        (raw_value.startswith('"')) and (raw_value.endswith('"')) or
        (raw_value.startswith("'")) and (raw_value.endswith("'"))
    ):
        return raw_value[1:-1]
    try:
        return int(raw_value)
    except ValueError:
        pass
    try:
        return float(raw_value)
    except ValueError:
        pass

    return raw_value


def _normalize_phase_name(name: str) -> str:
    return "".join(ch.lower() for ch in name if ch.isalnum())


def _section_matches_phase(section_name: str, phase: str) -> bool:
    section_norm = _normalize_phase_name(section_name)
    phase_norm = _normalize_phase_name(phase)

    return phase_norm in section_norm or section_norm in phase_norm

def read_ml_bhop_input(input_path: str) -> dict[str, dict[str, Any]]:
    """  
    Reads the input file for the ML-BHOP algorithm and returns a structured dictionary
    """
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")
    parsed: dict[str, dict[str, Any]] = {}
    current_section: str | None = None

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("-"):
            # This is just the header section
            continue
        if line.endswith(":") and not line.startswith("+"):
            # This is a new section header
            current_section = line[:-1].strip()
            parsed.setdefault(current_section, {})
            continue

        if line.startswith("+"):
            line = line[1:].strip()

        if "=" not in line:
            continue

        key, value = [part.strip() for part in line.split("=", 1)]
        coerced_value = _coerce_value(value)

        if current_section is None:
            parsed.setdefault("root", {})[key] = coerced_value
        else:
            parsed[current_section][key] = coerced_value

    return parsed

def read_operator_distribution(input_path: str, phase: str | None = None) -> dict[str, tuple[str, float]]:
    """ 
    Reads the operator distribution from a given input path for the BHMC algorithm.
    The operator distribution file has different sections for different phases of the algorithm

    e.g NL-OpsTraining:
          + OperatorA = 0.5
    
    and so on. If a phase is specified, only the operators from the matching section will be returned.
    If no phase is specified, all operators will be returned in a dictionary {'PhaseName': [(OperatorName, Weight), ...], ...}
    """
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Operator distribution file not found: {input_path}")
    
    sections: dict[str, list[tuple[str, float]]] = {}
    current_section: str | None = None

    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()

        if not line:
            continue

        if line.endswith(":"):
            current_section = line[:-1].strip()
            sections.setdefault(current_section, [])
            continue

        if line.startswith("+"):
            line = line[1:].strip()

        if "=" not in line or current_section is None:
            continue

        key, value = [part.strip() for part in line.split("=", 1)]
        try:
            weight = float(value)
            sections[current_section].append((key, weight))
        except ValueError:
            continue
    if phase is not None:
        for section_name, operators in sections.items():
            if _section_matches_phase(section_name, phase):
                return dict(operators)
        raise ValueError(f"No matching section found for phase: {phase}")
    else:
        return {section: dict(operators) for section, operators in sections.items()}


if __name__ == "__main__":
    import pprint
    print("Testing ML-BHOP Input Reader")
    input_data = read_ml_bhop_input("/media/storage_6/lme/master_thesis/cluster_generation/ml_bhop.inp")
    pprint.pprint(input_data) 

    # Testing operator distribution reader
    print("\nTesting Operator Distribution Reader")
    operator_dist = read_operator_distribution("/media/storage_6/lme/master_thesis/cluster_generation/operator_dist.inp")
    pprint.pprint(operator_dist)
