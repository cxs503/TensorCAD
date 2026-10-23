"""TensorCAD: deterministic, agent-callable engineering geometry."""

from .hull import HullParameters, HullModel, generate_hull, validate_hull

__all__ = ["HullParameters", "HullModel", "generate_hull", "validate_hull"]
__version__ = "0.1.0"
