from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable


@dataclass(frozen=True)
class MonteCarloCapabilities:
    backend_id: str
    supports_dynamic_mlc: bool
    supports_dynamic_jaws: bool
    supports_dynamic_gantry: bool
    supports_patient_ct: bool
    supports_dose_to_medium: bool
    supports_gpu: bool


@dataclass(frozen=True)
class MonteCarloRunConfig:
    histories: int
    random_seed: int
    output_directory: Path
    requested_statistical_uncertainty_percent: float | None = None

    def __post_init__(self) -> None:
        if self.histories <= 0:
            raise ValueError("histories must be positive.")
        if self.random_seed < 0:
            raise ValueError("random_seed must be non-negative.")
        if (
            self.requested_statistical_uncertainty_percent is not None
            and self.requested_statistical_uncertainty_percent <= 0
        ):
            raise ValueError(
                "requested_statistical_uncertainty_percent must be positive."
            )


@runtime_checkable
class MonteCarloBackend(Protocol):
    """Transport backend contract.

    EGSnrc/BEAMnrc will be the initial reference implementation. A future GPU
    backend must satisfy the same contract and be independently validated
    against the reference implementation.
    """

    @property
    def capabilities(self) -> MonteCarloCapabilities:
        ...

    def validate_environment(self) -> tuple[str, ...]:
        """Return environment/setup errors; an empty tuple means ready."""
        ...
