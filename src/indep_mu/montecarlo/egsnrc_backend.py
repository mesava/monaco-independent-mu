from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shutil

from .interface import MonteCarloCapabilities


def _resolve_executable(value: str | Path) -> Path | None:
    path = Path(value)
    if path.parent != Path(".") or path.is_absolute():
        return path if path.is_file() else None

    resolved = shutil.which(str(value))
    return Path(resolved) if resolved else None


@dataclass(frozen=True)
class EgsnrcReferenceBackend:
    """Environment descriptor for the CPU EGSnrc reference transport.

    This class deliberately does not launch simulations yet.  It establishes
    an auditable boundary between "EGSnrc is available" and "a validated
    machine transport input has been built".
    """

    beamnrc_executable: str | Path = "beamnrc"
    dosxyznrc_executable: str | Path = "dosxyznrc"
    hen_house: Path | None = None
    egs_home: Path | None = None

    @property
    def capabilities(self) -> MonteCarloCapabilities:
        return MonteCarloCapabilities(
            backend_id="egsnrc-reference",
            supports_dynamic_mlc=True,
            supports_dynamic_jaws=True,
            supports_dynamic_gantry=True,
            supports_patient_ct=True,
            supports_dose_to_medium=True,
            supports_gpu=False,
        )

    def resolved_beamnrc(self) -> Path | None:
        return _resolve_executable(self.beamnrc_executable)

    def resolved_dosxyznrc(self) -> Path | None:
        return _resolve_executable(self.dosxyznrc_executable)

    def validate_environment(self) -> tuple[str, ...]:
        errors: list[str] = []

        beamnrc = self.resolved_beamnrc()
        if beamnrc is None:
            errors.append(
                f"BEAMnrc executable not found: {self.beamnrc_executable!s}"
            )

        dosxyz = self.resolved_dosxyznrc()
        if dosxyz is None:
            errors.append(
                f"DOSXYZnrc executable not found: {self.dosxyznrc_executable!s}"
            )

        for label, path in (
            ("HEN_HOUSE", self.hen_house),
            ("EGS_HOME", self.egs_home),
        ):
            if path is None:
                continue
            if not Path(path).is_dir():
                errors.append(f"{label} directory does not exist: {path}")

        # On POSIX an existing direct path should also be executable.  PATH
        # resolution via shutil.which already applies the platform rule.
        if os.name != "nt":
            for label, path in (
                ("BEAMnrc", beamnrc),
                ("DOSXYZnrc", dosxyz),
            ):
                if path is not None and not os.access(path, os.X_OK):
                    errors.append(f"{label} is not executable: {path}")

        return tuple(errors)
