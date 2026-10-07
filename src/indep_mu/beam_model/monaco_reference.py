from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


_TAG_PATTERN = re.compile(
    r"<(?P<name>[A-Za-z0-9_]+)>\s*"
    r"(?P<value>[-+0-9.eE]+)\s*"
    r"</(?P=name)>"
)


@dataclass(frozen=True)
class MonacoLeafModelReference:
    """Raw Monaco leaf/jaw model parameters.

    These values are preserved for comparison/audit only. They are TPS beam
    model parameters and are NOT automatically interpreted as literal physical
    dimensions, probabilities or a complete Monte Carlo head geometry.
    """

    source_name: str
    values: dict[str, float]

    def require(self, name: str) -> float:
        try:
            return self.values[name]
        except KeyError as exc:
            raise ValueError(
                f"Monaco reference parameter {name!r} is missing."
            ) from exc


def parse_monaco_leaf_model_text(
    text: str,
    *,
    source_name: str = "<memory>",
) -> MonacoLeafModelReference:
    values: dict[str, float] = {}
    for match in _TAG_PATTERN.finditer(text):
        name = match.group("name")
        if name in values:
            raise ValueError(f"Duplicate Monaco parameter {name!r}.")
        values[name] = float(match.group("value"))

    if not values:
        raise ValueError("No <Name>numeric</Name> parameters found.")

    return MonacoLeafModelReference(
        source_name=source_name,
        values=values,
    )


def load_monaco_leaf_model(path: str | Path) -> MonacoLeafModelReference:
    input_path = Path(path)
    return parse_monaco_leaf_model_text(
        input_path.read_text(encoding="utf-8"),
        source_name=input_path.name,
    )
