from dataclasses import dataclass
from pathlib import Path
from typing import Any

from matplotlib.figure import Figure


@dataclass(frozen=True)
class FigureConfig:
    savefig: bool = False
    save_dir: str | Path | None = None
    prefix: str | None = None
    file_extension: str = "png"

    def __post_init__(self) -> None:
        if self.savefig and self.save_dir is not None:
            Path(self.save_dir).mkdir(parents=True, exist_ok=True)

    def get_path(self, figure_name: str) -> Path | None:
        if not self.savefig:
            return None

        if self.save_dir is None:
            raise ValueError("save_dir must be specified when savefig=True.")

        if self.prefix is None:
            raise ValueError("prefix must be specified when savefig=True.")

        if not figure_name:
            raise ValueError("figure_name must not be empty.")

        extension = self.file_extension.lstrip(".")
        prefix = "" if self.prefix == "" else f"{self.prefix}_"

        return Path(self.save_dir) / f"{prefix}{figure_name}.{extension}"


def save_figure(
    fig: Figure,
    config: FigureConfig,
    figure_name: str,
    **kwarg_savefig: Any,
) -> Path | None:
    path = config.get_path(figure_name)

    if path is None:
        return None

    fig.savefig(path, **kwarg_savefig)
    return path


__all__ = ["FigureConfig", "save_figure"]
