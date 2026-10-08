import logging
import os
from importlib import resources
from pathlib import Path


def resolve_database_path(
    db_argument: str = "", logger: logging.Logger | None = None
) -> str:
    """
    Determine which database directory to use.

    Preference order:
    1. Explicit --db argument
    2. CCTYPER_DB environment variable
    3. Packaged data shipped with the cctyper distribution
    """

    if db_argument:
        return db_argument

    env_db = os.environ.get("CCTYPER_DB")
    if env_db:
        return env_db

    try:
        data_root = resources.files(__package__).joinpath("data")
    except Exception as exc:
        raise RuntimeError("Could not locate packaged database files") from exc

    if not data_root.exists():
        raise RuntimeError(
            f"Packaged database directory missing at {data_root}"
        )

    if logger:
        logger.info("Using packaged database at %s", data_root)

    return str(Path(data_root))
