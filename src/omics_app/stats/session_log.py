"""
Shared helper for appending to store-analysis-log -- the Dash
equivalent of R's rv$analysis_log (initialized to
["Application started successfully"], R line 1846, then appended to
throughout the server from upload, comparisons, analysis, enrichment,
and color-mapping actions).

Since Dash requires every callback writing to the same Output to use
allow_duplicate=True and there's no single place all those actions
funnel through, this small helper is imported by each of those 5
UI modules so every append formats its entry the same way rather than
duplicating the format string five times.

Timestamps on each entry are an addition beyond R's plain message list
(R just appends bare strings, e.g. line 1974's
`sprintf("Data loaded: %d rows, %d columns", ...)` with no time
attached) -- a reasonable enhancement for a log meant to be read after
the fact, not a literal port requirement.
"""

from datetime import datetime


def append_log_entry(existing_log: list[str] | None, message: str) -> list[str]:
    """Returns a NEW list with a timestamped entry appended -- never
    mutates existing_log in place, since Dash Store callbacks should
    treat State values as read-only snapshots."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    log = list(existing_log or [])
    log.append(f"{timestamp} - {message}")
    return log
