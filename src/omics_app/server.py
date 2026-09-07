"""
Server-side cache, backing large/transient data that shouldn't round-trip
through client-side dcc.Store (which serializes to JSON in the browser on
every callback -- fine for small config values, much too slow for raw
uploaded file bytes).

Usage: app.py calls cache.init_app(app.server) once at startup. Any
module needing to stash something large imports `cache` from here
directly -- avoids a circular import with app.py.
"""

from flask_caching import Cache

cache = Cache(
    config={
        "CACHE_TYPE": "SimpleCache",  # in-memory, single-process -- fine for
        # local/dev use. Swap to "FileSystemCache"
        # or "RedisCache" before any multi-worker
        # deployment, since SimpleCache is not
        # shared across processes.
        "CACHE_DEFAULT_TIMEOUT": 3600,  # 1 hour; uploaded files older than
        # this get evicted automatically.
    }
)
