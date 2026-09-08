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
