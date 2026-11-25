"""Shared HTTP client: polite rate limits, retries, and a test hook.

Tests replace :data:`TRANSPORT` with an ``httpx.MockTransport`` that serves recorded
fixtures, so the unit suite never touches the network.
"""

from __future__ import annotations

import time
from typing import Any

import httpx

from genor import __version__
from genor.errors import GenorError

TRANSPORT: httpx.BaseTransport | None = None
USER_AGENT = f"genor/{__version__} (+https://github.com/amanzadi/GeNor)"
_last: dict[str, float] = {}

# minimum seconds between requests per host
RATE = {"rest.ensembl.org": 0.07, "string-db.org": 1.0, "version-12-0.string-db.org": 1.0,
        "api.platform.opentargets.org": 0.1, "biit.cs.ut.ee": 0.5}


def client() -> httpx.Client:
    return httpx.Client(transport=TRANSPORT, timeout=60, headers={"User-Agent": USER_AGENT},
                        follow_redirects=True)


def _wait(host: str) -> None:
    gap = RATE.get(host, 0.2)
    dt = time.monotonic() - _last.get(host, 0)
    if dt < gap and TRANSPORT is None:
        time.sleep(gap - dt)
    _last[host] = time.monotonic()


def request(method: str, url: str, *, retries: int = 3, **kw: Any) -> Any:
    host = httpx.URL(url).host
    err: Exception | None = None
    for attempt in range(retries):
        _wait(host)
        try:
            with client() as c:
                r = c.request(method, url, **kw)
            if r.status_code == 429 or r.status_code >= 500:
                retry_after = float(r.headers.get("Retry-After", 2 ** attempt))
                if TRANSPORT is None:
                    time.sleep(min(retry_after, 30))
                err = GenorError(f"{host} returned {r.status_code}", code="remote_unavailable")
                continue
            if r.status_code >= 400:
                raise GenorError(f"{host} returned {r.status_code}: {r.text[:300]}",
                                 code="remote_bad_request")
            return r.json()
        except httpx.HTTPError as e:
            err = e
    raise GenorError(f"Request to {host} failed: {err}", code="remote_unavailable",
                     hint="The service may be down or have changed its API. Retry later.")
