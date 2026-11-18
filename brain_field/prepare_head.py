"""Fetch and verify the native FieldTrip Subject01 head and anatomical parcellation."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import time
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
CACHE = HERE.parent / "tmp" / "brain_field" / "m2m_Subject01"
SOURCE = "https://download.fieldtriptoolbox.org/getting_started/simnibs/m2m_Subject01/"
MANIFEST = HERE / "head_inputs.json"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def fetch(filename, expected, direct=False):
    dest = CACHE / filename
    if dest.is_file() and dest.stat().st_size == expected["bytes"] and digest(dest) == expected["sha256"]:
        print(f"Verified {filename}", flush=True)
        return
    session = requests.Session()
    session.trust_env = not direct
    url = SOURCE + filename
    size = expected["bytes"]
    block_size = 1024 * 1024
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp = dest.with_name(dest.name + ".download")
    def part(start):
        end = min(start + block_size, size) - 1
        for attempt in range(4):
            try:
                with session.get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=(20, 90)) as response:
                    response.raise_for_status()
                    block = response.content
                    if response.status_code != 206 or len(block) != end - start + 1:
                        raise ValueError(f"Incomplete byte range for {filename}")
                    return start, block
            except (requests.RequestException, ValueError):
                if attempt == 3:
                    raise
                time.sleep(1 + attempt)
    with temp.open("wb") as stream, concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool:
        for start, data in pool.map(part, range(0, size, block_size)):
            stream.seek(start)
            stream.write(data)
    if digest(temp) != expected["sha256"]:
        raise ValueError(f"SHA-256 mismatch: {filename}")
    os.replace(temp, dest)
    print(f"Verified {filename}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direct", action="store_true", help="Connect without system HTTP proxies")
    args = parser.parse_args()
    for name, expected in json.loads(MANIFEST.read_text(encoding="utf-8")).items():
        fetch(name, expected, args.direct)
