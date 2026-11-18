"""Fetch the five CHARM heads from the versioned SimNIBS OSF group dataset.

Only the four requested members per person are transferred. OSF metadata pins
the release, ZIP CRC32 verifies the selected streams, and the member manifest
pins each decompressed byte sequence. Run under Python 3.11+ with requests.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import tempfile
import io
import json
from pathlib import Path
import struct
import threading
import time
import zipfile
import zlib

import requests


HERE = Path(__file__).resolve().parent
CACHE = HERE.parents[1] / "tmp" / "brain_cohort" / "heads"
MANIFEST = HERE / "head_sources.json"
DOWNLOAD = "https://osf.io/download/5ysve/?version=1"
LENGTH = 2702041298
SUBJECTS = ("sub01", "sub09", "sub10", "sub12", "sub15")


class Archive:
    def __init__(self, direct):
        self.direct = direct
        self.threads = threading.local()
        self.sign_lock = threading.Lock()
        self.signed_at = 0.0
        self.end_start = LENGTH - 131072
        tail = self.read(self.end_start, LENGTH - 1)
        # The ZIP reader rebases the central-directory offsets to this suffix.
        self.index = zipfile.ZipFile(io.BytesIO(tail))

    def session(self):
        if not hasattr(self.threads, "session"):
            client = requests.Session()
            client.trust_env = not self.direct
            self.threads.session = client
        return self.threads.session

    def signed_url(self):
        # OSF supplies a pre-signed object URL valid for roughly one minute.
        with self.sign_lock:
            if time.monotonic() - self.signed_at > 32 or not self.signed_at:
                params = {"timestamp": str(time.time_ns())}
                with self.session().get(
                    DOWNLOAD,
                    params=params,
                    headers={"Range": "bytes=0-0"},
                    timeout=(15, 45),
                    stream=True,
                ) as response:
                    response.raise_for_status()
                    self.url = response.url
                    self.signed_at = time.monotonic()
            return self.url

    def read(self, start, end):
        for attempt in range(5):
            try:
                with self.session().get(
                    self.signed_url(),
                    headers={"Range": f"bytes={start}-{end}"},
                    timeout=(18, 90),
                ) as resp:
                    resp.raise_for_status()
                    if (
                        resp.status_code != 206
                        or resp.headers.get("Content-Range")
                        != f"bytes {start}-{end}/{LENGTH}"
                    ):
                        raise OSError("Unexpected HTTP byte range")
                    data = resp.content
                    if len(data) != end - start + 1:
                        raise OSError("Short HTTP byte range")
                    return data
            except (OSError, requests.RequestException):
                if attempt == 4:
                    raise
                time.sleep(attempt + 1)

    def extract(self, name, path, expected):
        info = self.index.getinfo(name)
        start = self.end_start + info.header_offset
        header = self.read(start, start + 29)
        if header[:4] != b"PK\x03\x04" or info.compress_type != zipfile.ZIP_DEFLATED:
            raise ValueError(f"Unexpected zip header {name}")
        name_length, extra_length = struct.unpack_from("<HH", header, 26)
        first = start + 30 + name_length + extra_length
        step = 1024 * 1024
        ranges = [
            (i, min(first + info.compress_size - 1, i + step - 1))
            for i in range(first, first + info.compress_size, step)
        ]
        path.parent.mkdir(parents=True, exist_ok=True)
        # Each writer has its own file until both ZIP and manifest digests pass.
        handle = tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=path.name + ".", suffix=".part", delete=False
        )
        temp = Path(handle.name)
        decoder = zlib.decompressobj(-15)
        crc = 0
        sha = hashlib.sha256()

        def emit(file, part):
            nonlocal crc
            file.write(part)
            crc = zlib.crc32(part, crc)
            sha.update(part)

        try:
            with handle as out, ThreadPoolExecutor(max_workers=16) as pool:
                for block in pool.map(lambda item: self.read(*item), ranges):
                    emit(out, decoder.decompress(block))
                emit(out, decoder.flush())
            if (
                not decoder.eof
                or decoder.unused_data
                or temp.stat().st_size != info.file_size
                or crc != info.CRC
            ):
                raise ValueError(f"Corrupt ZIP member {name}")
            result = {
                "path": name,
                "bytes": info.file_size,
                "crc32": f"{crc:08x}",
                "sha256": sha.hexdigest(),
            }
            if result != expected:
                raise ValueError(f"Input differs from the member manifest: {name}")
            temp.replace(path)
        finally:
            temp.unlink(missing_ok=True)


def matches(path, expected):
    if not path.is_file() or path.stat().st_size != expected["bytes"]:
        return False
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest() == expected["sha256"]


def prepare(subjects, direct):
    with MANIFEST.open(encoding="utf-8") as stream:
        manifest = json.load(stream)["subjects"]
    archive = None
    for subject in subjects:
        for filename, expected in manifest[subject].items():
            path = CACHE / subject / filename
            if not matches(path, expected):
                if archive is None:
                    archive = Archive(direct)
                archive.extract(expected["path"], path, expected)
            print(subject, filename, "verified", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--direct",
        action="store_true",
        help="Connect without environment proxy settings",
    )
    parser.add_argument(
        "subjects",
        nargs="*",
        help="Subject IDs (default: all five published participants)",
    )
    args = parser.parse_args()
    subjects = tuple(args.subjects) or SUBJECTS
    if any(s not in SUBJECTS for s in subjects):
        parser.error("Unknown subject")
    prepare(subjects, args.direct)
