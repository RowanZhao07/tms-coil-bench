"""Fetch pinned public measurements; convert the MC-B70 NPZ to MATLAB format."""
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen
import numpy as np
from scipy.io import savemat

HERE = Path(__file__).resolve().parent
CACHE = HERE.parents[1] / "tmp" / "measured_coil_library"
DIPOLE = "https://raw.githubusercontent.com/simnibs/dipole-fit/63d8cba254e3ade29c3db58eca9f5397a23e07c4/coildata/"
LIBRARY = "https://raw.githubusercontent.com/simnibs/simnibs/fde509f6f8de38a4b03947944cf7ee93b8bebabe/simnibs/resources/coil_models/Drakaki_BrainStim_2022/"


def download(url):
    try:
        with urlopen(url, timeout=45) as response:
            return response.read()
    except (URLError, TimeoutError, ConnectionError):
        owner, repo, revision, filename = url.removeprefix(
            "https://raw.githubusercontent.com/"
        ).split("/", 3)
        mirror = f"https://cdn.jsdelivr.net/gh/{owner}/{repo}@{revision}/{filename}"
        with urlopen(mirror, timeout=45) as response:
            return response.read()


def checked(path, url, expected, git_blob=False):
    data = path.read_bytes() if path.exists() else download(url)
    actual = (hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
              if git_blob else hashlib.sha256(data).hexdigest())
    if actual.lower() != expected.lower():
        raise ValueError(f"Unexpected contents in {path.name}: {actual}")
    if not path.exists():
        path.write_bytes(data)


def fetch_waveform(row):
    name = row["filename"]
    checked(CACHE / name, LIBRARY + name, row["git_blob"], git_blob=True)


if __name__ == "__main__":
    CACHE.mkdir(parents=True, exist_ok=True)
    name = "MagVenture_MC-B70_dipole-fit"
    for ext, sha in [("npz", "01A90ADA7AD344D8C78CAE4927148DAA137623DBD5489532893E437CFC4CFC77"),
                     ("stl", "0A9D60A006B039FBC7C42DA4515E279AC0CE22C92AFCD35F6681E6392E1F7B5D")]:
        checked(CACHE / f"{name}.{ext}", DIPOLE + f"{name}.{ext}", sha)
    with np.load(CACHE / f"{name}.npz", allow_pickle=False) as data:
        assert data["xyz"].shape == data["data"].shape == (95773,3)
        savemat(CACHE / "mc_b70_scan.mat", {"xyz_mm": data["xyz"], "B_T_per_A": data["data"]})
    with (HERE / "waveform_manifest.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 25
    with ThreadPoolExecutor(max_workers=5) as pool:
        list(pool.map(fetch_waveform, rows))
    print("Verified measured MC-B70 scan, casing mesh and 25 waveform files.")
