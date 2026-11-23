"""Step 1. Download the public inputs.

I download:
  - six USGS 3DEP LiDAR point cloud tiles (LAZ) from the OK_Statewide_D22 project,
  - the matching USGS 3DEP 1 m DEM tile (used only as a comparison product),
  - the NHD high-resolution GeoPackage for hydrologic unit 11100301,
  - the USGS validation report and LAS QC report for the LiDAR work unit.

The National Map API is queried first so the tile list, sizes and publication dates are
saved in data/tnm_product_index.csv. If the API does not answer, the known tile URLs in
common.py are used. Nothing here needs an account or a key.
"""
import csv
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor

import requests

import common as C

HEAD = {"User-Agent": C.USER_AGENT}
BBOX = "-98.044,35.523,-98.013,35.562"   # lon/lat box inside the six tiles


def fetch(url, dest, min_bytes=1000):
    if dest.exists() and dest.stat().st_size > min_bytes:
        print(f"  have {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")
        return
    print(f"  downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(30):  # resume with a Range request if the connection drops
        done = tmp.stat().st_size if tmp.exists() else 0
        head = dict(HEAD, Range=f"bytes={done}-") if done else HEAD
        try:
            with requests.get(url, headers=head, stream=True, timeout=120) as r:
                if r.status_code == 416:
                    break
                r.raise_for_status()
                mode = "ab" if (done and r.status_code == 206) else "wb"
                with open(tmp, mode) as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            break
        except requests.RequestException as e:
            print(f"  retry {attempt + 1} for {dest.name}: {type(e).__name__}")
            time.sleep(30)  # the server answers 429 when there are too many connections
    else:
        raise RuntimeError(f"could not download {url}")
    tmp.rename(dest)
    print(f"  saved {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)")


def query_tnm(dataset):
    params = {"datasets": dataset, "bbox": BBOX, "max": 100, "outputFormat": "JSON"}
    r = requests.get(C.TNM_API, params=params, headers=HEAD, timeout=90)
    r.raise_for_status()
    return r.json().get("items", [])


def main():
    rows = []
    try:
        for ds in ("Lidar Point Cloud (LPC)", "Digital Elevation Model (DEM) 1 meter"):
            for it in query_tnm(ds):
                name = it["downloadURL"].rsplit("/", 1)[1]
                used = name in C.LAZ_NAMES or name == C.DEM_NAME
                bb = it["boundingBox"]
                rows.append([ds, it["title"], it["publicationDate"], it["sizeInBytes"],
                             bb["minX"], bb["minY"], bb["maxX"], bb["maxY"],
                             it["downloadURL"], it.get("metaUrl", ""), "yes" if used else "no"])
        with open(C.DATA / "tnm_product_index.csv", "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["dataset", "title", "publication_date", "size_bytes", "min_lon", "min_lat",
                        "max_lon", "max_lat", "download_url", "metadata_url", "used_here"])
            w.writerows(rows)
        print(f"TNM API returned {len(rows)} products; index saved to data/tnm_product_index.csv")
    except Exception as e:  # the API is sometimes slow; the direct URLs still work
        print(f"TNM API query failed ({e}); using the known tile URLs.", file=sys.stderr)

    print("LiDAR point cloud tiles:")
    # The USGS LiDAR server is slow per connection, so I fetch three tiles at a time.
    # More than about three connections at once gets HTTP 429 (too many requests).
    with ThreadPoolExecutor(max_workers=3) as ex:
        list(ex.map(lambda n: fetch(C.LPC_BASE + n, C.RAW / "laz" / n), C.LAZ_NAMES))
    print("USGS 1 m DEM:")
    fetch(C.DEM_URL, C.RAW / "dem" / C.DEM_NAME)
    print("NHD:")
    z = C.RAW / "nhd" / C.NHD_URL.rsplit("/", 1)[1]
    fetch(C.NHD_URL, z)
    if not (C.RAW / "nhd" / C.NHD_GPKG).exists():
        with zipfile.ZipFile(z) as zf:
            zf.extractall(C.RAW / "nhd")
    print("Project reports:")
    for name in C.META_FILES:
        fetch(C.META_BASE + name, C.RAW / "metadata" / name.rsplit("/", 1)[-1])
    print("Done.")


if __name__ == "__main__":
    main()
