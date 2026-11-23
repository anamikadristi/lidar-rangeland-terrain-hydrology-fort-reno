"""Shared settings and small helpers for the Fort Reno LiDAR terrain and hydrology workflow.

Author: Anamika Dristi, Ph.D.

All paths are relative to the repository root. Large raw downloads go to data/raw/
(git-ignored) unless the environment variable LRTH_RAW_DIR points somewhere else.
Full-resolution 1 m working rasters go to data/work/ (also git-ignored).
"""
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("LRTH_RAW_DIR", ROOT / "data" / "raw"))
WORK = Path(os.environ.get("LRTH_WORK_DIR", ROOT / "data" / "work"))
DATA = ROOT / "data"
FIG = ROOT / "outputs" / "figures"
TAB = ROOT / "outputs" / "tables"
RAS = ROOT / "outputs" / "rasters"
VEC = ROOT / "outputs" / "vectors"
for _d in (RAW / "laz", RAW / "dem", RAW / "nhd", RAW / "metadata", WORK, FIG, TAB, RAS, VEC):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "python-requests (public data download script)"

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------
TNM_API = "https://tnmaccess.nationalmap.gov/api/v1/products"
LPC_PROJECT = "OK_Statewide_D22"
LPC_WORK_UNIT = "OK_Statewide_2_D22"
LPC_BASE = ("https://rockyweb.usgs.gov/vdelivery/Datasets/Staged/Elevation/LPC/Projects/"
            f"{LPC_PROJECT}/{LPC_WORK_UNIT}/LAZ/")
# Six 5,000 ft by 5,000 ft tiles (2 columns by 3 rows). Tile names give the
# south-west corner in Oklahoma North State Plane US survey feet.
TILES = [f"{x}_{y}" for y in (190000, 195000, 200000) for x in (1955000, 1960000)]
LAZ_NAMES = [f"USGS_LPC_{LPC_PROJECT}_{t}.laz" for t in TILES]

DEM_URL = ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/1m/Projects/"
           f"{LPC_PROJECT}/TIFF/USGS_1M_14_x58y394_{LPC_PROJECT}.tif")
DEM_NAME = DEM_URL.rsplit("/", 1)[1]

NHD_URL = ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/GPKG/"
           "NHD_H_11100301_HU8_GPKG.zip")
NHD_GPKG = "NHD_H_11100301_HU8_GPKG.gpkg"

META_BASE = ("https://prd-tnm.s3.amazonaws.com/StagedProducts/Elevation/metadata/"
             f"{LPC_PROJECT}/{LPC_WORK_UNIT}/reports/")
META_FILES = ["USGS_OK_Statewide_2_D22_Summary_Report.pdf",
              "USGS_OK_Statewide_2_D22_FINAL_LPC_Report.txt",
              "vendor_provided_xml/OK_Statewide_2_D22_Classified_LiDAR_Point_Cloud_Metadata.xml",
              "vendor_provided_xml/OK_Statewide_2_D22_DEM_Metadata.xml"]

# ---------------------------------------------------------------------------
# Coordinate reference systems
# ---------------------------------------------------------------------------
SRC_EPSG = 6553    # NAD83(2011) / Oklahoma North (ftUS), the point cloud CRS
DST_EPSG = 26914   # NAD83 / UTM zone 14N (m), the CRS of the USGS 1 m DEM
CRS_LABEL = "NAD83 / UTM zone 14N (EPSG:26914), heights NAVD88 (GEOID18), m"
FTUS_TO_M = 1200.0 / 3937.0
CELL = 1.0         # working cell size, m

# LAS classification codes used by this project (ASPRS LAS 1.4, USGS Lidar Base Spec)
CLASS_NAMES = {
    1: "Processed, but unclassified",
    2: "Bare earth (ground)",
    7: "Low noise",
    9: "Water",
    17: "Bridge deck",
    18: "High noise",
    20: "Ignored ground (near a breakline)",
    21: "Snow",
    22: "Temporal exclusion",
}
NOISE_CLASSES = (7, 18)


# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------
def load_grid():
    """Return the working grid definition written by 02_point_cloud_qa.py."""
    import json
    from affine import Affine
    with open(DATA / "grid_definition.json") as f:
        g = json.load(f)
    g["transform"] = Affine(CELL, 0, g["xmin"], 0, -CELL, g["ymax"])
    return g


def read_raster(path, masked=True):
    import rasterio
    with rasterio.open(path) as src:
        a = src.read(1)
        if masked and src.nodata is not None:
            a = np.where(a == src.nodata, np.nan, a) if a.dtype.kind == "f" else a
        return a


def write_raster(path, arr, transform=None, nodata=None, cog=False, epsg=DST_EPSG, scale=None, offset=None):
    """Write a single-band GeoTIFF (DEFLATE). With cog=True write a Cloud-Optimized GeoTIFF."""
    import rasterio
    if transform is None:
        transform = load_grid()["transform"]
    arr = np.asarray(arr)
    if arr.dtype == np.float64:
        arr = arr.astype("float32")
    if arr.dtype == bool:
        arr = arr.astype("uint8")
    if nodata is None and arr.dtype.kind == "f":
        nodata = -9999.0
    if arr.dtype.kind == "f":
        arr = np.where(np.isfinite(arr), arr, nodata).astype(arr.dtype)
    profile = dict(driver="COG" if cog else "GTiff", height=arr.shape[0], width=arr.shape[1],
                   count=1, dtype=arr.dtype, crs=f"EPSG:{epsg}", transform=transform,
                   nodata=nodata, compress="deflate")
    if arr.dtype.kind == "f":
        profile["predictor"] = 3 if not cog else "FLOATING_POINT"
    if not cog:
        profile.update(tiled=True, blockxsize=512, blockysize=512)
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr, 1)
        if scale is not None:
            dst.scales = (scale,)
            dst.offsets = (offset,)
            dst.units = ("metre",)


def block_reduce(a, f, func=np.nanmean):
    """Aggregate a 2-D array by an integer factor f (edges are trimmed)."""
    ny, nx = (a.shape[0] // f) * f, (a.shape[1] // f) * f
    b = a[:ny, :nx].reshape(ny // f, f, nx // f, f)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return func(func(b, axis=3), axis=1)


def hillshade(z, cell=CELL, azimuth=315.0, altitude=45.0, zfac=1.0):
    """Hillshade (0 to 1) with the light from the given azimuth (degrees clockwise from north)."""
    dzdy, dzdx = np.gradient(z * zfac, cell)
    dzdy = -dzdy  # row index grows to the south, so flip the sign to get d/d(north)
    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(-dzdx, -dzdy)  # downslope direction, clockwise from north
    az, alt = np.radians(azimuth), np.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(hs, 0, 1)


# ---------------------------------------------------------------------------
# Map figure helpers
# ---------------------------------------------------------------------------
def map_extent(g=None):
    g = g or load_grid()
    return [g["xmin"], g["xmax"], g["ymin"], g["ymax"]]


def finish_map(ax, g=None, scalebar_m=500, note=True, color="k", box=True):
    """Add a scale bar, a north arrow and tidy axes. All maps are north-up in UTM 14N."""
    g = g or load_grid()
    x0 = g["xmin"] + 0.05 * (g["xmax"] - g["xmin"])
    y0 = g["ymin"] + 0.035 * (g["ymax"] - g["ymin"])
    bbox = dict(facecolor="white", alpha=0.75, edgecolor="none", pad=1.5) if box else None
    ax.plot([x0, x0 + scalebar_m], [y0, y0], color=color, lw=3, solid_capstyle="butt", zorder=20)
    ax.text(x0 + scalebar_m / 2, y0 + 45, f"{scalebar_m} m", ha="center", va="bottom",
            fontsize=8, color=color, bbox=bbox, zorder=20)
    xa = g["xmax"] - 0.07 * (g["xmax"] - g["xmin"])
    ya = g["ymax"] - 0.11 * (g["ymax"] - g["ymin"])
    ax.annotate("N", xy=(xa, ya + 260), xytext=(xa, ya), ha="center", va="center",
                fontsize=10, fontweight="bold", color=color, bbox=bbox, zorder=20,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.8))
    ax.set_xlim(g["xmin"], g["xmax"])
    ax.set_ylim(g["ymin"], g["ymax"])
    ax.set_aspect("equal")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=7)
    for lab in ax.get_yticklabels():
        lab.set_rotation(90)
        lab.set_va("center")
    if note:
        ax.set_xlabel("Easting (m)", fontsize=8)
        ax.set_ylabel("Northing (m)", fontsize=8)


def crs_note(fig, extra=""):
    fig.text(0.5, 0.005, f"North is up. CRS: {CRS_LABEL}. {extra}".strip(),
             ha="center", va="bottom", fontsize=7.5, color="0.25")


# ---------------------------------------------------------------------------
# Point cloud reading and gridding
# ---------------------------------------------------------------------------
def read_tile(path):
    """Read one LAZ tile and return a dict of numpy arrays in the working CRS (UTM 14N, m).

    The source tiles are in Oklahoma North State Plane with US survey feet for both the
    horizontal and the vertical axis. I reproject x and y with pyproj and convert z to m.
    """
    import laspy
    from pyproj import Transformer
    las = laspy.read(path)
    tr = Transformer.from_crs(SRC_EPSG, DST_EPSG, always_xy=True)
    x, y = tr.transform(np.asarray(las.x), np.asarray(las.y))
    return dict(
        x=x, y=y, z=np.asarray(las.z) * FTUS_TO_M,
        cls=np.asarray(las.classification).astype(np.uint8),
        rn=np.asarray(las.return_number).astype(np.uint8),
        nr=np.asarray(las.number_of_returns).astype(np.uint8),
        overlap=np.asarray(las.overlap).astype(bool),
        withheld=np.asarray(las.withheld).astype(bool),
        synthetic=np.asarray(las.synthetic).astype(bool),
        key_point=np.asarray(las.key_point).astype(bool),
        scan_angle=np.asarray(las.scan_angle) * 0.006,  # LAS 1.4 stores 0.006 degree units
        psid=np.asarray(las.point_source_id),
        gps=np.asarray(las.gps_time),
        intensity=np.asarray(las.intensity),
        header=las.header,
    )


def cell_index(x, y, g):
    """Flat cell index of each point in the working grid, and a mask of points inside it."""
    col = np.floor((x - g["xmin"]) / CELL).astype(np.int64)
    row = np.floor((g["ymax"] - y) / CELL).astype(np.int64)
    ok = (col >= 0) & (col < g["nx"]) & (row >= 0) & (row < g["ny"])
    return row * g["nx"] + col, ok


def grid_count(idx, g):
    return np.bincount(idx, minlength=g["nx"] * g["ny"]).reshape(g["ny"], g["nx"])


def grid_sum(idx, v, g):
    return np.bincount(idx, weights=v, minlength=g["nx"] * g["ny"]).reshape(g["ny"], g["nx"])


def grid_extreme(idx, v, g, kind="max"):
    """Per-cell maximum or minimum of v. Cells with no point are NaN."""
    out = np.full(g["nx"] * g["ny"], np.nan, dtype="float32")
    if idx.size:
        order = np.argsort(idx, kind="stable")
        si, sv = idx[order], v[order]
        starts = np.flatnonzero(np.r_[True, si[1:] != si[:-1]])
        red = np.maximum.reduceat(sv, starts) if kind == "max" else np.minimum.reduceat(sv, starts)
        out[si[starts]] = red
    return out.reshape(g["ny"], g["nx"])


# ---------------------------------------------------------------------------
# Hydrologic conditioning
# ---------------------------------------------------------------------------
# pysheds 0.5 still calls numpy.in1d, which newer numpy versions removed.
if not hasattr(np, "in1d"):
    np.in1d = np.isin


import heapq  # noqa: E402

import numba  # noqa: E402


@numba.njit(cache=True)
def _pf(dem, eps):
    ny, nx = dem.shape
    out = dem.copy()
    closed = np.zeros((ny, nx), dtype=np.bool_)
    heap = [(out[0, 0], 0)]
    closed[0, 0] = True
    for c in range(nx):
        for r in (0, ny - 1):
            if not closed[r, c]:
                closed[r, c] = True
                heapq.heappush(heap, (out[r, c], r * nx + c))
    for r in range(ny):
        for c in (0, nx - 1):
            if not closed[r, c]:
                closed[r, c] = True
                heapq.heappush(heap, (out[r, c], r * nx + c))
    while len(heap) > 0:
        z, k = heapq.heappop(heap)
        r0, c0 = k // nx, k % nx
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                r, c = r0 + dr, c0 + dc
                if r < 0 or c < 0 or r >= ny or c >= nx or closed[r, c]:
                    continue
                closed[r, c] = True
                if out[r, c] <= z:
                    out[r, c] = z + eps
                heapq.heappush(heap, (out[r, c], r * nx + c))
    return out


def priority_flood(dem, eps=0.0):
    """Fill closed depressions with the priority-flood method (Barnes and others, 2014).

    The grid edge is treated as open. Cells are visited from the lowest edge cell upward,
    and any cell lower than the cell it is reached from is raised to that level. With
    eps > 0 each raised cell is set eps higher than its parent, so filled areas keep a
    very small gradient toward their spill point and every cell has a downhill neighbour.
    """
    return _pf(np.ascontiguousarray(dem, dtype=np.float64), float(eps))


@numba.njit(cache=True)
def _breach(z, pr, pc, max_len, radius):
    ny, nx = z.shape
    r0, r1 = max(pr - radius, 0), min(pr + radius + 1, ny)
    c0, c1 = max(pc - radius, 0), min(pc + radius + 1, nx)
    h, w = r1 - r0, c1 - c0
    zp = z[pr, pc]
    cost = np.full((h, w), np.inf)
    dist = np.zeros((h, w))
    prev = np.full((h, w), -1, dtype=np.int64)
    done = np.zeros((h, w), dtype=np.bool_)
    s = (pr - r0) * w + (pc - c0)
    cost[pr - r0, pc - c0] = 0.0
    heap = [(0.0, s)]
    end = -1
    while len(heap) > 0:
        cst, k = heapq.heappop(heap)
        r, c = k // w, k % w
        if done[r, c]:
            continue
        done[r, c] = True
        if z[r + r0, c + c0] < zp - 0.01:
            end = k
            break
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                rr, cc = r + dr, c + dc
                if rr < 0 or cc < 0 or rr >= h or cc >= w or done[rr, cc]:
                    continue
                step = 1.41421356 if (dr != 0 and dc != 0) else 1.0
                d = dist[r, c] + step
                if d > max_len:
                    continue
                cut = z[rr + r0, cc + c0] - zp
                if cut < 0.0:
                    cut = 0.0
                nc = cst + cut * step + 0.01 * step
                if nc < cost[rr, cc]:
                    cost[rr, cc] = nc
                    dist[rr, cc] = d
                    prev[rr, cc] = k
                    heapq.heappush(heap, (nc, rr * w + cc))
    if end < 0:
        return False, 0.0, 0.0, -1, -1
    er, ec = end // w, end % w
    zt = z[er + r0, ec + c0]
    total = dist[er, ec]
    max_cut = 0.0
    k = end
    while k != s:
        r, c = k // w, k % w
        target = zp - (zp - zt) * dist[r, c] / total
        cut = z[r + r0, c + c0] - target
        if cut > 0.0:
            if cut > max_cut:
                max_cut = cut
            z[r + r0, c + c0] = target
        k = prev[r, c]
    return True, total, max_cut, er + r0, ec + c0


def breach_path(z, pit_r, pit_c, max_len=150.0, radius=200):
    """Find the cheapest cut from a depression bottom to lower ground and carve it into z.

    This is a small least-cost breach in the spirit of Lindsay (2016). From the lowest cell
    of a depression I run Dijkstra's search. The cost of stepping onto a cell is the height
    that would have to be cut away there to keep a channel at the pit level, plus a small
    cost per metre. The search stops at the first cell lower than the pit. The cells on
    the path are then lowered to a straight ramp from the pit down to that cell.

    Returns (found, path_length_m, max_cut_m, end_row, end_col). z is changed in place.
    """
    return _breach(z, int(pit_r), int(pit_c), float(max_len), int(radius))


# ---------------------------------------------------------------------------
# Terrain and vector helpers
# ---------------------------------------------------------------------------
def horn_gradients(z, cell):
    """dz/dx (east) and dz/dy (north) with Horn's weights. Edges use nearest-value padding."""
    p = np.pad(z, 1, mode="edge")
    a, b, c = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]
    d, f = p[1:-1, :-2], p[1:-1, 2:]
    g_, h, i = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]
    dzdx = ((c + 2 * f + i) - (a + 2 * d + g_)) / (8 * cell)
    dzdy = ((a + 2 * b + c) - (g_ + 2 * h + i)) / (8 * cell)   # north minus south
    return dzdx, dzdy


def nhd_layer(layer, clip=True):
    """Read one NHD layer for the working area, reprojected to UTM 14N (2-D geometry)."""
    import warnings
    import geopandas as gpd
    import shapely
    g = load_grid()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gdf = gpd.read_file(RAW / "nhd" / NHD_GPKG, layer=layer,
                            bbox=(g["lon_min"] - 0.01, g["lat_min"] - 0.01, g["lon_max"] + 0.01, g["lat_max"] + 0.01))
    gdf = gdf.to_crs(DST_EPSG)
    gdf["geometry"] = shapely.force_2d(gdf.geometry.values)
    if clip:
        gdf = gpd.clip(gdf, shapely.box(g["xmin"], g["ymin"], g["xmax"], g["ymax"]))
        gdf = gdf[~gdf.geometry.is_empty].reset_index(drop=True)
    return gdf


def rasterize(geoms, g=None, all_touched=False, values=None):
    """Burn geometries into the working grid. Returns uint8 (or int32 when values are given)."""
    from rasterio import features
    g = g or load_grid()
    geoms = list(geoms)
    if not geoms:
        return np.zeros((g["ny"], g["nx"]), np.uint8)
    shapes = geoms if values is None else list(zip(geoms, values))
    return features.rasterize(shapes, out_shape=(g["ny"], g["nx"]), transform=g["transform"],
                              all_touched=all_touched, fill=0,
                              dtype="uint8" if values is None else "int32")


def error_stats(d):
    """Standard elevation-difference statistics for an array of differences (m)."""
    d = np.asarray(d, dtype="float64")
    d = d[np.isfinite(d)]
    if d.size == 0:
        return dict(n=0)
    med = np.median(d)
    return dict(n=int(d.size), bias_m=d.mean(), median_m=med, std_m=d.std(),
                rmse_m=np.sqrt((d ** 2).mean()), mae_m=np.abs(d).mean(),
                nmad_m=1.4826 * np.median(np.abs(d - med)),
                p95_abs_m=np.percentile(np.abs(d), 95),
                within_5cm_pct=100 * (np.abs(d) <= 0.05).mean(),
                within_10cm_pct=100 * (np.abs(d) <= 0.10).mean(),
                within_19p6cm_pct=100 * (np.abs(d) <= 0.196).mean())


@numba.njit(cache=True)
def label_catchments(fdir, out_r, out_c):
    """Label every cell with the number (1, 2, ...) of the outlet cell its D8 flow path reaches.

    fdir uses the pysheds / ESRI codes (64 N, 128 NE, 1 E, 2 SE, 4 S, 8 SW, 16 W, 32 NW).
    Cells whose path ends somewhere else get 0. I wrote this because the pysheds catchment
    function missed diagonal inflow to outlet cells on the grid edge.
    """
    ny, nx = fdir.shape
    lab = np.zeros((ny, nx), dtype=np.int32)      # 0 unvisited, -1 no outlet, k outlet number
    for k in range(out_r.size):
        lab[out_r[k], out_c[k]] = k + 1
    path = np.empty(ny * nx, dtype=np.int64)
    for r0 in range(ny):
        for c0 in range(nx):
            if lab[r0, c0] != 0:
                continue
            n = 0
            r, c = r0, c0
            val = -1
            while True:
                if lab[r, c] != 0:
                    val = lab[r, c]
                    break
                path[n] = r * nx + c
                n += 1
                lab[r, c] = -2                     # on the current path (guards against loops)
                d = fdir[r, c]
                if d == 64:
                    r -= 1
                elif d == 128:
                    r -= 1
                    c += 1
                elif d == 1:
                    c += 1
                elif d == 2:
                    r += 1
                    c += 1
                elif d == 4:
                    r += 1
                elif d == 8:
                    r += 1
                    c -= 1
                elif d == 16:
                    c -= 1
                elif d == 32:
                    r -= 1
                    c -= 1
                else:
                    break
                if r < 0 or c < 0 or r >= ny or c >= nx:
                    break
                if lab[r, c] == -2:
                    break
            for i in range(n):
                lab[path[i] // nx, path[i] % nx] = val
    for r in range(ny):
        for c in range(nx):
            if lab[r, c] < 0:
                lab[r, c] = 0
    return lab


def shrink_png(path, limit_mb=1.4):
    """If a figure is larger than limit_mb, store it with a 256-colour palette to save space."""
    from pathlib import Path
    from PIL import Image
    path = Path(path)
    if path.stat().st_size <= limit_mb * 1e6:
        return
    dpi = Image.open(path).info.get("dpi", (150, 150))
    Image.open(path).convert("RGB").quantize(colors=256, method=Image.Quantize.MEDIANCUT,
                                             dither=Image.Dither.NONE).save(path, optimize=True, dpi=dpi)
