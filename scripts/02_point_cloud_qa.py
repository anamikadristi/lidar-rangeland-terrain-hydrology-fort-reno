"""Step 2. Point cloud quality checks and the working grid.

For each LAZ tile I count points by classification code, return number, number of
returns and flag (overlap, withheld, synthetic, key point). I read the acquisition dates from the GPS
times, and I map the point density on the 1 m working grid. The working grid is the
largest north-up rectangle in UTM zone 14N that fits inside the six State Plane tiles,
snapped to whole metres so it lines up with the USGS 1 m DEM pixels.

Outputs: data/grid_definition.json, outputs/tables/qa_*.csv,
outputs/figures/fig01_point_density_qa.png, and density rasters in data/work/.
"""
import datetime as dt
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyproj import Transformer

import common as C


def define_grid():
    """Largest UTM rectangle inside the tile block (tile names give SW corners in ftUS)."""
    xs = sorted({int(t.split("_")[0]) for t in C.TILES})
    ys = sorted({int(t.split("_")[1]) for t in C.TILES})
    x0, x1, y0, y1 = xs[0], xs[-1] + 5000, ys[0], ys[-1] + 5000
    tr = Transformer.from_crs(C.SRC_EPSG, C.DST_EPSG, always_xy=True)
    n = 200
    west = tr.transform(np.full(n, x0), np.linspace(y0, y1, n))[0].max()
    east = tr.transform(np.full(n, x1), np.linspace(y0, y1, n))[0].min()
    south = tr.transform(np.linspace(x0, x1, n), np.full(n, y0))[1].max()
    north = tr.transform(np.linspace(x0, x1, n), np.full(n, y1))[1].min()
    g = dict(xmin=float(np.ceil(west)), xmax=float(np.floor(east)),
             ymin=float(np.ceil(south)), ymax=float(np.floor(north)), cell=C.CELL, epsg=C.DST_EPSG)
    g["nx"] = int(g["xmax"] - g["xmin"])
    g["ny"] = int(g["ymax"] - g["ymin"])
    g["area_km2"] = g["nx"] * g["ny"] / 1e6
    ll = Transformer.from_crs(C.DST_EPSG, 4269, always_xy=True)
    g["lon_min"], g["lat_min"] = ll.transform(g["xmin"], g["ymin"])
    g["lon_max"], g["lat_max"] = ll.transform(g["xmax"], g["ymax"])
    with open(C.DATA / "grid_definition.json", "w") as f:
        json.dump(g, f, indent=1)
    return C.load_grid()


def gps_to_date(t):
    """Adjusted standard GPS time (seconds since 1980-01-06 minus 1e9) to a UTC date."""
    return (dt.datetime(1980, 1, 6) + dt.timedelta(seconds=float(t) + 1e9 - 18)).date()


def main():
    g = define_grid()
    print(f"Working grid: {g['nx']} x {g['ny']} cells of {C.CELL} m = {g['area_km2']:.2f} km2")
    shape = (g["ny"], g["nx"])
    n_all = np.zeros(shape, np.int32)
    n_first = np.zeros(shape, np.int32)
    n_ground = np.zeros(shape, np.int32)
    n_overlap = np.zeros(shape, np.int32)
    tile_rows, cls_rows, ret_rows = [], [], []
    psids = set()
    for name in C.LAZ_NAMES:
        p = C.read_tile(C.RAW / "laz" / name)
        h = p["header"]
        n = p["z"].size
        valid = ~p["withheld"] & ~np.isin(p["cls"], C.NOISE_CLASSES)
        days = sorted(set(np.floor((p["gps"] + 1e9 - 18) / 86400.0).tolist()))
        dates = sorted({gps_to_date(d * 86400.0 - 1e9 + 18) for d in days})
        lines_here = set(np.flatnonzero(np.bincount(p["psid"])).tolist())
        psids |= lines_here
        tile_rows.append(dict(
            tile=name, las_version=str(h.version), point_format=h.point_format.id, points=n,
            first_returns=int((p["rn"] == 1).sum()), ground_points=int((p["cls"] == 2).sum()),
            overlap_flag=int(p["overlap"].sum()), withheld_flag=int(p["withheld"].sum()),
            synthetic_flag=int(p["synthetic"].sum()), key_point_flag=int(p["key_point"].sum()),
            noise_points=int(np.isin(p["cls"], C.NOISE_CLASSES).sum()),
            scan_angle_min_deg=float(p["scan_angle"].min()), scan_angle_max_deg=float(p["scan_angle"].max()),
            flight_lines=len(lines_here),
            z_min_m=float(p["z"][valid].min()), z_max_m=float(p["z"][valid].max()),
            acquisition_dates_utc=" ".join(str(d) for d in dates)))
        for c, k in enumerate(np.bincount(p["cls"])):
            if k:
                cls_rows.append(dict(tile=name, code=int(c), points=int(k)))
        for r, k in enumerate(np.bincount(p["rn"])):
            if k:
                ret_rows.append(dict(tile=name, kind="return_number", value=int(r), points=int(k)))
        for r, k in enumerate(np.bincount(p["nr"])):
            if k:
                ret_rows.append(dict(tile=name, kind="number_of_returns", value=int(r), points=int(k)))
        idx, ok = C.cell_index(p["x"], p["y"], g)
        ok &= valid
        n_all += C.grid_count(idx[ok], g).astype(np.int32)
        n_first += C.grid_count(idx[ok & (p["rn"] == 1)], g).astype(np.int32)
        n_ground += C.grid_count(idx[ok & (p["cls"] == 2)], g).astype(np.int32)
        n_overlap += C.grid_count(idx[ok & p["overlap"]], g).astype(np.int32)
        print(f"  {name}: {n:,} points, dates {tile_rows[-1]['acquisition_dates_utc']}")

    tiles = pd.DataFrame(tile_rows)
    tiles.to_csv(C.TAB / "qa_tile_summary.csv", index=False)
    cls = pd.DataFrame(cls_rows).groupby("code", as_index=False)["points"].sum()
    cls["class_name"] = cls["code"].map(C.CLASS_NAMES).fillna("Other")
    cls["percent"] = 100 * cls["points"] / cls["points"].sum()
    cls[["code", "class_name", "points", "percent"]].to_csv(C.TAB / "qa_classification_counts.csv",
                                                           index=False, float_format="%.4f")
    ret = pd.DataFrame(ret_rows).groupby(["kind", "value"], as_index=False)["points"].sum()
    ret["percent"] = 100 * ret["points"] / ret.groupby("kind")["points"].transform("sum")
    ret.to_csv(C.TAB / "qa_return_counts.csv", index=False, float_format="%.4f")

    # Density statistics on the working grid (valid points only: not withheld, not noise)
    area = g["nx"] * g["ny"] * C.CELL ** 2
    d10_first = C.block_reduce(n_first.astype("float32"), 10, np.mean)
    d10_ground = C.block_reduce(n_ground.astype("float32"), 10, np.mean)
    dens = pd.DataFrame([
        dict(metric="Area of working grid (km2)", value=area / 1e6),
        dict(metric="Valid points inside grid", value=int(n_all.sum())),
        dict(metric="All-return density (points/m2)", value=n_all.sum() / area),
        dict(metric="First-return (pulse) density (points/m2)", value=n_first.sum() / area),
        dict(metric="Nominal pulse spacing from first returns (m)", value=(area / n_first.sum()) ** 0.5),
        dict(metric="Ground-point density (points/m2)", value=n_ground.sum() / area),
        dict(metric="Share of 1 m cells with no valid return (%)", value=100 * (n_all == 0).mean()),
        dict(metric="Share of 1 m cells with no ground point (%)", value=100 * (n_ground == 0).mean()),
        dict(metric="Share of 10 m cells with pulse density below 2 per m2 (%)", value=100 * (d10_first < 2).mean()),
        dict(metric="5th percentile of 10 m pulse density (points/m2)", value=float(np.percentile(d10_first, 5))),
        dict(metric="Median of 10 m pulse density (points/m2)", value=float(np.median(d10_first))),
        dict(metric="95th percentile of 10 m pulse density (points/m2)", value=float(np.percentile(d10_first, 95))),
        dict(metric="Share of valid points with the overlap flag (%)", value=100 * n_overlap.sum() / n_all.sum()),
        dict(metric="Flight lines (point source IDs) in the six tiles", value=len(psids)),
    ])
    dens.to_csv(C.TAB / "qa_density_summary.csv", index=False, float_format="%.4f")
    print(dens.to_string(index=False))

    C.write_raster(C.WORK / "density_all_1m.tif", n_all.astype("uint16"), g["transform"], nodata=None)
    C.write_raster(C.WORK / "density_first_1m.tif", n_first.astype("uint16"), g["transform"], nodata=None)
    C.write_raster(C.WORK / "density_ground_1m.tif", n_ground.astype("uint16"), g["transform"], nodata=None)

    # Figure
    fig = plt.figure(figsize=(12, 8.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.95], hspace=0.38, wspace=0.42,
                          left=0.05, right=0.98, top=0.93, bottom=0.09)
    ext = C.map_extent(g)
    ax = fig.add_subplot(gs[:, 0])
    im = ax.imshow(d10_first, extent=ext, cmap="viridis", vmin=0, vmax=np.percentile(d10_first, 99))
    C.finish_map(ax, g)
    ax.set_title("(a) Pulse density (first returns), 10 m cells", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.07, shrink=0.9, label="first returns per m$^2$")
    ax = fig.add_subplot(gs[:, 1])
    im = ax.imshow(d10_ground, extent=ext, cmap="cividis", vmin=0, vmax=np.percentile(d10_ground, 99))
    C.finish_map(ax, g)
    ax.set_ylabel("")
    ax.set_title("(b) Ground-point density, 10 m cells", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.07, shrink=0.9, label="ground points per m$^2$")
    ax = fig.add_subplot(gs[0, 2])
    lab = [f"{c} {C.CLASS_NAMES.get(c, 'Other').split(' (')[0].split(',')[0]}" for c in cls["code"]]
    ax.barh(lab, cls["points"] / 1e6, color="#33658a")
    ax.set_xscale("log")
    ax.set_xlabel("points (millions, log scale)")
    ax.set_title("(c) Classification codes", fontsize=10)
    ax.invert_yaxis()
    ax.tick_params(labelsize=8)
    ax = fig.add_subplot(gs[1, 2])
    rn = ret[ret["kind"] == "return_number"]
    ax.bar(rn["value"].astype(str), rn["percent"], color="#e69f00")
    for v, pct in zip(rn["value"].astype(str), rn["percent"]):
        ax.text(v, pct + 1, f"{pct:.1f}%" if pct >= 0.05 else f"{pct:.3f}%", ha="center", fontsize=7.5)
    ax.set_ylim(0, 108)
    ax.set_xlabel("return number")
    ax.set_ylabel("share of points (%)")
    ax.set_title("(d) Return numbers", fontsize=10)
    C.crs_note(fig, "Withheld and noise points are left out of the density maps.")
    fig.savefig(C.FIG / "fig01_point_density_qa.png", dpi=170)
    C.shrink_png(C.FIG / "fig01_point_density_qa.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
