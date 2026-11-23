"""Step 4. Compare my DTM with the USGS 3DEP 1 m DEM.

The USGS DEM for this area was made by the LiDAR contractor from the same ground points
plus hydro breaklines, on a TIN, and then resampled by USGS to a 1 m UTM grid. So this
comparison shows how much my gridding differs from the official product. It is not a
check against surveyed ground points, and it cannot show the absolute vertical accuracy.

I report the difference (my DTM minus USGS DEM) for all cells, for land cells outside NHD
waterbodies, by how the cell was filled, by slope class and by vegetation class.

Outputs: outputs/tables/accuracy_*.csv, outputs/figures/fig04_dtm_vs_usgs_dem.png,
data/work/usgs_dem_1m.tif and data/work/dtm_minus_usgs_1m.tif.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.windows import from_bounds

import common as C


def main():
    g = C.load_grid()
    dtm = C.read_raster(C.WORK / "dtm_1m.tif").astype("float64")
    with rasterio.open(C.RAW / "dem" / C.DEM_NAME) as src:
        assert src.crs.to_epsg() == C.DST_EPSG
        w = from_bounds(g["xmin"], g["ymin"], g["xmax"], g["ymax"], src.transform)
        off = (w.col_off - round(w.col_off), w.row_off - round(w.row_off))
        print(f"USGS DEM pixel grid offset from my grid: {off[0]:.5f}, {off[1]:.5f} pixels")
        usgs = src.read(1, window=w.round_offsets().round_lengths()).astype("float64")
        usgs[usgs == src.nodata] = np.nan
    assert usgs.shape == dtm.shape
    C.write_raster(C.WORK / "usgs_dem_1m.tif", usgs.astype("float32"), g["transform"])
    diff = dtm - usgs
    C.write_raster(C.WORK / "dtm_minus_usgs_1m.tif", diff.astype("float32"), g["transform"])

    fill_type = C.read_raster(C.WORK / "dtm_fill_type_1m.tif")
    chm = C.read_raster(C.WORK / "chm_1m.tif")
    woody = C.read_raster(C.WORK / "woody_mask_1m.tif").astype(bool)
    bldg = C.read_raster(C.WORK / "building_like_mask_1m.tif").astype(bool)
    wb = C.nhd_layer("NHDWaterbody")
    water = C.rasterize(wb.geometry, g, all_touched=True).astype(bool)
    dzdx, dzdy = C.horn_gradients(usgs, C.CELL)
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))
    land = ~water

    rows = []

    def add(group, name, mask):
        s = C.error_stats(diff[mask])
        s.update(group=group, subset=name, share_of_area_pct=100 * mask.mean())
        rows.append(s)

    add("overall", "All cells", np.ones(dtm.shape, bool))
    add("overall", "Land cells (outside NHD waterbodies)", land)
    add("overall", "Cells inside NHD waterbodies", water)
    add("overall", "Land cells that have a ground point", land & (fill_type == 0))
    add("fill type", "Ground point in cell", fill_type == 0)
    add("fill type", "Small gap filled from neighbours", fill_type == 1)
    add("fill type", "Large gap filled by TIN", fill_type == 2)
    add("fill type", "Large gap filled by TIN, land only", (fill_type == 2) & land)
    for lo, hi in [(0, 1), (1, 2), (2, 5), (5, 10), (10, 20), (20, 91)]:
        nm = f"Slope {lo} to {hi} degrees" if hi < 90 else f"Slope above {lo} degrees"
        add("slope (land)", nm, land & (slope >= lo) & (slope < hi))
    veg_classes = [
        ("Open, under 0.5 m", land & (chm <= 0.5)),
        ("Low vegetation, 0.5 to 2 m", land & (chm > 0.5) & (chm <= 2.0)),
        ("Woody cover, 2 to 5 m", land & woody & (chm <= 5)),
        ("Woody cover, 5 to 10 m", land & woody & (chm > 5) & (chm <= 10)),
        ("Woody cover, above 10 m", land & woody & (chm > 10)),
        ("Building-like patches", land & bldg),
    ]
    for nm, m in veg_classes:
        add("vegetation (land)", nm, m)
    df = pd.DataFrame(rows)
    cols = ["group", "subset", "n", "share_of_area_pct", "bias_m", "median_m", "std_m", "rmse_m", "mae_m",
            "nmad_m", "p95_abs_m", "within_5cm_pct", "within_10cm_pct", "within_19p6cm_pct"]
    df[cols].to_csv(C.TAB / "accuracy_dtm_vs_usgs_dem.csv", index=False, float_format="%.4f")
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ---- figure ---------------------------------------------------------------------------
    fig = plt.figure(figsize=(12.5, 7.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1, 1], left=0.06, right=0.98, top=0.94, bottom=0.15,
                          wspace=0.32, hspace=0.55)
    ax = fig.add_subplot(gs[:, 0])
    im = ax.imshow(C.block_reduce(diff, 3), extent=C.map_extent(g), cmap="RdBu_r", vmin=-0.15, vmax=0.15)
    wb.boundary.plot(ax=ax, color="k", lw=0.5)
    C.finish_map(ax, g)
    ax.set_title("(a) My DTM minus USGS 1 m DEM\n(3 m mean; black lines: NHD waterbodies)", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, extend="both",
                 label="elevation difference (m)")
    ax = fig.add_subplot(gs[0, 1])
    d = diff[land]
    ax.hist(d, bins=np.linspace(-0.3, 0.3, 121), color="#33658a")
    ax.set_yscale("log")
    ax.axvline(0, color="k", lw=0.6)
    s = df[df["subset"] == "Land cells (outside NHD waterbodies)"].iloc[0]
    ax.set_title("(b) Differences on land cells", fontsize=10)
    ax.text(0.03, 0.96, f"n = {int(s['n']):,}\nbias = {s['bias_m'] * 100:+.1f} cm\nRMSE = {s['rmse_m'] * 100:.1f} cm\n"
            f"NMAD = {s['nmad_m'] * 100:.1f} cm", transform=ax.transAxes, va="top", fontsize=8.5)
    ax.set_xlabel("difference (m)")
    ax.set_ylabel("cells (log scale)")

    def bars(ax, group, title):
        sub = df[df["group"] == group]
        x = np.arange(len(sub))
        ax.bar(x - 0.2, sub["rmse_m"] * 100, 0.4, color="#d55e00", label="RMSE")
        ax.bar(x + 0.2, sub["bias_m"] * 100, 0.4, color="#0072b2", label="bias")
        ax.axhline(0, color="k", lw=0.6)
        labs = [t.replace("Slope ", "").replace(" degrees", "°").replace("Woody cover, ", "woody ")
                .replace("Low vegetation, ", "low veg ").replace("Building-like patches", "building-like")
                .replace("Open, ", "open ") for t in sub["subset"]]
        ax.set_xticks(x, labs, rotation=30, ha="right", fontsize=7.5)
        ax.set_ylabel("difference (cm)")
        ax.set_title(title, fontsize=10)
        top = float(sub["rmse_m"].max() * 100)
        for xi, n, v in zip(x, sub["share_of_area_pct"], sub["rmse_m"] * 100):
            ax.text(xi - 0.2, v + 0.02 * top, f"{n:.1f}%", ha="center", va="bottom", fontsize=7, color="0.3")
        ax.set_ylim(top=top * 1.3)
        ax.legend(fontsize=8, loc="upper left", ncol=2)

    bars(fig.add_subplot(gs[0, 2]), "slope (land)", "(c) By slope class (labels: share of area)")
    bars(fig.add_subplot(gs[1, 1]), "vegetation (land)", "(d) By vegetation class")
    ax = fig.add_subplot(gs[1, 2])
    sub = df[df["group"] == "fill type"].iloc[[0, 1, 3]]
    x = np.arange(len(sub))
    ax.bar(x - 0.2, sub["rmse_m"] * 100, 0.4, color="#d55e00", label="RMSE")
    ax.bar(x + 0.2, sub["bias_m"] * 100, 0.4, color="#0072b2", label="bias")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x, ["ground point\nin cell", "small gap\n(neighbour mean)", "large gap\n(TIN, land)"], fontsize=8)
    ax.set_ylabel("difference (cm)")
    ax.set_title("(e) By how the DTM cell was made", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    C.crs_note(fig, "Both surfaces come from the same LiDAR survey.")
    fig.savefig(C.FIG / "fig04_dtm_vs_usgs_dem.png", dpi=160)
    C.shrink_png(C.FIG / "fig04_dtm_vs_usgs_dem.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
