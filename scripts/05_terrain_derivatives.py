"""Step 5. Terrain derivatives from my 1 m DTM: slope, aspect, hillshade and curvature.

Slope and aspect use Horn's 3 by 3 finite differences on the 1 m DTM. Curvature is a
second derivative, so it amplifies small elevation noise. I compute it on a DTM smoothed
with a Gaussian filter (sigma 2 m) and use the Zevenbergen and Thorne (1987) formulas.
Sign convention: positive profile curvature is convex (flow speeds up), negative is
concave (flow slows down). Positive plan curvature is divergent (ridges, noses), negative
is convergent (hollows, swales).

Outputs: 1 m rasters in data/work/, outputs/tables/terrain_summary.csv,
outputs/figures/fig05_terrain_derivatives.png.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

import common as C


def curvatures(z, cell):
    """Profile and plan curvature (1/m) from second derivatives of a smoothed surface."""
    zy, zx = np.gradient(z, cell)
    zy = -zy
    zxy_, zxx = np.gradient(zx, cell)
    zyy_, zyx = np.gradient(zy, cell)
    zyy = -zyy_
    zxy = -zxy_
    p2, q2 = zx ** 2, zy ** 2
    s = p2 + q2
    s = np.where(s < 1e-10, np.nan, s)
    prof = -(p2 * zxx + 2 * zx * zy * zxy + q2 * zyy) / (s * (1 + s) ** 1.5)
    plan = -(q2 * zxx - 2 * zx * zy * zxy + p2 * zyy) / (s ** 1.5)
    return np.nan_to_num(prof).astype("float32"), np.nan_to_num(plan).astype("float32")


def main():
    g = C.load_grid()
    dtm = C.read_raster(C.WORK / "dtm_1m.tif").astype("float64")
    dzdx, dzdy = C.horn_gradients(dtm, C.CELL)
    grad = np.hypot(dzdx, dzdy)
    slope_deg = np.degrees(np.arctan(grad)).astype("float32")
    slope_pct = (100 * grad).astype("float32")
    aspect = (np.degrees(np.arctan2(-dzdx, -dzdy)) % 360).astype("float32")  # downslope, clockwise from north
    aspect[grad < 1e-4] = -1
    hs = C.hillshade(dtm).astype("float32")
    zs = ndi.gaussian_filter(dtm, 2.0, mode="nearest")
    prof, plan = curvatures(zs, C.CELL)

    C.write_raster(C.WORK / "slope_deg_1m.tif", slope_deg, g["transform"])
    C.write_raster(C.WORK / "aspect_deg_1m.tif", aspect, g["transform"], nodata=-9999.0)
    C.write_raster(C.WORK / "hillshade_1m.tif", (hs * 255).astype("uint8"), g["transform"])
    C.write_raster(C.WORK / "profile_curvature_1m.tif", prof, g["transform"])
    C.write_raster(C.WORK / "plan_curvature_1m.tif", plan, g["transform"])

    rows = [("Elevation minimum (m NAVD88)", dtm.min()), ("Elevation maximum (m NAVD88)", dtm.max()),
            ("Elevation mean (m NAVD88)", dtm.mean()), ("Relief (m)", dtm.max() - dtm.min()),
            ("Slope mean (degrees)", slope_deg.mean()), ("Slope median (degrees)", np.median(slope_deg)),
            ("Slope 95th percentile (degrees)", np.percentile(slope_deg, 95)),
            ("Slope mean (percent)", slope_pct.mean()), ("Slope median (percent)", np.median(slope_pct))]
    for lo, hi in [(0, 2), (2, 5), (5, 10), (10, 20), (20, 1e9)]:
        lab = f"Area with slope {lo} to {hi} percent (%)" if hi < 1e8 else f"Area with slope above {lo} percent (%)"
        rows.append((lab, 100 * ((slope_pct >= lo) & (slope_pct < hi)).mean()))
    names = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    sect = (np.floor(((aspect + 22.5) % 360) / 45)).astype(int)
    valid = aspect >= 0
    for k, nm in enumerate(names):
        rows.append((f"Area facing {nm} (%)", 100 * (sect[valid] == k).mean()))
    rows.append(("Area with convergent plan curvature below -0.05 per m (%)", 100 * (plan < -0.05).mean()))
    rows.append(("Area with divergent plan curvature above 0.05 per m (%)", 100 * (plan > 0.05).mean()))
    pd.DataFrame(rows, columns=["metric", "value"]).to_csv(C.TAB / "terrain_summary.csv", index=False,
                                                          float_format="%.3f")
    print(pd.DataFrame(rows, columns=["metric", "value"]).to_string(index=False))

    f = 5  # plot at 5 m to keep the figure small
    ext = C.map_extent(g)
    fig, axes = plt.subplots(1, 4, figsize=(13.5, 6.3))
    ax = axes[0]
    ax.imshow(C.block_reduce(hs, f), extent=ext, cmap="gray", vmin=0.15, vmax=1)
    im = ax.imshow(C.block_reduce(dtm, f), extent=ext, cmap="terrain", alpha=0.55)
    ax.set_title("(a) Elevation over hillshade\n(light from NW, 45 degrees)", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, label="elevation (m NAVD88)")
    ax = axes[1]
    im = ax.imshow(C.block_reduce(slope_deg, f), extent=ext, cmap="magma_r", vmin=0, vmax=15)
    ax.set_title("(b) Slope", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, label="slope (degrees)", extend="max")
    ax = axes[2]
    asp3 = aspect[5::10, 5::10].astype("float32")   # one cell in ten, to keep the picture readable
    asp3[asp3 < 0] = np.nan
    im = ax.imshow(asp3, extent=ext, cmap="twilight", vmin=0, vmax=360, interpolation="nearest")
    ax.set_title("(c) Aspect (downslope direction)\n(sampled every 10 m)", fontsize=10)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9,
                      label="aspect (degrees clockwise from north)", ticks=[0, 90, 180, 270, 360])
    cb.ax.set_xticklabels(["N 0", "E 90", "S 180", "W 270", "N 360"])
    ax = axes[3]
    im = ax.imshow(C.block_reduce(np.clip(plan, -0.1, 0.1), 10), extent=ext, cmap="BrBG_r", vmin=-0.03, vmax=0.03)
    ax.set_title("(d) Plan curvature\n(sigma 2 m, drawn as 10 m mean)", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, extend="both",
                 label="plan curvature (1/m)\nnegative convergent, positive divergent")
    for k, ax in enumerate(axes):
        C.finish_map(ax, g)
        if k:
            ax.set_ylabel("")
    fig.subplots_adjust(left=0.04, right=0.99, top=0.93, bottom=0.07, wspace=0.16)
    C.crs_note(fig)
    fig.savefig(C.FIG / "fig05_terrain_derivatives.png", dpi=150)
    C.shrink_png(C.FIG / "fig05_terrain_derivatives.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
