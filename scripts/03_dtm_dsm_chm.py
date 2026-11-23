"""Step 3. Grid the point cloud into a DTM, a DSM, a canopy height model and woody cover.

DTM  = mean elevation of ground-classified points (class 2) in each 1 m cell.
       Cells with no ground point are filled in two stages: small gaps take the mean of
       their filled neighbours, larger gaps (water, buildings, dense canopy) are filled by
       linear interpolation on a triangulation (TIN) of the cells around the gap.
DSM  = highest first return in each 1 m cell (noise and withheld points left out).
CHM  = DSM minus DTM, with a simple noise filter.
Woody cover = CHM above 2 m, after removing features narrower than 2 m (wires, fences,
       poles) and patches that look like buildings (few multiple returns, smooth top).

Outputs: 1 m rasters in data/work/, outputs/tables/vegetation_cover_summary.csv,
outputs/tables/dtm_gridding_summary.csv, fig02 and fig03.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

import common as C

K8 = np.ones((3, 3), bool)
HEIGHT_CUTS = (0.5, 1.0, 2.0, 5.0, 10.0)   # m above ground
WOODY_MIN_H = 2.0                          # m, woody plant height threshold
CHM_MAX_H = 40.0                           # m, anything higher is treated as noise


def fill_small_gaps(a, passes=2):
    """Fill NaN cells that touch valid cells with the mean of their valid 3x3 neighbours."""
    a = a.copy()
    filled = np.zeros(a.shape, bool)
    for _ in range(passes):
        nan = np.isnan(a)
        if not nan.any():
            break
        v = np.where(nan, 0.0, a)
        s = ndi.uniform_filter(v, 3, mode="constant") * 9
        n = ndi.uniform_filter((~nan).astype("float64"), 3, mode="constant") * 9
        ok = nan & (n >= 3)          # need at least 3 valid neighbours
        a[ok] = s[ok] / n[ok]
        filled |= ok
    return a, filled


def fill_tin(a):
    """Fill the remaining NaN cells by linear interpolation on a TIN of the surrounding cells."""
    nan = np.isnan(a)
    if not nan.any():
        return a, nan
    ring = ~nan & ndi.binary_dilation(nan, K8, iterations=2)
    rr, cc = np.nonzero(ring)
    vr, vc = np.nonzero(nan)
    interp = LinearNDInterpolator(np.c_[cc, rr], a[rr, cc])
    vals = interp(np.c_[vc, vr])
    bad = np.isnan(vals)             # outside the convex hull of the ring cells
    if bad.any():
        vals[bad] = NearestNDInterpolator(np.c_[cc, rr], a[rr, cc])(np.c_[vc[bad], vr[bad]])
    out = a.copy()
    out[vr, vc] = vals
    return out, nan


def main():
    g = C.load_grid()
    shape = (g["ny"], g["nx"])
    sum_g = np.zeros(shape)
    n_g = np.zeros(shape, np.int32)
    n_all = np.zeros(shape, np.int32)
    n_multi = np.zeros(shape, np.int32)
    n_water = np.zeros(shape, np.int32)
    dsm = np.full(shape, np.nan, "float32")
    for name in C.LAZ_NAMES:
        p = C.read_tile(C.RAW / "laz" / name)
        idx, ok = C.cell_index(p["x"], p["y"], g)
        ok &= ~p["withheld"] & ~np.isin(p["cls"], C.NOISE_CLASSES)
        m = ok & (p["cls"] == 2)
        sum_g += C.grid_sum(idx[m], p["z"][m], g)
        n_g += C.grid_count(idx[m], g).astype(np.int32)
        n_all += C.grid_count(idx[ok], g).astype(np.int32)
        n_multi += C.grid_count(idx[ok & (p["nr"] > 1)], g).astype(np.int32)
        n_water += C.grid_count(idx[ok & (p["cls"] == 9)], g).astype(np.int32)
        m = ok & (p["rn"] == 1)
        dsm = np.fmax(dsm, C.grid_extreme(idx[m], p["z"][m].astype("float32"), g, "max"))
        print(f"  gridded {name}")

    # ---- DTM -------------------------------------------------------------------------
    dtm_obs = np.where(n_g > 0, sum_g / np.maximum(n_g, 1), np.nan)
    dtm1, small = fill_small_gaps(dtm_obs)
    dtm, tin = fill_tin(dtm1)
    fill_type = np.zeros(shape, np.uint8)        # 0 observed, 1 small-gap mean, 2 TIN
    fill_type[small] = 1
    fill_type[tin] = 2

    # ---- DSM -------------------------------------------------------------------------
    dsm_obs_void = np.isnan(dsm)
    dsm1, _ = fill_small_gaps(dsm.astype("float64"))
    dsm_f = np.where(np.isnan(dsm1), dtm, dsm1)  # no first return at all: use the ground
    dsm_f = np.maximum(dsm_f, dtm)

    # ---- CHM with a simple noise filter ----------------------------------------------
    chm_raw = dsm_f - dtm
    med = ndi.median_filter(chm_raw, 3)
    high_nb = ndi.uniform_filter((chm_raw > WOODY_MIN_H).astype("float32"), 3) * 9
    too_high = chm_raw > CHM_MAX_H
    spike = (chm_raw - med > 5.0) & (high_nb <= 2)     # lone high cell: bird, wire, pole top
    chm = np.where(too_high | spike, med, chm_raw)
    chm = np.clip(chm, 0, CHM_MAX_H).astype("float32")

    # ---- woody cover ------------------------------------------------------------------
    tall = chm > WOODY_MIN_H
    tall_open = ndi.binary_opening(tall, np.ones((2, 2), bool))   # drop features under 2 m wide
    lab, n = ndi.label(tall_open, K8)
    ids = np.arange(1, n + 1)
    area = ndi.sum(tall_open, lab, ids)
    multi_frac = ndi.mean(n_multi > 0, lab, ids)
    rough = ndi.mean(np.abs(chm - ndi.median_filter(chm, 3)), lab, ids)
    is_bldg = (area >= 20) & (multi_frac < 0.7) & (rough < 0.2)
    bldg = np.r_[False, is_bldg][lab]
    woody = tall_open & ~bldg

    # water-classified cells (class 9), closed and cleaned
    water9 = ndi.binary_opening(ndi.binary_closing(n_water > 0, K8, iterations=2), K8)

    # ---- save -------------------------------------------------------------------------
    t = g["transform"]
    C.write_raster(C.WORK / "dtm_1m.tif", dtm.astype("float32"), t)
    C.write_raster(C.WORK / "dsm_1m.tif", dsm_f.astype("float32"), t)
    C.write_raster(C.WORK / "chm_1m.tif", chm, t)
    C.write_raster(C.WORK / "dtm_fill_type_1m.tif", fill_type, t)
    C.write_raster(C.WORK / "woody_mask_1m.tif", woody, t)
    C.write_raster(C.WORK / "building_like_mask_1m.tif", bldg, t)
    C.write_raster(C.WORK / "water_class9_mask_1m.tif", water9, t)
    C.write_raster(C.WORK / "n_multi_1m.tif", np.minimum(n_multi, 255).astype("uint8"), t)

    # ---- tables -----------------------------------------------------------------------
    ncell = dtm.size
    grid_rows = [
        ("Cells in the working grid", ncell),
        ("Cells with at least one ground point (%)", 100 * (fill_type == 0).mean()),
        ("Cells filled from neighbours, small gaps (%)", 100 * (fill_type == 1).mean()),
        ("Cells filled by TIN interpolation, large gaps (%)", 100 * (fill_type == 2).mean()),
        ("Cells with no first return before filling (%)", 100 * dsm_obs_void.mean()),
        ("Cells with water-classified points, class 9 (%)", 100 * (n_water > 0).mean()),
        ("CHM cells replaced by the spike filter", int(spike.sum())),
        ("CHM cells above the 40 m limit", int(too_high.sum())),
        ("DTM minimum (m NAVD88)", float(dtm.min())), ("DTM maximum (m NAVD88)", float(dtm.max())),
    ]
    pd.DataFrame(grid_rows, columns=["metric", "value"]).to_csv(
        C.TAB / "dtm_gridding_summary.csv", index=False, float_format="%.4f")

    veg_rows = []
    for h in HEIGHT_CUTS:
        veg_rows.append((f"Surface cover above {h:g} m, all objects (CHM, unfiltered)",
                         (chm > h).sum() / 1e4, 100 * (chm > h).mean()))
    veg_rows += [
        ("Features above 2 m and under 2 m wide, removed (wires, fences, poles)",
         (tall & ~tall_open).sum() / 1e4, 100 * (tall & ~tall_open).mean()),
        ("Building-like patches above 2 m, removed", bldg.sum() / 1e4, 100 * bldg.mean()),
        ("Woody cover above 2 m (filtered)", woody.sum() / 1e4, 100 * woody.mean()),
        ("Woody cover above 5 m (filtered)", (woody & (chm > 5)).sum() / 1e4, 100 * (woody & (chm > 5)).mean()),
        ("Woody cover above 10 m (filtered)", (woody & (chm > 10)).sum() / 1e4, 100 * (woody & (chm > 10)).mean()),
    ]
    veg = pd.DataFrame(veg_rows, columns=["cover_class", "area_ha", "percent_of_area"])
    hw = chm[woody]
    extra = pd.DataFrame([("Woody canopy height, mean (m)", np.nan, np.nan, hw.mean()),
                          ("Woody canopy height, median (m)", np.nan, np.nan, np.median(hw)),
                          ("Woody canopy height, 95th percentile (m)", np.nan, np.nan, np.percentile(hw, 95)),
                          ("Woody canopy height, maximum (m)", np.nan, np.nan, hw.max()),
                          ("Building-like patches (count)", np.nan, np.nan, int(is_bldg.sum())),
                          ("Woody patches (count)", np.nan, np.nan, int((~is_bldg).sum()))],
                         columns=["cover_class", "area_ha", "percent_of_area", "value"])
    veg = pd.concat([veg, extra], ignore_index=True)
    veg.to_csv(C.TAB / "vegetation_cover_summary.csv", index=False, float_format="%.3f")
    print(pd.DataFrame(grid_rows).to_string(index=False, header=False))
    print(veg.to_string(index=False))

    # ---- figures ----------------------------------------------------------------------
    f = 3
    ext = C.map_extent(g)
    hs_t = C.block_reduce(C.hillshade(dtm), f)
    hs_s = C.block_reduce(C.hillshade(dsm_f), f)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 7.0))
    ax = axes[0]
    ax.imshow(hs_t, extent=ext, cmap="gray", vmin=0.2, vmax=1)
    im = ax.imshow(C.block_reduce(dtm, f), extent=ext, cmap="terrain", alpha=0.55)
    ax.set_title("(a) DTM: bare earth from ground points", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, label="elevation (m NAVD88)")
    ax = axes[1]
    ax.imshow(hs_s, extent=ext, cmap="gray", vmin=0.2, vmax=1)
    im = ax.imshow(C.block_reduce(dsm_f, f, np.nanmax), extent=ext, cmap="terrain", alpha=0.55)
    ax.set_title("(b) DSM: highest first return", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, label="elevation (m NAVD88)")
    ax = axes[2]
    c3 = C.block_reduce(chm, f, np.nanmax)
    im = ax.imshow(np.where(c3 > 0.5, c3, np.nan), extent=ext, cmap="YlGn", vmin=0, vmax=20)
    ax.set_facecolor("0.92")
    ax.set_title("(c) CHM: DSM minus DTM (above 0.5 m)", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, extend="max",
                 label="height above ground (m)")
    for k, ax in enumerate(axes):
        C.finish_map(ax, g)
        if k:
            ax.set_ylabel("")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.07, wspace=0.18)
    C.crs_note(fig, "Maps drawn at 3 m (maximum per cell for DSM and CHM).")
    fig.savefig(C.FIG / "fig02_dtm_dsm_chm.png", dpi=160)
    C.shrink_png(C.FIG / "fig02_dtm_dsm_chm.png")
    plt.close(fig)

    wf = C.block_reduce(woody.astype("float32"), 10, np.mean) * 100
    fig = plt.figure(figsize=(12, 7.0))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.05], left=0.06, right=0.98, top=0.94, bottom=0.09,
                          wspace=0.3, hspace=0.38)
    ax = fig.add_subplot(gs[:, 0])
    ax.imshow(hs_t, extent=ext, cmap="gray", vmin=0.2, vmax=1, alpha=0.6)
    im = ax.imshow(np.where(wf > 0, wf, np.nan), extent=ext, cmap="YlGn", vmin=0, vmax=100)
    ax.set_title("(a) Woody cover above 2 m, 10 m cells", fontsize=10)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, label="woody cover (% of cell)")
    C.finish_map(ax, g)
    ax = fig.add_subplot(gs[:, 1])
    ax.imshow(hs_t, extent=ext, cmap="gray", vmin=0.2, vmax=1, alpha=0.6)
    cls = np.zeros(shape, np.uint8)
    cls[tall & ~tall_open] = 1
    cls[bldg] = 2
    cls3 = C.block_reduce(cls, f, np.max).astype(float)
    cls3[cls3 == 0] = np.nan
    from matplotlib.colors import ListedColormap
    ax.imshow(cls3, extent=ext, cmap=ListedColormap(["#0072b2", "#d55e00"]), vmin=1, vmax=2,
              interpolation="nearest")
    ax.set_title("(b) Tall objects left out of woody cover", fontsize=10)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#0072b2", label="narrow features (under 2 m wide)"),
                       Patch(color="#d55e00", label="building-like patches")],
              loc="upper left", fontsize=7.5, framealpha=0.9)
    C.finish_map(ax, g)
    ax.set_ylabel("")
    ax = fig.add_subplot(gs[0, 2])
    ax.hist(hw, bins=np.arange(2, 31, 1), color="#2e7d32", weights=np.full(hw.size, 1e-4))
    ax.set_xlabel("canopy height (m)")
    ax.set_ylabel("area (ha)")
    ax.set_title("(c) Height of woody canopy cells", fontsize=10)
    ax = fig.add_subplot(gs[1, 2])
    labs = [f">{h:g} m" for h in HEIGHT_CUTS]
    raw_pct = [100 * (chm > h).mean() for h in HEIGHT_CUTS]
    wood_pct = [100 * (woody & (chm > h)).mean() if h >= WOODY_MIN_H else np.nan for h in HEIGHT_CUTS]
    x = np.arange(len(labs))
    ax.bar(x - 0.2, raw_pct, 0.4, color="0.6", label="all objects (CHM)")
    ax.bar(x + 0.2, wood_pct, 0.4, color="#2e7d32", label="woody (filtered)")
    ax.set_xticks(x, labs)
    ax.set_ylabel("share of area (%)")
    ax.set_xlabel("height threshold")
    ax.set_title("(d) Cover above height thresholds", fontsize=10)
    ax.legend(fontsize=8)
    C.crs_note(fig, "LiDAR flown leaf-off on 25 December 2022.")
    fig.savefig(C.FIG / "fig03_woody_cover.png", dpi=160)
    C.shrink_png(C.FIG / "fig03_woody_cover.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
