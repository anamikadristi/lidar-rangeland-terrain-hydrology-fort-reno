"""Step 7. Compare the derived streams and depressions with NHD, and relate woody cover to drainage.

1. Streams: distance from points on NHD flowlines to the nearest derived stream cell, and
   from derived stream cells to the nearest NHD flowline.
2. Ponds: how many NHD waterbodies are picked up by the depression inventory, and how many
   depressions have no NHD waterbody.
3. Woody cover by distance from the derived streams, by TWI class and by catchment.

Outputs: outputs/tables/stream_comparison_nhd.csv, nhd_waterbody_detection.csv,
pond_detection_summary.csv, woody_cover_by_stream_distance.csv, woody_cover_by_twi.csv,
woody_cover_by_catchment.csv and fig08.
"""
import warnings

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import ndimage as ndi

import common as C

warnings.filterwarnings("ignore")
STEP = 5.0   # m between sample points on NHD lines
FCODE = {46003: "Stream, intermittent", 46006: "Stream, perennial", 55800: "Artificial path (through waterbody)",
         46000: "Stream", 33600: "Canal or ditch", 33400: "Connector"}


def sample_points(gdf, step=STEP):
    xs, ys, codes = [], [], []
    for geom, code in zip(gdf.geometry, gdf["fcode"]):
        for part in getattr(geom, "geoms", [geom]):
            if part.length == 0:
                continue
            for d in np.arange(step / 2, part.length, step):
                p = part.interpolate(d)
                xs.append(p.x), ys.append(p.y), codes.append(code)
    return np.array(xs), np.array(ys), np.array(codes)


def dist_stats(d, label, extra=None):
    row = dict(subset=label, n_samples=int(d.size), length_km=d.size * STEP / 1e3,
               median_m=np.median(d), mean_m=d.mean(), p90_m=np.percentile(d, 90),
               within_10m_pct=100 * (d <= 10).mean(), within_20m_pct=100 * (d <= 20).mean(),
               within_50m_pct=100 * (d <= 50).mean())
    if extra:
        row.update(extra)
    return row


def main():
    g = C.load_grid()
    t = g["transform"]
    streams = C.read_raster(C.WORK / "streams_1m.tif").astype(bool)
    acc = C.read_raster(C.WORK / "flowacc_cells_1m.tif")
    woody = C.read_raster(C.WORK / "woody_mask_1m.tif").astype(bool)
    dep_id = C.read_raster(C.WORK / "depression_id_1m.tif")
    depth = C.read_raster(C.WORK / "depression_depth_1m.tif")
    catch = C.read_raster(C.WORK / "catchments_1m.tif")
    pond_mask = C.read_raster(C.WORK / "pond_mask_1m.tif").astype(bool)
    fl = C.nhd_layer("NHDFlowline")
    wb = C.nhd_layer("NHDWaterbody")

    # ---- 1. streams -------------------------------------------------------------------------
    d_to_derived = ndi.distance_transform_edt(~streams) * C.CELL
    x, y, code = sample_points(fl)
    col = np.clip(((x - t.c) / C.CELL).astype(int), 0, g["nx"] - 1)
    row = np.clip(((t.f - y) / C.CELL).astype(int), 0, g["ny"] - 1)
    d1 = d_to_derived[row, col]
    a_here = ndi.maximum_filter(acc, 21)[row, col] * 1e-4      # largest contributing area within 10 m, ha
    rows = [dist_stats(d1, "NHD flowlines to derived streams: all NHD")]
    for c in sorted(set(code.tolist())):
        rows.append(dist_stats(d1[code == c], f"NHD flowlines to derived streams: {FCODE.get(int(c), int(c))}"))
    big = a_here >= 5.0
    rows.append(dist_stats(d1[big], "NHD flowlines to derived streams: NHD points where my DTM has 5 ha or more within 10 m"))
    rows.append(dist_stats(d1[~big], "NHD flowlines to derived streams: NHD points where my DTM has under 5 ha within 10 m"))
    nhd_r = C.rasterize(fl.geometry, g, all_touched=True).astype(bool)
    d_to_nhd = ndi.distance_transform_edt(~nhd_r) * C.CELL
    d2 = d_to_nhd[streams]
    r2 = dist_stats(d2, "Derived stream cells to NHD flowlines: all derived")
    r2["length_km"] = np.nan
    rows.append(r2)
    sc = pd.DataFrame(rows)
    sc.to_csv(C.TAB / "stream_comparison_nhd.csv", index=False, float_format="%.2f")
    print(sc.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    # ---- 2. ponds -----------------------------------------------------------------------------
    wb = wb.copy()
    wb["area_in_study_m2"] = wb.area
    wb["touches_edge"] = ~wb.within(gpd.GeoSeries.from_wkt(
        [f"POLYGON(({g['xmin'] + 2} {g['ymin'] + 2},{g['xmax'] - 2} {g['ymin'] + 2},{g['xmax'] - 2} {g['ymax'] - 2},"
         f"{g['xmin'] + 2} {g['ymax'] - 2},{g['xmin'] + 2} {g['ymin'] + 2}))"], crs=C.DST_EPSG).iloc[0])
    wid = C.rasterize(wb.geometry, g, all_touched=True, values=np.arange(1, len(wb) + 1))
    ids = np.arange(1, len(wb) + 1)
    n_cells = ndi.sum(wid > 0, wid, ids)
    wb["depression_cover_frac"] = ndi.sum(dep_id > 0, wid, ids) / np.maximum(n_cells, 1)
    wb["max_depression_depth_m"] = ndi.maximum(np.where(dep_id > 0, depth, 0), wid, ids)
    wb["detected"] = wb["depression_cover_frac"] >= 0.25
    wb["ftype_name"] = wb["ftype"].map({390: "Lake or pond", 436: "Reservoir", 466: "Swamp or marsh"}).fillna("Other")
    cols = ["permanent_identifier", "gnis_name", "ftype_name", "fcode", "area_in_study_m2", "touches_edge",
            "depression_cover_frac", "max_depression_depth_m", "detected"]
    wb[cols].sort_values("area_in_study_m2", ascending=False).to_csv(
        C.TAB / "nhd_waterbody_detection.csv", index=False, float_format="%.3f")
    dep = pd.read_csv(C.TAB / "depressions_inventory.csv")
    prow = []
    for lab, m in [("All NHD waterbodies in the study area", np.ones(len(wb), bool)),
                   ("NHD lake or pond", (wb["ftype"] == 390).values),
                   ("NHD reservoir (treatment and storage basins)", (wb["ftype"] == 436).values),
                   ("NHD waterbodies of 1,000 m2 or more", (wb["area_in_study_m2"] >= 1000).values),
                   ("NHD waterbodies under 1,000 m2", (wb["area_in_study_m2"] < 1000).values)]:
        if m.sum() == 0:
            continue
        prow.append(dict(group=lab, count=int(m.sum()), area_ha=wb.loc[m, "area_in_study_m2"].sum() / 1e4,
                         matched_by_depression=int((m & wb["detected"].values).sum()),
                         matched_pct=100 * (m & wb["detected"].values).sum() / max(m.sum(), 1),
                         any_overlap=int((m & (wb["depression_cover_frac"] > 0).values).sum())))
    for lo in (100, 500, 2000):
        m = dep["area_m2"] >= lo
        prow.append(dict(group=f"Depressions of {lo:,} m2 or more", count=int(m.sum()),
                         area_ha=dep.loc[m, "area_m2"].sum() / 1e4,
                         matched_by_depression=np.nan, matched_pct=np.nan,
                         any_overlap=int((m & dep["overlaps_nhd_waterbody"]).sum())))
    pond = pd.DataFrame(prow)
    pond.to_csv(C.TAB / "pond_detection_summary.csv", index=False, float_format="%.2f")
    print(pond.to_string(index=False))

    # ---- 3. woody cover and drainage ----------------------------------------------------------
    bins = [0, 10, 25, 50, 100, 200, 1e9]
    wrows = []
    tot_w = woody.sum()
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (d_to_derived >= lo) & (d_to_derived < hi)
        wrows.append(dict(distance_to_derived_stream_m=f"{lo:g} to {hi:g}" if hi < 1e8 else f"above {lo:g}",
                          area_ha=m.sum() / 1e4, share_of_study_area_pct=100 * m.mean(),
                          woody_cover_pct=100 * woody[m].mean(), woody_area_ha=woody[m].sum() / 1e4,
                          share_of_all_woody_cover_pct=100 * woody[m].sum() / tot_w))
    wd = pd.DataFrame(wrows)
    wd.to_csv(C.TAB / "woody_cover_by_stream_distance.csv", index=False, float_format="%.2f")
    print(wd.to_string(index=False))

    twi = C.read_raster(C.WORK / "twi_3m.tif")
    w3 = C.block_reduce(woody.astype("float32"), 3, np.mean)[:twi.shape[0], :twi.shape[1]]
    twi = twi[:w3.shape[0], :w3.shape[1]]
    trows = []
    for lo, hi in [(-99, 6), (6, 8), (8, 10), (10, 12), (12, 99)]:
        m = (twi >= lo) & (twi < hi)
        trows.append(dict(twi_class=f"{lo:g} to {hi:g}" if lo > -90 and hi < 90 else (f"under {hi:g}" if lo < -90 else f"{lo:g} or more"),
                          share_of_study_area_pct=100 * m.mean(), woody_cover_pct=100 * w3[m].mean()))
    tw = pd.DataFrame(trows)
    tw.to_csv(C.TAB / "woody_cover_by_twi.csv", index=False, float_format="%.2f")
    print(tw.to_string(index=False))

    outl = pd.read_csv(C.TAB / "outlets_catchments.csv")
    crow = []
    for _, o in outl.iterrows():
        m = catch == o["outlet_id"]
        crow.append(dict(outlet_id=int(o["outlet_id"]), edge=o["edge"], catchment_area_km2=o["catchment_area_km2"],
                         woody_cover_pct=100 * woody[m].mean(),
                         woody_cover_within_25m_of_stream_pct=100 * woody[m & (d_to_derived < 25)].mean(),
                         pond_storage_m3=float(depth[m & pond_mask].sum()),
                         pond_storage_mm=float(depth[m & pond_mask].sum() / m.sum() * 1000)))
    cw = pd.DataFrame(crow)
    cw.to_csv(C.TAB / "woody_cover_by_catchment.csv", index=False, float_format="%.2f")
    print(cw.to_string(index=False))

    # ---- figure -------------------------------------------------------------------------------
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.05, 1, 1], left=0.06, right=0.98, top=0.94, bottom=0.10,
                          wspace=0.33, hspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    hs = C.block_reduce(C.hillshade(C.read_raster(C.WORK / "dtm_1m.tif").astype("float64")), 3)
    ax.imshow(hs, extent=C.map_extent(g), cmap="gray", vmin=0.2, vmax=1, alpha=0.6)
    w3m = C.block_reduce(woody.astype("float32"), 3, np.max)
    ax.imshow(np.where(w3m > 0, 1.0, np.nan), extent=C.map_extent(g), cmap="Greens", vmin=0, vmax=1.3,
              interpolation="nearest")
    sg = gpd.read_file(C.WORK / "hydro_vectors.gpkg", layer="derived_streams")
    sg.plot(ax=ax, color="#0033aa", lw=0.8)
    det = wb[wb["detected"]]
    und = wb[~wb["detected"]]
    if len(det):
        det.plot(ax=ax, facecolor="#56b4e9", edgecolor="k", lw=0.4)
    if len(und):
        und.plot(ax=ax, facecolor="#e69f00", edgecolor="k", lw=0.4)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#2e7d32", label="woody cover above 2 m"),
                       Line2D([], [], color="#0033aa", lw=1, label="derived streams (5 ha)"),
                       Patch(facecolor="#56b4e9", edgecolor="k", label="NHD waterbody, matched"),
                       Patch(facecolor="#e69f00", edgecolor="k", label="NHD waterbody, not matched")],
              loc="upper left", fontsize=7.5, framealpha=0.9)
    C.finish_map(ax, g)
    ax.set_title("(a) Woody cover, derived streams and NHD waterbodies", fontsize=10)

    ax = fig.add_subplot(gs[0, 1])
    colors = {46003: "#e69f00", 46006: "#0072b2", 55800: "#009e73"}
    xs = np.linspace(0, 150, 301)
    for c in sorted(set(code.tolist())):
        dd = np.sort(d1[code == c])
        ax.plot(xs, 100 * np.searchsorted(dd, xs, side="right") / dd.size, color=colors.get(int(c), "k"),
                label=f"{FCODE.get(int(c), int(c))} ({dd.size * STEP / 1e3:.1f} km)")
    dd = np.sort(d1)
    ax.plot(xs, 100 * np.searchsorted(dd, xs, side="right") / dd.size, color="k", lw=2, label="all NHD flowlines")
    ax.set_xlabel("distance from NHD flowline to nearest derived stream (m)")
    ax.set_ylabel("share of NHD length (%)")
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower right")
    ax.set_title("(b) NHD flowlines against derived streams", fontsize=10)

    ax = fig.add_subplot(gs[0, 2])
    ax.bar(np.arange(len(wd)), wd["woody_cover_pct"], color="#2e7d32")
    for i, v in enumerate(wd["share_of_all_woody_cover_pct"]):
        ax.text(i, wd["woody_cover_pct"][i] + 0.5, f"{v:.0f}%", ha="center", fontsize=8)
    ax.set_xticks(np.arange(len(wd)), [s.replace(" to ", "-").replace("above ", ">") for s in wd["distance_to_derived_stream_m"]],
                  fontsize=8)
    ax.set_xlabel("distance to nearest derived stream (m)")
    ax.set_ylabel("woody cover (% of zone)")
    ax.set_ylim(0, wd["woody_cover_pct"].max() * 1.18)
    ax.set_title("(c) Woody cover by distance to streams\n(labels: share of all woody cover)", fontsize=10)

    ax = fig.add_subplot(gs[1, 1])
    m = dep["treatment"] == "kept as pond or basin"
    ax.scatter(dep.loc[~m, "area_m2"], dep.loc[~m, "volume_m3"], s=12, c="#d55e00", alpha=0.6,
               label="dry depression (breached)")
    ax.scatter(dep.loc[m, "area_m2"], dep.loc[m, "volume_m3"], s=22, c="#0072b2", alpha=0.9, label="pond or basin (kept)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("depression area (m$^2$)")
    ax.set_ylabel("storage to spill level (m$^3$)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7.5, loc="upper left")
    ax.set_title("(d) Depression area and storage", fontsize=10)

    ax = fig.add_subplot(gs[1, 2])
    ax.bar(np.arange(len(tw)), tw["woody_cover_pct"], color="#33658a")
    for i, v in enumerate(tw["share_of_study_area_pct"]):
        ax.text(i, tw["woody_cover_pct"][i] + 0.3, f"{v:.0f}%", ha="center", fontsize=8)
    ax.set_xticks(np.arange(len(tw)), tw["twi_class"], fontsize=8)
    ax.set_xlabel("TWI class (3 m cells)")
    ax.set_ylabel("woody cover (% of class)")
    ax.set_ylim(0, tw["woody_cover_pct"].max() * 1.18)
    ax.set_title("(e) Woody cover by wetness index\n(labels: share of study area)", fontsize=10)
    C.crs_note(fig)
    fig.savefig(C.FIG / "fig08_nhd_comparison_woody_drainage.png", dpi=160)
    C.shrink_png(C.FIG / "fig08_nhd_comparison_woody_drainage.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
