"""Step 6. Hydrologic analysis of my 1 m DTM.

1. Depression filling with my own priority-flood routine (common.priority_flood). The
   depth of fill gives the closed depressions (ponds and other hollows).
2. D8 flow direction and flow accumulation with pysheds on the filled surface.
3. Stream network: cells with a contributing area of STREAM_THRESHOLD_HA or more.
4. Outlets where the largest flow paths leave the area, and their catchments.
5. Topographic wetness index (TWI) on a 3 m DTM with D-infinity flow accumulation.
6. A wider-area check: the same outlets on the USGS DEM tile, which extends beyond my
   area, to see how much of each catchment lies outside the LiDAR tiles I processed.

Outputs: rasters in data/work/, vectors in data/work/hydro_vectors.gpkg, tables in
outputs/tables/, fig06 and fig07.
"""
import warnings

import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from affine import Affine
from rasterio import features
from scipy import ndimage as ndi
from shapely.geometry import LineString, Point, shape

import common as C
from pysheds.grid import Grid
from pysheds.sview import Raster, ViewFinder

warnings.filterwarnings("ignore")

STREAM_THRESHOLD_HA = 5.0       # contributing area that starts a stream
SENSITIVITY_HA = (1.0, 2.0, 5.0, 10.0, 20.0, 30.0)
OUTLET_MIN_HA = 25.0            # smallest edge outlet that gets its own catchment
EPS = 1e-5                      # m per cell, gradient added across filled areas
DEP_MIN_AREA = 100.0            # m2, smallest depression kept in the inventory
DEP_MIN_DEPTH = 0.15            # m, smallest maximum depth kept (above LiDAR noise)
WET_FRAC = 0.30                 # share of a depression that must be water to call it a pond or basin
BREACH_MAX_LEN = 300.0          # m, longest cut allowed through a barrier
K8 = np.ones((3, 3), bool)


def route(z, transform, routing="d8"):
    """Fill z with a small gradient and return (grid, flow direction, accumulation in cells)."""
    filled = C.priority_flood(z, EPS)
    vf = ViewFinder(affine=transform, shape=z.shape, nodata=np.nan,
                    crs=rasterio.crs.CRS.from_epsg(C.DST_EPSG))
    grid = Grid(viewfinder=vf)
    r = Raster(filled, viewfinder=vf)
    fdir = grid.flowdir(r, routing=routing)
    acc = grid.accumulation(fdir, routing=routing)
    return grid, fdir, np.asarray(acc, dtype="float64")


def mask_raster(mask, transform):
    vf = ViewFinder(affine=transform, shape=mask.shape, nodata=False, crs=rasterio.crs.CRS.from_epsg(C.DST_EPSG))
    return Raster(mask.astype(bool), viewfinder=vf)


def edge_outlets(acc, min_cells, sep=50):
    """Largest accumulation cells on the grid edge, at least `sep` cells apart."""
    ny, nx = acc.shape
    rr = np.r_[np.zeros(nx, int), np.full(nx, ny - 1), np.arange(ny), np.arange(ny)]
    cc = np.r_[np.arange(nx), np.arange(nx), np.zeros(ny, int), np.full(ny, nx - 1)]
    side = np.r_[np.full(nx, "north"), np.full(nx, "south"), np.full(ny, "west"), np.full(ny, "east")]
    a = acc[rr, cc]
    out = []
    for k in np.argsort(-a):
        if a[k] < min_cells:
            break
        if all(max(abs(rr[k] - o[0]), abs(cc[k] - o[1])) > sep for o in out):
            out.append((int(rr[k]), int(cc[k]), side[k], float(a[k])))
    return out


def main():
    g = C.load_grid()
    t = g["transform"]
    cell_ha = C.CELL ** 2 / 1e4
    dtm = C.read_raster(C.WORK / "dtm_1m.tif").astype("float64")
    area_km2 = dtm.size * C.CELL ** 2 / 1e6

    # ---- 1. depressions ----------------------------------------------------------------
    # Inventory on the DTM as gridded (fill only). Each depression is then classed as a
    # pond or basin (kept as storage) or as a dry hollow behind a barrier (breached).
    fill_type = C.read_raster(C.WORK / "dtm_fill_type_1m.tif")
    chm = C.read_raster(C.WORK / "chm_1m.tif")
    bldg = C.read_raster(C.WORK / "building_like_mask_1m.tif").astype(bool)
    water9 = C.read_raster(C.WORK / "water_class9_mask_1m.tif").astype(bool)
    # water seen by the LiDAR: class 9 points, or large gaps in the ground returns with nothing above them
    lidar_water = water9 | ((fill_type == 2) & (chm < 0.5) & ~bldg)
    wb = C.nhd_layer("NHDWaterbody")
    nhd_water = C.rasterize(wb.geometry, g, all_touched=True).astype(bool)

    def inventory(z):
        filled = C.priority_flood(z, 0.0)
        depth = filled - z
        lab, n = ndi.label(depth > 0, K8)
        ids = np.arange(1, n + 1)
        a_m2 = ndi.sum(depth > 0, lab, ids) * C.CELL ** 2
        dmax = ndi.maximum(depth, lab, ids)
        vol = ndi.sum(depth, lab, ids) * C.CELL ** 2
        keep = (a_m2 >= DEP_MIN_AREA) & (dmax >= DEP_MIN_DEPTH)
        relab = np.zeros(n + 1, np.int32)
        order = ids[keep][np.argsort(-vol[keep])]
        relab[order] = np.arange(1, order.size + 1)
        dep_lab = relab[lab]
        k = np.arange(1, order.size + 1)
        com = ndi.center_of_mass(dep_lab > 0, dep_lab, k)
        pit = ndi.minimum_position(z, dep_lab, k)
        df = pd.DataFrame(dict(
            dep_id=k,
            x_utm=[t.c + (c[1] + 0.5) * C.CELL for c in com],
            y_utm=[t.f - (c[0] + 0.5) * C.CELL for c in com],
            area_m2=ndi.sum(dep_lab > 0, dep_lab, k) * C.CELL ** 2,
            max_depth_m=ndi.maximum(depth, dep_lab, k),
            mean_depth_m=ndi.mean(depth, dep_lab, k),
            volume_m3=ndi.sum(depth, dep_lab, k) * C.CELL ** 2,
            spill_elevation_m=ndi.maximum(filled, dep_lab, k),
            lidar_water_frac=ndi.mean(lidar_water, dep_lab, k),
            nhd_waterbody_frac=ndi.mean(nhd_water, dep_lab, k),
            pit_row=[p[0] for p in pit], pit_col=[p[1] for p in pit]))
        df["is_pond_or_basin"] = (df["lidar_water_frac"] >= WET_FRAC) | (df["nhd_waterbody_frac"] >= WET_FRAC)
        return filled, depth, dep_lab, df, int(n)

    filled, depth, dep_lab, dep, n = inventory(dtm)
    dep["overlaps_nhd_waterbody"] = dep["nhd_waterbody_frac"] > 0

    # Breach the dry depressions. A few passes are needed because a breach can expose a new one.
    zb = dtm.copy()
    breaches = []
    cur = dep
    for it in range(1, 7):
        todo = cur[~cur["is_pond_or_basin"]]
        n_ok = 0
        for _, d in todo.sort_values("spill_elevation_m").iterrows():
            ok, length, cut, er, ec = C.breach_path(zb, d["pit_row"], d["pit_col"], BREACH_MAX_LEN, radius=int(BREACH_MAX_LEN))
            if ok:
                n_ok += 1
            breaches.append(dict(pass_no=it, dep_id_in_pass=int(d["dep_id"]), breached=bool(ok),
                                 x_utm=t.c + (d["pit_col"] + 0.5) * C.CELL, y_utm=t.f - (d["pit_row"] + 0.5) * C.CELL,
                                 depression_area_m2=d["area_m2"], depression_volume_m3=d["volume_m3"],
                                 breach_length_m=length, max_cut_m=cut))
        print(f"breach pass {it}: {len(todo)} dry depressions, {n_ok} breached")
        if n_ok == 0:
            break
        _, _, _, cur, _ = inventory(zb)
    br = pd.DataFrame(breaches)
    br.to_csv(C.TAB / "breaches.csv", index=False, float_format="%.3f")
    first = br[br["pass_no"] == 1].set_index("dep_id_in_pass")
    dep["treatment"] = np.where(dep["is_pond_or_basin"], "kept as pond or basin",
                                np.where(dep["dep_id"].map(first["breached"]).fillna(False).astype(bool),
                                         "breached (dry, behind a barrier)", "not breached (no outlet within 300 m)"))
    dep["breach_length_m"] = dep["dep_id"].map(first["breach_length_m"])
    dep["breach_max_cut_m"] = dep["dep_id"].map(first["max_cut_m"])
    dep.drop(columns=["pit_row", "pit_col"]).to_csv(C.TAB / "depressions_inventory.csv", index=False, float_format="%.3f")
    polys = [(shape(s), int(v)) for s, v in features.shapes(dep_lab, mask=dep_lab > 0, transform=t, connectivity=8)]
    dep_gdf = gpd.GeoDataFrame(dict(dep_id=[v for _, v in polys]), geometry=[p for p, _ in polys],
                               crs=C.DST_EPSG).dissolve("dep_id").reset_index().merge(
                                   dep.drop(columns=["pit_row", "pit_col"]), on="dep_id")

    # what is left after breaching
    filled_b, depth_b, lab_b, dep_b, _ = inventory(zb)
    size_rows = []
    for name, m in [("All depressions on the unbreached DTM", np.ones(len(dep), bool)),
                    ("Kept as pond or basin", (dep["treatment"] == "kept as pond or basin").values),
                    ("Breached (dry, behind a barrier)", dep["treatment"].str.startswith("breached").values),
                    ("Not breached (no outlet within 300 m)", dep["treatment"].str.startswith("not").values)]:
        size_rows.append(dict(group=name, count=int(m.sum()), area_ha=dep.loc[m, "area_m2"].sum() / 1e4,
                              volume_m3=dep.loc[m, "volume_m3"].sum(),
                              max_depth_m=dep.loc[m, "max_depth_m"].max() if m.any() else np.nan,
                              overlapping_nhd_waterbody=int((m & dep["overlaps_nhd_waterbody"].values).sum())))
    for lo, hi in [(100, 500), (500, 2000), (2000, 10000), (10000, 1e12)]:
        m = ((dep["area_m2"] >= lo) & (dep["area_m2"] < hi) & (dep["treatment"] == "kept as pond or basin")).values
        size_rows.append(dict(group=f"Pond or basin, {lo:g} to {hi:g} m2" if hi < 1e11 else f"Pond or basin, above {lo:g} m2",
                              count=int(m.sum()), area_ha=dep.loc[m, "area_m2"].sum() / 1e4,
                              volume_m3=dep.loc[m, "volume_m3"].sum(),
                              max_depth_m=dep.loc[m, "max_depth_m"].max() if m.any() else np.nan,
                              overlapping_nhd_waterbody=int((m & dep["overlaps_nhd_waterbody"].values).sum())))
    size_rows.append(dict(group="Depressions left on the breached DTM (used for routing)", count=len(dep_b),
                          area_ha=dep_b["area_m2"].sum() / 1e4, volume_m3=dep_b["volume_m3"].sum(),
                          max_depth_m=dep_b["max_depth_m"].max(), overlapping_nhd_waterbody=int((dep_b["nhd_waterbody_frac"] > 0).sum())))
    size_rows.append(dict(group="All filled cells on the unbreached DTM, including small pits", count=n,
                          area_ha=(depth > 0).sum() * cell_ha, volume_m3=depth.sum() * C.CELL ** 2,
                          max_depth_m=depth.max(), overlapping_nhd_waterbody=np.nan))
    pd.DataFrame(size_rows).to_csv(C.TAB / "depressions_summary.csv", index=False, float_format="%.2f")
    print(pd.DataFrame(size_rows).to_string(index=False))
    pond = dep["treatment"] == "kept as pond or basin"
    pond_lab = np.r_[False, pond.values][dep_lab]

    # ---- 2. flow direction and accumulation (D8) ---------------------------------------------
    grid, fdir, acc = route(zb, t, "d8")
    F = np.asarray(fdir)
    interior_bad = (F[1:-1, 1:-1] < 0).sum()
    print(f"cells with no flow direction away from the edge: {interior_bad}")

    # ---- 3. streams ---------------------------------------------------------------------------
    thr_cells = STREAM_THRESHOLD_HA / cell_ha
    streams = acc >= thr_cells
    net = grid.extract_river_network(fdir, mask_raster(streams, t))
    lines = [LineString(f["geometry"]["coordinates"]) for f in net["features"] if len(f["geometry"]["coordinates"]) > 1]
    stream_gdf = gpd.GeoDataFrame(dict(seg_id=np.arange(1, len(lines) + 1)), geometry=lines, crs=C.DST_EPSG)
    stream_gdf["length_m"] = stream_gdf.length
    # contributing area at the downstream end of each segment
    ends = [ln.coords[-1] for ln in lines]
    starts = [ln.coords[0] for ln in lines]

    def acc_at(xy):
        c = int((xy[0] - t.c) / C.CELL)
        r = int((t.f - xy[1]) / C.CELL)
        return acc[min(max(r, 0), g["ny"] - 1), min(max(c, 0), g["nx"] - 1)]

    stream_gdf["contrib_area_ha"] = [max(acc_at(a), acc_at(b)) * cell_ha for a, b in zip(starts, ends)]
    nhd_fl = C.nhd_layer("NHDFlowline")
    sens = []
    for ha in SENSITIVITY_HA:
        s = acc >= ha / cell_ha
        nt = grid.extract_river_network(fdir, mask_raster(s, t))
        L = sum(LineString(f["geometry"]["coordinates"]).length for f in nt["features"]
                if len(f["geometry"]["coordinates"]) > 1)
        sens.append(dict(threshold_ha=ha, stream_cells=int(s.sum()), stream_length_km=L / 1e3,
                         drainage_density_km_per_km2=L / 1e3 / area_km2))
    sens.append(dict(threshold_ha=np.nan, stream_cells=np.nan, stream_length_km=nhd_fl.length.sum() / 1e3,
                     drainage_density_km_per_km2=nhd_fl.length.sum() / 1e3 / area_km2))
    sens = pd.DataFrame(sens)
    sens.insert(0, "network", [f"Derived, {h:g} ha threshold" for h in SENSITIVITY_HA] + ["NHD flowlines (all types)"])
    sens.to_csv(C.TAB / "stream_threshold_sensitivity.csv", index=False, float_format="%.3f")
    print(sens.to_string(index=False))

    # ---- 4. outlets and catchments ------------------------------------------------------------
    outs = edge_outlets(acc, OUTLET_MIN_HA / cell_ha)
    catch = C.label_catchments(F.astype(np.int64), np.array([o[0] for o in outs]),
                               np.array([o[1] for o in outs])).astype(np.int16)
    orows = []
    for k, (r, c, side, a) in enumerate(outs, 1):
        m = catch == k
        on_edge = np.r_[m[0], m[-1], m[:, 0], m[:, -1]].sum()
        ms = m & streams
        orows.append(dict(outlet_id=k, edge=side, x_utm=t.c + (c + 0.5) * C.CELL, y_utm=t.f - (r + 0.5) * C.CELL,
                          catchment_area_km2=m.sum() * C.CELL ** 2 / 1e6,
                          share_of_study_area_pct=100 * m.mean(),
                          catchment_cells_on_area_edge_m=float(on_edge) * C.CELL,
                          stream_length_km=float(stream_gdf.length[stream_gdf.representative_point().apply(
                              lambda p: m[min(int((t.f - p.y) / C.CELL), g["ny"] - 1),
                                          min(int((p.x - t.c) / C.CELL), g["nx"] - 1)])].sum() / 1e3),
                          mean_elevation_m=float(dtm[m].mean()),
                          pond_storage_m3=float(depth[m & pond_lab].sum() * C.CELL ** 2)))
    outlets = pd.DataFrame(orows)

    # ---- 6. wider-area check on the USGS DEM tile (3 m) ---------------------------------------
    with rasterio.open(C.RAW / "dem" / C.DEM_NAME) as src:
        big = src.read(1).astype("float64")
        big[big == src.nodata] = np.nan
        bt = src.transform
    f = 3
    off_c = int(round((g["xmin"] - bt.c))) % f      # keep the 3 m cells aligned with my grid
    off_r = int(round((bt.f - g["ymax"]))) % f
    big = big[off_r:, off_c:]
    big3 = C.block_reduce(big, f, np.nanmean)
    big3 = np.where(np.isfinite(big3), big3, np.nanmin(big3) - 50.0)   # nodata acts as an open low edge
    t3 = Affine(f, 0, bt.c + off_c, 0, -f, bt.f - off_r)
    _, _, acc_big = route(big3, t3, "d8")
    wide = []
    for _, o in outlets.iterrows():
        c3 = int((o["x_utm"] - t3.c) / f)
        r3 = int((t3.f - o["y_utm"]) / f)
        win = acc_big[max(r3 - 3, 0):r3 + 4, max(c3 - 3, 0):c3 + 4]    # snap within 9 m
        wide.append(win.max() * f * f / 1e6)
    wide = np.array(wide)
    # The wider DEM is only filled, not breached, so a flow path can leave it at a different
    # place near embankments. I keep the wider value only where it is at least 90 % of the
    # area found inside the study area; otherwise it is left blank (not resolved).
    wide[wide < 0.9 * outlets["catchment_area_km2"].values] = np.nan
    outlets["area_on_wider_usgs_dem_km2"] = wide
    outlets["share_inside_study_area_pct"] = (100 * outlets["catchment_area_km2"] /
                                              outlets["area_on_wider_usgs_dem_km2"]).clip(upper=100)
    outlets.to_csv(C.TAB / "outlets_catchments.csv", index=False, float_format="%.3f")
    print(outlets.to_string(index=False))
    out_gdf = gpd.GeoDataFrame(outlets, geometry=[Point(x, y) for x, y in zip(outlets.x_utm, outlets.y_utm)],
                               crs=C.DST_EPSG)
    cpolys = [(shape(s), int(v)) for s, v in features.shapes(catch, mask=catch > 0, transform=t, connectivity=8)]
    catch_gdf = gpd.GeoDataFrame(dict(outlet_id=[v for _, v in cpolys]), geometry=[p for p, _ in cpolys],
                                 crs=C.DST_EPSG).dissolve("outlet_id").reset_index()
    catch_gdf["geometry"] = catch_gdf.buffer(0)
    catch_gdf = catch_gdf.merge(outlets[["outlet_id", "edge", "catchment_area_km2",
                                         "area_on_wider_usgs_dem_km2", "share_inside_study_area_pct"]], on="outlet_id")

    # ---- 5. topographic wetness index on a 3 m DTM, D-infinity ---------------------------------
    dtm3 = C.block_reduce(zb, f, np.mean)
    ta = Affine(f, 0, t.c, 0, -f, t.f)
    _, _, acc_inf = route(dtm3, ta, "dinf")
    dzdx, dzdy = C.horn_gradients(dtm3, f)
    tanb = np.maximum(np.hypot(dzdx, dzdy), 0.001)
    sca = np.maximum(np.where(np.isfinite(acc_inf), acc_inf, 1.0), 1.0) * f          # specific catchment area, m2 per m
    twi = np.log(sca / tanb).astype("float32")
    twi_rows = [("TWI cell size (m)", f), ("TWI minimum", twi.min()), ("TWI 10th percentile", np.percentile(twi, 10)),
                ("TWI median", np.median(twi)), ("TWI 90th percentile", np.percentile(twi, 90)),
                ("TWI 99th percentile", np.percentile(twi, 99)), ("TWI maximum", twi.max())]
    for v in (8, 10, 12, 14):
        twi_rows.append((f"Area with TWI of {v} or more (%)", 100 * (twi >= v).mean()))
        twi_rows.append((f"Area with TWI of {v} or more (ha)", (twi >= v).sum() * f * f / 1e4))
    pd.DataFrame(twi_rows, columns=["metric", "value"]).to_csv(C.TAB / "twi_summary.csv", index=False,
                                                              float_format="%.3f")
    print(pd.DataFrame(twi_rows).to_string(index=False, header=False))

    # ---- save ---------------------------------------------------------------------------------
    C.write_raster(C.WORK / "dtm_breached_1m.tif", zb.astype("float32"), t)
    C.write_raster(C.WORK / "pond_mask_1m.tif", pond_lab, t)
    C.write_raster(C.WORK / "depression_depth_1m.tif", depth.astype("float32"), t)
    C.write_raster(C.WORK / "depression_id_1m.tif", dep_lab.astype("int32"), t, nodata=0)
    C.write_raster(C.WORK / "flowdir_d8_1m.tif", F.astype("int16"), t, nodata=0)
    C.write_raster(C.WORK / "flowacc_cells_1m.tif", acc.astype("float32"), t)
    C.write_raster(C.WORK / "streams_1m.tif", streams, t)
    C.write_raster(C.WORK / "catchments_1m.tif", catch, t, nodata=0)
    C.write_raster(C.WORK / "twi_3m.tif", twi, ta)
    gp = C.WORK / "hydro_vectors.gpkg"
    if gp.exists():
        gp.unlink()
    stream_gdf.to_file(gp, layer="derived_streams", driver="GPKG")
    catch_gdf.to_file(gp, layer="catchments", driver="GPKG")
    out_gdf.to_file(gp, layer="outlets", driver="GPKG")
    dep_gdf.to_file(gp, layer="depressions", driver="GPKG")

    pv = dep.loc[pond, "volume_m3"].sum()
    summ = [("Study area (km2)", area_km2),
            ("Stream threshold (ha)", STREAM_THRESHOLD_HA),
            ("Derived stream length (km)", stream_gdf.length.sum() / 1e3),
            ("Derived drainage density (km per km2)", stream_gdf.length.sum() / 1e3 / area_km2),
            ("Edge outlets with 25 ha or more", len(outlets)),
            ("Area draining to the largest outlet (km2)", outlets["catchment_area_km2"].max()),
            ("Area draining to the largest outlet (% of study area)", outlets["share_of_study_area_pct"].max()),
            ("Area in delineated catchments (% of study area)", 100 * (catch > 0).mean()),
            ("Depressions on the unbreached DTM (count)", len(dep)),
            ("Apparent storage of all depressions on the unbreached DTM (m3)", dep["volume_m3"].sum()),
            ("Ponds and basins kept (count)", int(pond.sum())),
            ("Pond and basin area at spill level (ha)", dep.loc[pond, "area_m2"].sum() / 1e4),
            ("Pond and basin storage between survey surface and spill level (m3)", pv),
            ("Pond and basin storage as depth over the study area (mm)", pv / (area_km2 * 1e6) * 1e3),
            ("Dry depressions breached in pass 1 (count)", int(dep["treatment"].str.startswith("breached").sum())),
            ("Apparent storage removed by breaching, pass 1 (m3)",
             dep.loc[dep["treatment"].str.startswith("breached"), "volume_m3"].sum()),
            ("Breaches in all passes (count)", int(br["breached"].sum())),
            ("Median breach length (m)", br.loc[br["breached"], "breach_length_m"].median()),
            ("Longest breach (m)", br.loc[br["breached"], "breach_length_m"].max()),
            ("Deepest breach cut (m)", br.loc[br["breached"], "max_cut_m"].max()),
            ("Dry depressions that could not be breached (count)", int(dep["treatment"].str.startswith("not").sum()))]
    pd.DataFrame(summ, columns=["metric", "value"]).to_csv(C.TAB / "hydrology_summary.csv", index=False,
                                                           float_format="%.3f")
    print(pd.DataFrame(summ).to_string(index=False, header=False))

    # ---- figures ------------------------------------------------------------------------------
    ext = C.map_extent(g)
    hs = C.block_reduce(C.hillshade(dtm), 3)
    fig, axes = plt.subplots(1, 3, figsize=(13, 7.2))
    ax = axes[0]
    a3 = C.block_reduce(acc, 6, np.max) * cell_ha
    im = ax.imshow(np.log10(np.maximum(a3, 1e-4)), extent=ext, cmap="Blues", vmin=-3, vmax=3)
    cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, ticks=[-3, -2, -1, 0, 1, 2, 3])
    cb.ax.set_xticklabels(["0.001", "0.01", "0.1", "1", "10", "100", "1000"])
    cb.set_label("D8 contributing area (ha, log scale)")
    ax.set_title("(a) Flow accumulation", fontsize=10)
    ax = axes[1]
    ax.imshow(hs, extent=ext, cmap="gray", vmin=0.2, vmax=1, alpha=0.7)
    wb.plot(ax=ax, color="#56b4e9", edgecolor="none", alpha=0.9)
    nhd_fl.plot(ax=ax, color="#d55e00", lw=1.6, label="NHD flowlines")
    stream_gdf.plot(ax=ax, color="#0033aa", lw=0.8, label=f"derived streams ({STREAM_THRESHOLD_HA:g} ha)")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    ax.legend(handles=[Line2D([], [], color="#0033aa", lw=1, label=f"derived streams ({STREAM_THRESHOLD_HA:g} ha)"),
                       Line2D([], [], color="#d55e00", lw=1.6, label="NHD flowlines"),
                       Patch(color="#56b4e9", label="NHD waterbodies")], loc="upper left", fontsize=7.5, framealpha=0.9)
    ax.set_title("(b) Derived streams and NHD", fontsize=10)
    ax = axes[2]
    ax.imshow(hs, extent=ext, cmap="gray", vmin=0.2, vmax=1, alpha=0.7)
    cm = plt.get_cmap("tab20", max(len(outlets), 1))
    c3 = catch[1::3, 1::3].astype(float)
    c3[c3 == 0] = np.nan
    ax.imshow(c3, extent=ext, cmap=cm, vmin=0.5, vmax=len(outlets) + 0.5, alpha=0.55, interpolation="nearest")
    stream_gdf.plot(ax=ax, color="#0033aa", lw=0.7)
    ax.scatter(outlets.x_utm, outlets.y_utm, s=28, c="k", marker="v", zorder=10, clip_on=False)
    for _, o in outlets.iterrows():
        dx = -70 if o["edge"] == "east" else (70 if o["edge"] == "west" else 0)
        dy = -90 if o["edge"] == "north" else (90 if o["edge"] == "south" else 0)
        ax.text(o["x_utm"] + dx, o["y_utm"] + dy, f"{int(o['outlet_id'])}", fontsize=8, ha="center", va="center",
                fontweight="bold", bbox=dict(facecolor="w", alpha=0.8, edgecolor="none", pad=1))
    ax.set_title(f"(c) Catchments of edge outlets ({OUTLET_MIN_HA:g} ha or more)", fontsize=10)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0.5, len(outlets) + 0.5)), ax=ax,
                 orientation="horizontal", pad=0.08, shrink=0.9, label="outlet number (see outlets_catchments.csv)",
                 ticks=range(1, len(outlets) + 1))
    for k, ax in enumerate(axes):
        C.finish_map(ax, g)
        if k:
            ax.set_ylabel("")
    fig.subplots_adjust(left=0.055, right=0.99, top=0.95, bottom=0.07, wspace=0.18)
    C.crs_note(fig, "Uncoloured areas in (c) drain to smaller edge outlets.")
    fig.savefig(C.FIG / "fig06_flow_streams_catchments.png", dpi=150)
    C.shrink_png(C.FIG / "fig06_flow_streams_catchments.png")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 7.4))
    ax = axes[0]
    im = ax.imshow(C.block_reduce(twi, 2, np.mean), extent=ext, cmap="YlGnBu", vmin=4, vmax=14)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, extend="both",
                 label="topographic wetness index, ln(a / tan b)")
    ax.set_title("(a) TWI, D-infinity (computed at 3 m, drawn at 6 m)", fontsize=10)
    ax = axes[1]
    ax.imshow(hs, extent=ext, cmap="gray", vmin=0.2, vmax=1, alpha=0.7)
    d3 = C.block_reduce(np.where(pond_lab, depth, 0.0), 3, np.max)
    b3 = C.block_reduce(((dep_lab > 0) & ~pond_lab).astype("float32"), 3, np.max)
    ax.imshow(np.where(b3 > 0, 1.0, np.nan), extent=ext, cmap="Oranges", vmin=0, vmax=2, interpolation="nearest")
    im = ax.imshow(np.where(d3 > 0, d3, np.nan), extent=ext, cmap="viridis_r", vmin=0, vmax=2)
    wb.boundary.plot(ax=ax, color="k", lw=0.5)
    ax.legend(handles=[Line2D([], [], color="k", lw=0.8, label="NHD waterbody outline"),
                       Patch(color="#fdae6b", label="dry depression, breached for routing")],
              loc="upper left", fontsize=7.5, framealpha=0.9)
    fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.08, shrink=0.9, extend="max",
                 label="pond or basin: depth below spill level (m)")
    ax.set_title("(b) Closed depressions (100 m$^2$ and 0.15 m or more)", fontsize=10)
    for k, ax in enumerate(axes):
        C.finish_map(ax, g)
        if k:
            ax.set_ylabel("")
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.07, wspace=0.2)
    C.crs_note(fig)
    fig.savefig(C.FIG / "fig07_twi_depressions.png", dpi=150)
    C.shrink_png(C.FIG / "fig07_twi_depressions.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
