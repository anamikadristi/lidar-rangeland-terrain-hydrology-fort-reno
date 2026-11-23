"""Step 8. Export the key products at a size that fits in the repository.

The full 1 m float rasters are too large for git, so I export:
  - the DTM at 2 m as a Cloud-Optimized GeoTIFF with elevations stored as whole centimetres
    (unsigned 16-bit, scale 0.01 and offset written in the file),
  - the other rasters as Cloud-Optimized GeoTIFFs at 3 m or 10 m, also as scaled integers,
  - one GeoPackage with the study area, derived streams, catchments, outlets,
    depressions and the NHD features clipped to the study area.

Everything is in NAD83 / UTM zone 14N (EPSG:26914), heights in m NAVD88.
"""
import geopandas as gpd
import numpy as np
import pandas as pd
from affine import Affine
from shapely.geometry import box

import common as C


def cog(name, arr, transform, nodata, desc, units, res, method, rows, scale=None, offset=None):
    path = C.RAS / name
    C.write_raster(path, arr, transform, nodata=nodata, cog=True, scale=scale, offset=offset)
    rows.append(dict(file=f"outputs/rasters/{name}", content=desc, units=units, cell_size_m=res,
                     resampling_from_1m=method, dtype=str(arr.dtype), size_mb=round(path.stat().st_size / 1e6, 2),
                     crs=f"EPSG:{C.DST_EPSG}"))


def main():
    g = C.load_grid()
    t1 = g["transform"]
    t3 = Affine(3, 0, t1.c, 0, -3, t1.f)
    t10 = Affine(10, 0, t1.c, 0, -10, t1.f)
    rows = []
    for p in C.RAS.glob("*.tif"):
        p.unlink()

    t2 = Affine(2, 0, t1.c, 0, -2, t1.f)
    dtm = C.read_raster(C.WORK / "dtm_1m.tif")
    off = float(np.floor(dtm.min() / 10) * 10)
    cog("dtm_2m_cm.tif", np.round((C.block_reduce(dtm, 2, np.mean) - off) * 100).astype("uint16"), t2, 65535,
        "Digital terrain model (bare earth). Elevation in m = value * 0.01 + offset", "m NAVD88 (stored as cm)", 2,
        "mean", rows, scale=0.01, offset=off)
    chm = C.read_raster(C.WORK / "chm_1m.tif")
    cog("chm_3m_max_cm.tif", np.round(C.block_reduce(chm, 3, np.max) * 100).astype("uint16"), t3, 65535,
        "Canopy height model (DSM minus DTM, filtered). Height in m = value * 0.01", "m above ground (stored as cm)",
        3, "maximum", rows, scale=0.01, offset=0.0)
    woody = C.read_raster(C.WORK / "woody_mask_1m.tif").astype("float32")
    cog("woody_cover_gt2m_10m_pct.tif", np.round(C.block_reduce(woody, 10, np.mean) * 100).astype("uint8"), t10, 255,
        "Woody cover above 2 m (filtered)", "percent of cell", 10, "mean", rows)
    slope = C.read_raster(C.WORK / "slope_deg_1m.tif")
    cog("slope_3m_decideg.tif", np.round(C.block_reduce(slope, 3, np.mean) * 10).astype("uint16"), t3, 65535,
        "Slope (Horn method on the 1 m DTM). Degrees = value * 0.1", "degrees (stored as tenths)", 3, "mean", rows,
        scale=0.1, offset=0.0)
    diff = C.read_raster(C.WORK / "dtm_minus_usgs_1m.tif")
    d3 = np.clip(np.round(C.block_reduce(diff, 3, np.nanmean) * 100), -32000, 32000).astype("int16")
    cog("dtm_minus_usgs_dem_3m_cm.tif", d3, t3, -32768,
        "My DTM minus the USGS 3DEP 1 m DEM. Difference in m = value * 0.01", "m (stored as cm)", 3, "mean", rows,
        scale=0.01, offset=0.0)
    acc = C.read_raster(C.WORK / "flowacc_cells_1m.tif")
    cog("flow_accumulation_3m_m2.tif", C.block_reduce(acc, 3, np.max).astype("uint32"), t3, 0,
        "D8 contributing area on the breached and filled 1 m DTM", "m2", 3, "maximum", rows)
    twi = C.read_raster(C.WORK / "twi_3m.tif")
    cog("twi_3m_x10.tif", np.round(twi * 10).astype("int16"), t3, -32768,
        "Topographic wetness index, D-infinity, computed at 3 m. TWI = value * 0.1", "dimensionless (stored x10)", 3,
        "computed at 3 m", rows, scale=0.1, offset=0.0)
    depth = C.read_raster(C.WORK / "depression_depth_1m.tif")
    pond = C.read_raster(C.WORK / "pond_mask_1m.tif").astype(bool)
    cog("pond_depth_to_spill_3m_cm.tif", np.round(C.block_reduce(np.where(pond, depth, 0.0), 3, np.max) * 100).astype("uint16"),
        t3, 65535, "Depth below spill level inside the kept ponds and basins. Depth in m = value * 0.01",
        "m (stored as cm)", 3, "maximum", rows, scale=0.01, offset=0.0)
    cat = pd.DataFrame(rows)
    cat.to_csv(C.TAB / "exported_rasters.csv", index=False)
    print(cat[["file", "cell_size_m", "dtype", "size_mb"]].to_string(index=False))

    # ---- vectors ------------------------------------------------------------------------------
    out = C.VEC / "fort_reno_lidar_hydrology.gpkg"
    if out.exists():
        out.unlink()
    aoi = gpd.GeoDataFrame(dict(name=["Study area (working grid)"], area_km2=[g["area_km2"]]),
                           geometry=[box(g["xmin"], g["ymin"], g["xmax"], g["ymax"])], crs=C.DST_EPSG)
    aoi.to_file(out, layer="study_area", driver="GPKG")
    src = C.WORK / "hydro_vectors.gpkg"
    for layer in ("derived_streams", "catchments", "outlets", "depressions"):
        gdf = gpd.read_file(src, layer=layer)
        if layer == "catchments":
            gdf["geometry"] = gdf.simplify(1.0)
        if layer == "depressions":
            gdf["geometry"] = gdf.simplify(0.5)
        gdf.to_file(out, layer=layer, driver="GPKG")
    keep_fl = ["permanent_identifier", "gnis_name", "ftype", "fcode", "fcode_description", "geometry"]
    fl = C.nhd_layer("NHDFlowline")
    fl[[c for c in keep_fl if c in fl.columns]].to_file(out, layer="nhd_flowlines_clip", driver="GPKG")
    wb = C.nhd_layer("NHDWaterbody")
    pond_match = C.TAB / "nhd_waterbody_detection.csv"
    if pond_match.exists():
        wb["permanent_identifier"] = wb["permanent_identifier"].astype(str)
        wb = wb.merge(pd.read_csv(pond_match, dtype={"permanent_identifier": str})[["permanent_identifier", "area_in_study_m2", "depression_cover_frac",
                                               "detected"]], on="permanent_identifier", how="left")
    keep_wb = ["permanent_identifier", "gnis_name", "ftype", "fcode", "areasqkm", "area_in_study_m2",
               "depression_cover_frac", "detected", "geometry"]
    wb[[c for c in keep_wb if c in wb.columns]].to_file(out, layer="nhd_waterbodies_clip", driver="GPKG")
    import pyogrio
    layers = pyogrio.list_layers(out)
    print("GeoPackage layers:", ", ".join(l[0] for l in layers), f"({out.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
