# Run guide

These steps rerun the whole project on a new computer. I ran it on macOS (Apple silicon) with Python 3.12. It should also work on Linux and Windows, because it uses only Python packages.

## 1. What you need

- Python 3.12 (3.11 should also work).
- About 1.5 GB of free disk space: about 0.65 GB of downloads and about 0.4 GB of working rasters.
- About 8 GB of memory. The largest step holds a few 13.6 million cell arrays at once.
- An internet connection for step 1 only.
- No accounts, no API keys, no GIS software.

## 2. Install

```bash
git clone https://github.com/anamikadristi/lidar-rangeland-terrain-hydrology-fort-reno.git
cd lidar-rangeland-terrain-hydrology-fort-reno
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Notes:

- `lazrs` is the LAZ decompression backend for `laspy`. If `laspy` says it has no LAZ backend, run `pip install "laspy[lazrs]"`.
- `pysheds` 0.5 calls `numpy.in1d`, which newer numpy versions removed. `scripts/common.py` adds it back as an alias of `numpy.isin`, so nothing has to be patched by hand.
- If pip cannot build `rasterio` or `pyogrio` on your system, install the same packages with conda from the conda-forge channel.

## 3. Run everything

```bash
bash scripts/run_all.sh
```

This runs the eight numbered scripts in order. To reuse downloads that are already in `data/raw/`:

```bash
SKIP_DOWNLOAD=1 bash scripts/run_all.sh
```

To use a specific interpreter: `PYTHON=/path/to/python bash scripts/run_all.sh`.

To keep the raw downloads or the working rasters somewhere else, set `LRTH_RAW_DIR` and `LRTH_WORK_DIR` to folders of your choice before running.

## 4. What each step does and how long it takes

Times are from my run on a laptop with 11 CPU cores and 36 GB of memory.

| Step | Script | What it does | Time |
|---|---|---|---|
| 1 | `01_download_data.py` | Queries The National Map API, downloads 6 LAZ tiles (364 MB), the USGS 1 m DEM tile (281 MB), the NHD GeoPackage (12 MB zipped) and 4 metadata files | 15 to 30 minutes |
| 2 | `02_point_cloud_qa.py` | Defines the working grid, counts classes and returns, maps density | under 1 minute |
| 3 | `03_dtm_dsm_chm.py` | Grids the DTM, DSM, CHM and woody cover | about 1 minute |
| 4 | `04_accuracy_vs_usgs_dem.py` | Compares my DTM with the USGS 1 m DEM | under 1 minute |
| 5 | `05_terrain_derivatives.py` | Slope, aspect, hillshade, curvature | under 1 minute |
| 6 | `06_hydrology.py` | Depressions, breaching, flow direction, accumulation, streams, catchments, TWI | 2 to 3 minutes |
| 7 | `07_compare_nhd_interpret.py` | NHD comparison, woody cover against drainage | under 1 minute |
| 8 | `08_export_products.py` | Small GeoTIFFs and the GeoPackage | under 1 minute |

About the download: the USGS LiDAR server (`rockyweb.usgs.gov`) gave me only about 80 to 100 kB/s per connection and answered HTTP 429 when I opened more than about three connections. The script downloads three tiles at a time, waits 30 seconds after a refusal, and resumes partial files. If it stops, run it again. Finished files are skipped.

## 5. How to check the result

After the run, compare your numbers with mine. Small differences in the last digit are possible on other systems.

| Check | File | My value |
|---|---|---|
| Points in the six tiles | `outputs/tables/qa_tile_summary.csv` (sum of `points`) | 33,999,458 |
| Working grid | `data/grid_definition.json` | 3000 by 4539 cells, 13.617 km2 |
| Pulse density | `outputs/tables/qa_density_summary.csv` | 2.36 first returns per m2 |
| DTM against USGS DEM, land cells | `outputs/tables/accuracy_dtm_vs_usgs_dem.csv`, row "Land cells (outside NHD waterbodies)" | bias +0.001 m, RMSE 0.026 m |
| Woody cover above 2 m | `outputs/tables/vegetation_cover_summary.csv` | 34.05 ha, 2.50 % |
| Largest outlet | `outputs/tables/outlets_catchments.csv` | 9.584 km2 at the north edge |
| Derived stream length at 5 ha | `outputs/tables/hydrology_summary.csv` | 46.03 km |
| Ponds and basins kept | `outputs/tables/hydrology_summary.csv` | 25, with 365,347 m3 to spill level |

Then open the eight figures in `outputs/figures/` and compare them with the ones in the README.

## 6. Looking at the outputs in a GIS

Everything is in NAD83 / UTM zone 14N (EPSG:26914) with heights in metres above NAVD88.

1. In QGIS or ArcGIS Pro, add the rasters from `outputs/rasters/`. They are Cloud-Optimized GeoTIFFs stored as scaled integers. The scale and offset are written in each file, and `outputs/tables/exported_rasters.csv` lists them. For example, elevation in metres is `dtm_2m_cm.tif` value times 0.01 plus 400.
2. Add the layers of `outputs/vectors/fort_reno_lidar_hydrology.gpkg`: `study_area`, `derived_streams`, `catchments`, `outlets`, `depressions`, `nhd_flowlines_clip`, `nhd_waterbodies_clip`.
3. The full 1 m rasters (DTM, DSM, CHM, slope, flow direction, flow accumulation and others) are in `data/work/` after a run. They are not in the repository because they are too large.

I did not use QGIS or ArcGIS Pro for this project. All processing and all figures were made in Python.

## 7. Changing the study area

The six tiles are listed in `scripts/common.py` (`TILES`). To move the study area, change that list, and change `DEM_URL` and `NHD_URL` if the new tiles fall in another 1 m DEM tile or another HUC8. `data/tnm_product_index.csv` lists the tiles that The National Map API returned around this site. Other LiDAR projects may use other units and coordinate systems, so check `SRC_EPSG` and the feet-to-metre conversion in `read_tile` before using the scripts elsewhere.
