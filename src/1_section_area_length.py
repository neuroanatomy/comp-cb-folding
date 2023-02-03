#---------------------------------------
# Compute all section areas and lengths
# 1 Feb 2023
#---------------------------------------

import os
import json
import pandas as pd
import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import Polygon, MultiPolygon

data = pd.read_csv("../data/raw/01_cb_data.csv")
scale = np.array(pd.read_csv("../data/raw/02_scale.csv")["Scale"])

def bbox(polys):
  '''Compute bounding box of a list of polygons
  Parameters
  ----------
  polys : list of shapely.geometry.Polygon
    List of polygons
  Returns
  -------
  list of float
    Bounding box of the list of polygons
  '''

  b = []
  for poly in polys:
    b = [*poly.bounds]
    if len(b) == 4:
      break
  if len(b) == 0:
    return []

  for poly in polys:
    bb = [*poly.bounds]
    if len(bb) == 0:
      continue
    b[0] = min(b[0], bb[0])
    b[1] = min(b[1], bb[1])
    b[2] = max(b[2], bb[2])
    b[3] = max(b[3], bb[3])
  return b

def convert_polygons_to_shapely_multipolygons(polygons):
  '''Convert a list of polygons to a shapely multipolygon
  Parameters
  ----------
  polygons : list of list of float
    List of polygons. The first polygon is the exterior, the rest are holes.
  Returns
  -------
  shapely.geometry.MultiPolygon
    Shapely multipolygon
  '''

  shapely_multipolygon = []
  for poly in [p[0] for p in polygons]:
    shapely_polygon = Polygon(poly[0], holes=poly[1:])
    shapely_multipolygon.append(shapely_polygon)
  return MultiPolygon(shapely_multipolygon)

# cerebellum
cb = []
for row in range(len(data)):
  try:
    if scale[row] == 0:
      print(row, "no scale")
      continue
    source = data.iloc[row]["URL"]
    dic = json.load(open(f"../data/raw/json/cb/{source.split('/')[-1]}.cb-50%.json", "r", encoding="utf-8"))
    cb_polys = dic["slice_polygons"]
    name = dic["name"]
    sm = convert_polygons_to_shapely_multipolygons(cb_polys)
    area = sm.area
    length = sm.length
    cb.append([name, np.log10(area), np.log10(length)])
  except BaseException as err:
    print(row, source, err)
    raise ValueError(err)
cb_df = pd.DataFrame({
  "Log10Area": [r[1] for r in cb],
  "Log10Length": [r[2] for r in cb]
})
# save as csv
cb_df.index = [r[0] for r in cb]
cb_df.to_csv(os.path.join("../data/derived/csv/01_cb_area_length.csv"))

# cerebrum
ctx = []
for row in range(len(data)):
  try:
    if scale[row] == 0:
      print(row, "no scale")
      continue
    source = data.iloc[row]["URL"]
    dic = json.load(open(f"../data/raw/json/ctx/{source.split('/')[-1]}.ctx-50%.json", "r", encoding="utf-8"))
    ctx_polys = dic["slice_polygons"]
    name = dic["name"]
    sm = convert_polygons_to_shapely_multipolygons(ctx_polys)
    area = sm.area
    length = sm.length
    ctx.append([name, np.log10(area), np.log10(length)])
  except BaseException as err:
    print(row, source, err)
    raise ValueError(err)
ctx_df = pd.DataFrame({
  "Log10Area": [r[1] for r in ctx],
  "Log10Length": [r[2] for r in ctx]
})
# save as csv
ctx_df.index = [r[0] for r in ctx]
ctx_df.to_csv(os.path.join("../data/derived/csv/02_ctx_area_length.csv"))
