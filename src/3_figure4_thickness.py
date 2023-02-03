#----------------------------------------
# Compute all molecular layer thicknesses
# 1 February 2023
#----------------------------------------

import json
import numpy as np
import pandas as pd
import thickness as th
from convert_polygons_to_shapely_multipolygons import convert_polygons_to_shapely_multipolygons

data = pd.read_csv("../data/raw/01_cb_data.csv")
scale = np.array(pd.read_csv("../data/raw/02_scale.csv")["Scale"])

def compute_all_thicknesses():
  '''compute thickness of the molecular layer for all subjects'''

  # loop to process all subjects
  with open("../data/derived/csv/05_cb_thickness.csv", "w", encoding="utf-8") as file:
    file.write(",ThicknessMedian,ThicknessMean,ThicknessStd\n")
    for row in range(len(data)):
      name = data.iloc[row]["Name"]
      source = data.iloc[row]["URL"]
      print(row, name, source)

      if scale[row] == 0:
        continue

      dic = json.load(
        open(f"../data/raw/json/cb/{source.split('/')[-1]}.cb-50%.json",
        "r", encoding="utf-8"))
      cb_polys = dic["slice_polygons"]
      name = dic["name"]
      cb_mid = convert_polygons_to_shapely_multipolygons(cb_polys)

      img_path = "../data/raw/img/cb/" + source.split("/")[-1] + ".cb-50%.png"

      result_data = th.compute_thickness(
          scale[row], cb_mid, name,
          img_path
      )
      if result_data:
        file.write(result_data)

compute_all_thicknesses()
