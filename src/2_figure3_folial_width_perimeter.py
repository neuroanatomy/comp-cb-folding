#------------------------------------------
# Compute all folial widths and perimeters
# 1 Feb 2023
#------------------------------------------

'''Compute and save gyral measurements for all subjects in the dataset'''

import json
import numpy as np
import pandas as pd
import gyri as gy
from convert_polygons_to_shapely_multipolygons import convert_polygons_to_shapely_multipolygons

data = pd.read_csv("../data/raw/01_cb_data.csv")
scale = np.array(pd.read_csv("../data/raw/02_scale.csv")["Scale"])

def compute_gyral_measurements_for_all():
  '''compute gyral measurements for all subjects in the dataset'''

  results_period = []
  results_width = []
  for row in range(len(data)):
    # get scale. Skip subject if scale is unavailable
    if scale[row] == 0:
      continue

    source = data.iloc[row]["URL"]
    dic = json.load(
      open(f"../data/raw/json/cb/{source.split('/')[-1]}.cb-50%.json",
      "r", encoding="utf-8"))
    cb_polys = dic["slice_polygons"]
    name = dic["name"]
    cb_mid = convert_polygons_to_shapely_multipolygons(cb_polys)

    polys, min_length = gy.resample_cerebellum_contour(cb_mid)
    labels, sulci_index, _, _ = gy.label_contour(polys)
    labels = gy.filter_sulci(labels, sulci_index)

    period = gy.compute_gyral_period(labels, 1, min_length)
    width = gy.compute_gyral_width(labels, 1, polys)

    results_period.append((name, np.median(period), np.mean(period), np.std(period)))
    results_width.append((name, np.median(width), np.mean(width), np.std(width)))

    print(row, name, np.median(period), np.median(width))

  # save results
  data_frame = pd.DataFrame(
    [col[1:] for col in results_width],
    index=[col[0] for col in results_width],
    columns = ["WidthMedian", "WidthMean", "WidthStd"])
  data_frame.to_csv("../data/derived/csv/03_cb_width.csv")

  data_frame = pd.DataFrame(
    [col[1:] for col in results_period],
    index=[col[0] for col in results_period],
    columns = ["PeriodMedian", "PeriodMean", "PeriodStd"])
  data_frame.to_csv("../data/derived/csv/04_cb_period.csv")


compute_gyral_measurements_for_all()
