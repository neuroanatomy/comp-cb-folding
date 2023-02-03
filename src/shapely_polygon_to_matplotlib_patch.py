import numpy as np
from matplotlib.patches import PathPatch
from matplotlib.path import Path

def polygon_patch(polygon, **kwargs):
  '''
  Make a matplotlib patch from a shapely polygon
  replaces descartes
  adapted from https://github.com/geopandas/geopandas/issues/1039#issuecomment-748625852
  Parameters
  ----------
  polygon : shapely.geometry.polygon.Polygon
    polygon to convert
  Returns
  -------
  patch : matplotlib.patches.PathPatch
    matplotlib patch
  '''
  patch = PathPatch(
    Path.make_compound_path(
      Path(np.asarray(polygon.exterior.coords)[:, :2]),
      *[Path(np.asarray(ring.coords)[:, :2]) for ring in polygon.interiors]
    ), **kwargs)
  return patch
