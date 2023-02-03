'''Convert a list of polygons to a shapely multipolygon'''

from shapely.geometry import Polygon, MultiPolygon

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
