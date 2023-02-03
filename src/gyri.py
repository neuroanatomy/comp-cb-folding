'''compute gyral measurements for a contour'''

import numpy as np
import shapely
from shapely import affinity
from shapely.geometry import LinearRing
from sklearn.cluster import KMeans
from skimage.measure import subdivide_polygon
from scipy.signal import find_peaks

def smooth_polygon(poly, iters=1):
  '''Smooth polygon by averaging with its neighbours
  a certain number of iterations
  Parameters
  ----------
  poly : np.array
    polygon coordinates
  iters : int
    number of iterations
  Returns
  -------
  poly : np.array
    smoothed polygon coordinates
  '''
  for _ in range(iters):
    poly = (poly + np.roll(poly, -1, axis=0) + np.roll(poly, 1, axis=0))/3
  return poly

def curvature_features(poly):
  '''Compute curvature measures: as 2nd derivative of poly, as cross product of
  tangent vectors, and a signed version of 2nd derivative
  Parameters
  ----------
  poly : np.array
    polygon coordinates
  Returns
  -------
  curv : np.array
    curvature measures
  scurv : np.array
    signed curvature measures
  cross : np.array
    cross product of tangent vectors
  '''
  deltap = np.diff(poly, append=[poly[0]], axis=0)
  ndeltap = deltap/np.sqrt(np.sum(deltap*deltap, axis=1))[:, np.newaxis]
  deltap2 = 0.5 * (deltap - np.roll(deltap, 1, axis=0))
  curv = np.sqrt(np.sum(deltap2*deltap2, axis=1))
  cross = np.cross(ndeltap, deltap2)
  scurv = (curv * np.sign(cross))
  return curv, scurv, cross

def resample_cerebellum_contour(cb_mid_row, n_points=None):
  '''Resample cerebellum contour taking into account its scaling with size
  Parameters
  ----------
  cb_mid_row : shapely.geometry.multipolygon.MultiPolygon
    cerebellum contour
  n_points : int
    number of points to resample to
  Returns
  -------
  poly2 : np.array
    resampled polygon coordinates
  n_points : int
    number of points to resample to
  '''

  scaled_mpoly = affinity.scale(cb_mid_row, xfact=1, yfact=1, origin=(0,0))
  print("Polygon length:", scaled_mpoly.length)

  if isinstance(scaled_mpoly, shapely.geometry.polygon.Polygon):
    poly = np.array(scaled_mpoly.exterior.coords)
  else:
    p_index = np.argmax([np.array(p.exterior.coords).shape[0] for p in scaled_mpoly.geoms])
    poly = np.array(scaled_mpoly.geoms[p_index].exterior.coords)
  poly2 = subdivide_polygon(poly, degree=3)
  ring = LinearRing(poly2)
  ring_length = ring.length
  print("Ring length:", ring_length)

  if n_points:
    min_length = ring_length/n_points
  else:
    better_min_length = -3 + (np.log10(scaled_mpoly.length)+0.56)/3.1
    min_length = 10**better_min_length
    n_points = int(np.ceil(ring_length/min_length))

  polys = np.zeros((n_points, 2))
  for i in range(n_points):
    polys[i, :] = np.array(ring.interpolate(i/n_points, normalized=True).coords)

  return polys, min_length

def label_contour(polys, iters=10):
  '''Label contour vertices in 3 classes: sulci, gyri and wall
  Parameters
  ----------
  polys : np.array
    polygon coordinates
  iters : int
    number of iterations for smoothing
  Returns
  -------
  labels : np.array
    labels for each vertex
  sulci_index : int
    index of sulci label
  gyri_index : int
    index of gyri label
  wall_index : int
    index of wall label
  '''
  # compute features
  features = np.zeros((len(polys), 3*(iters+1)))
  features[:, 0], features[:, 1], features[:, 2] = curvature_features(polys)
  for i in range(1, iters):
    polys1 = smooth_polygon(polys, 10*2**i)
    features[:, 3*i], features[:, 3*i+1], features[:, 3*i+2] = curvature_features(polys1)
  kmclustering = KMeans(n_clusters=3)
  kmclustering.fit(features)

  # label vertices with kmeans clustering
  labels = kmclustering.labels_
  sulci_index = np.argmax([np.mean(features[labels==label,2]) for label in [0, 1, 2]])
  gyri_index = np.argmin([np.mean(features[labels==label,2]) for label in [0, 1, 2]])
  tmp = [0, 1, 2]
  tmp.remove(sulci_index)
  tmp.remove(gyri_index)
  wall_index = tmp[0]

  return labels, sulci_index, gyri_index, wall_index

def filter_sulci(labels, sulci_index):
  '''Filter sulci labels to keep only the peaks
  Parameters
  ----------
  labels : np.array
    labels for each vertex
  sulci_index : int
    index of sulci label
  Returns
  -------
  filtered_labels : np.array
    filtered labels for each vertex
  '''
  labels = np.array([1 if x==sulci_index else 0 for x in labels])
  filtered_labels = labels + np.roll(labels, -1) + np.roll(labels, 1)
  filtered_labels = smooth_polygon(filtered_labels, iters=100)
  peaks = find_peaks(filtered_labels)
  filtered_labels[:] = 0
  filtered_labels[peaks[0]] = 1
  return filtered_labels

def compute_gyral_period(labels, sulci_index, min_length):
  '''Compute the period of gyri, defined as the length along the
  contour from one sulcus to the next
  Parameters
  ----------
  labels : np.array
    labels for each vertex
  sulci_index : int
    index of sulci label
  min_length : float
    minimum length between two vertices
  Returns
  -------
  period : list
    list of gyral periods
  '''
  pos = np.argmax(labels == sulci_index)
  pos += np.argmax(labels[pos:] != sulci_index)
  period = []
  while True:
    pos0 = pos
    pos1 = np.argmax(labels[pos:] == sulci_index)
    pos2 = np.argmax(labels[pos+pos1:] != sulci_index)
    pos += pos1 + pos2
    period.append((pos1+pos2) * min_length)
    if pos == pos0:
      break
    pos0 = pos
  return period

def compute_gyral_width(labels, sulci_index, polys):
  '''Compute the width of gyri, defined as the euclidean distance
  from one sulcus to the next
  Parameters
  ----------
  labels : np.array
    labels for each vertex
  sulci_index : int
    index of sulci label
  polys : np.array
    polygon coordinates
  Returns
  -------
  width : list
    list of gyral widths
  '''
  pos = np.argmax(labels == sulci_index)
  pos += np.argmax(labels[pos:] != sulci_index)
  width = []
  while True:
    pos0 = pos
    pos1 = np.argmax(labels[pos:] == sulci_index)
    pos2 = np.argmax(labels[pos+pos1:] != sulci_index)
    width.append(np.linalg.norm(polys[pos + pos1+pos2] - polys[pos0]))
    pos += pos1 + pos2
    if pos == pos0:
      break
    pos0 = pos
  return width
