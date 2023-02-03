'''
Compute the thickness of the molecular layer of the cerebellum
'''

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import shapely
import skimage
from scipy.interpolate import interp2d
from scipy.signal import find_peaks
from shapely import affinity
from shapely.geometry import LinearRing, LineString
from skimage import exposure
from skimage import io, filters
from skimage.measure import subdivide_polygon
from shapely_polygon_to_matplotlib_patch import polygon_patch

def mid_vector(v1, v2):
  '''compute the angle and the vector between two vectors'''
  v1 = v1/np.linalg.norm(v1)
  v2 = v2/np.linalg.norm(v2)
  v1t = np.array([-v1[1], v1[0]])
  x = v2.dot(v1)
  y = v2.dot(v1t)
  ang12 = np.arctan2(y, x)
  r = v1*np.cos(ang12/2) + v1t*np.sin(ang12/2)
  return ang12, r

def polygon_normals(poly, profile_length = 5):
  '''compute normal vectors for each point in the polygon'''
  normals = np.zeros(poly.shape)
  for i, b in enumerate(poly):
    a = poly[(i-1+len(poly))%len(poly)]
    c = poly[(i+1)%len(poly)]
    x1 = a[1] - b[1]
    y1 = -(a[0] - b[0])
    x2 = -(c[1] - b[1])
    y2 = (c[0] - b[0])
    _, normal = mid_vector(np.array([x1, y1]), np.array([x2, y2]))
    x, y = profile_length * normal
    normals[i] = [b[0] - x, b[1] - y]
  return normals

def polygon_resample(poly, normals, profile_length=5, max_ang = np.pi/6):
  '''improve sampling by comparing the angle between consecutive normal vectors'''

  where = []
  what_start = []
  what_end = []
  for i, _ in enumerate(poly):
    prev = (i - 1 + len(poly))%len(poly)
    ang, new_normal = mid_vector(normals[prev]-poly[prev], normals[i]-poly[i])
    if np.abs(ang) > max_ang:
      new_poly_point = (poly[i] + poly[prev])/2
      where.append(i)
      what_start.append(new_poly_point)
      what_end.append(new_poly_point + profile_length * new_normal)
  if where:
    return (
      np.insert(poly, where, what_start, axis=0),
      np.insert(normals, where, what_end, axis=0)
    )
  return poly, normals

def load_image(img_path):
  '''Load image and apply preprocessing'''
  img = io.imread(img_path)
  if len(img.shape) == 2:
    img_gray = img
  else:
    img_gray = skimage.color.rgb2gray(img)
  img_gray = img_gray/np.max(img_gray)
  img = filters.median(img_gray, skimage.morphology.disk(1))
  img = exposure.equalize_adapthist(img)
  return img

def scale_contour_to_image(img, sub_scale, sub_cb_mid):
  '''Scale the contour to the image size'''
  g = (img.shape[1]/1000)/sub_scale
  return affinity.scale(sub_cb_mid, xfact=g, yfact=g, origin=(0,0))

def resample_contour(sm, img, min_length):
  '''Resample the contour to a given minimum length between vertices'''
  p = np.array(sm.exterior.coords)
  p2 = subdivide_polygon(p, degree=3)
  ring = LinearRing(p2)
  ring_length = ring.length * (1000/img.shape[1])
  n_points = int(np.ceil(ring_length/min_length))
  pp = np.zeros((n_points, 2))
  for i in range(n_points):
    pp[i, :] = np.array(*ring.interpolate(i/n_points, normalized=True).coords)
  return pp

def make_mask(img, scaled_mpoly):
  '''Make a mask from the contour'''
  fig = plt.figure(figsize=(img.shape[1],img.shape[0]), dpi=1)
  ax = fig.add_subplot()
  plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, hspace = 0, wspace = 0)
  ax.set_axis_off()
  ax.set_xlim(0, img.shape[1])
  ax.set_ylim(img.shape[0], 0)

  for sm in scaled_mpoly:
    patch = polygon_patch(sm, alpha=1, ec="r", fc="k", lw=2)
    ax.add_patch(patch)
    patch = polygon_patch(sm, alpha=1, ec="w", fill=False, lw=200)
    ax.add_patch(patch)

  fig.canvas.draw()
  fig.canvas.figure.set_facecolor("white")
  mask = fig.canvas.buffer_rgba()
  mask = np.array(mask)
  plt.close()
  mask = np.invert(mask)[:,:,0]
  return mask

def compute_image_gradients(img, mask, sigma=3, smooth_iterations=10):
  '''Compute image gradients'''
  smo = img.copy()
  for _ in range(smooth_iterations):
    smo[mask==0] = 1
    smo = skimage.filters.gaussian(smo, sigma=3)
  DxW, DyW = np.gradient(smo, 1)

  return DxW, DyW, smo

def interp_functions(img, smo, DxW, DyW):
  '''Interpolate function for image and gradients'''
  X, Y = np.arange(0, img.shape[1], 1), np.arange(0, img.shape[0], 1)
  fni = interp2d(X, Y, img)
  fng = interp2d(X, Y, smo)
  fnx = interp2d(X, Y, DyW)
  fny = interp2d(X, Y, DxW)
  return fni, fng, fnx, fny

def get_profile_lines(pp, fng, fnx, fny, profile_length=30, data_steps=20, total_steps=40, step_length=0.5):
  '''Compute profile lines'''
  end = polygon_normals(pp, profile_length=profile_length)
  qq = pp.copy()

  # 0.25 px inset
  for i, start in enumerate(qq):
    norm = np.sum((end[i] - start)**2)**0.5

    inset = 0.25/norm

    qq[i] = start*(1-inset) + end[i]*inset

  # profile lines
  profile_indices = []
  profile_lines = []
  for index, q_new in enumerate(qq):
    n_samples = 0
    xy = np.zeros((total_steps - 1, 2))
    val0 = fng(*q_new)[0]
    for _ in range(data_steps):
      d = np.array([fnx(*q_new)[0], fny(*q_new)[0]])
      d = step_length * d / np.linalg.norm(d)
      q_new = np.array([q_new[0]-d[0], q_new[1]-d[1]])
      val = fng(*q_new)[0]
      if val > val0:
        break
      val0 = val
      xy[n_samples, :] = q_new
      n_samples += 1
    if n_samples < 3:
      continue

    # extend
    x, y = xy[n_samples-1]
    x0, y0 = xy[n_samples-2]
    new_steps = total_steps - n_samples - 1
    for i in range(2, new_steps + 2):
      xy[n_samples+i-2] = (x-x0)*i + x0, (y-y0)*i + y0

    # add back start point
    xy = np.array([pp[index], *xy])

    # append to the list of profiles
    profile_lines.append(xy)
    profile_indices.append(index)

  return profile_lines, profile_indices

def get_profile_levels(profile_lines, fni, total_steps):
  '''Compute grey levels along profiles'''
  profile_levels = []
  for xy in profile_lines:
    gr = np.zeros(total_steps)
    for i,(x,y) in enumerate(xy):
      gr[i] = fni(x,y)[0]
    profile_levels.append(gr)
  return profile_levels

def extract_image_profiles(img, smo, DxW, DyW, pp):
  '''
  fni: interpolating function for image
  fng: interpolating function for gaussian smoothed image
  fnx: interpolating function for x gradient
  fny: interpolating function for y gradient
  '''
  total_steps = 40
  fni, fng, fnx, fny = interp_functions(img, smo, DxW, DyW)
  profile_lines, profile_indices = get_profile_lines(pp, fng, fnx, fny, total_steps=total_steps)
  profile_levels = get_profile_levels(profile_lines, fni, total_steps)
  return profile_lines, profile_levels, profile_indices

def molecular_layer_thickness(profile_levels):
  '''
  Returns an array of thicknesses, and an array of indices
  for the profiles to which those thicknesses correspond.
  The thickness value is the number of profile samples
  from the beginning of the profile to the point where
  the molecular layer boundary occurs.
  '''
  ind = []
  th = []
  for i, prl in enumerate(profile_levels):
    gprl = np.gradient(prl)
    jj = find_peaks(-gprl, prominence=0.05)[0]
    if len(jj) == 0:
      continue
    j = jj[0]
    if prl[j]>0.5 and len(jj)>1:
      j = jj[1]
    ind.append(i)
    th.append(j)

  return np.array(th), ind

#--------------------------------------------------------#
#                                                        #
#                       Main script                      #
#                                                        #
#--------------------------------------------------------#

def compute_thickness(
  scale_row,
  cb_mid_row,
  name,
  img_path=None
):
  '''Compute thickness of the molecular layer
  from the image and the cerebellum contour
  Parameters
  ----------
  scale_row : int
    scale of the image
  cb_mid_row : int
    cerebellum contour
  name : str
    name of the image
  img_path : str
    path to the image
  Returns
  -------
  thickness : float
    thickness of the molecular layer
  '''

  # set matplotlib background (required for making the mask)
  matplotlib.use('agg')
  if os.path.exists(img_path):
    img = load_image(img_path)
  else:
    print("WARNING: No image file at path", img_path)
    return

  if scale_row == 0:
    print("No scale. Skipping")
    return None

  try:
    # scale the contour to fit the dimensions of the image
    scaled_mpoly = scale_contour_to_image(img, scale_row, cb_mid_row)
  except BaseException as err:
    print("ERR2:", err)
    return None

  # length in svg dimensions
  scaled_mpoly_length = scaled_mpoly.length * (1000/img.shape[1])
  min_length = scaled_mpoly_length**0.3/10

  # compute cb mask

  # form an array of polygons
  if isinstance(scaled_mpoly, shapely.geometry.Polygon):
    scaled_mpoly = [scaled_mpoly]
  else:
    scaled_mpoly = scaled_mpoly.geoms
  mask = make_mask(img, scaled_mpoly)

  # compute image gradients
  DxW, DyW, smo = compute_image_gradients(img, mask)

  # compute profiles
  profile_lines = []
  profile_levels = []
  profile_indices = []
  pps = []
  for sm in scaled_mpoly:
    try:
      pp = resample_contour(sm, img, min_length)
    except BaseException as err:
      print("ERR3:", err)
      continue

    # compute normal direction profiles
    profile_length = 50
    end = polygon_normals(pp, profile_length=profile_length)
    pps.append(pp)

    # extract grey level profiles
    ind0 = 0
    plin, plev, pind = extract_image_profiles(img, smo, DxW, DyW, pp)
    profile_lines.extend(plin)
    profile_levels.extend(plev)
    profile_indices.extend(np.array(pind) + ind0)
    ind0 += len(pind)

  # estimate molecular layer thickness
  th, ind = molecular_layer_thickness(profile_levels)

  thickness_array = []
  try:
    for i, _ in enumerate(th):
      xy = profile_lines[ind[i]]
      total_length = LineString(xy).length # length in image dimensions (px)
      thickness_array.append(total_length*th/len(xy))
  except BaseException as err:
    print("ERROR:", err)
    return

  print(
    np.median(thickness_array) * (1000/img.shape[1]) * scale_row,
    np.mean(thickness_array) * (1000/img.shape[1]) * scale_row,
    np.std(thickness_array) * (1000/img.shape[1]) * scale_row
  )

  csv = "%s,%g,%g,%g\n"%(
    name,
    np.median(thickness_array) * (1000/img.shape[1]) * scale_row,
    np.mean(thickness_array) * (1000/img.shape[1]) * scale_row,
    np.std(thickness_array) * (1000/img.shape[1]) * scale_row
  )
  return csv
