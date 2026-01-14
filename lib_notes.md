# Library Notes for neurothresh Development

## neuroim2 Library (Primary Dependency)

### Package Overview
- **Version**: 0.8.3
- **Purpose**: Data structures and I/O for volumetric brain imaging (focus on fMRI)
- **Architecture**: S4 OOP with Rcpp for performance-critical operations
- **Core Dependencies**: Matrix, RNifti, Rcpp, RcppParallel

### Key Data Structures for neurothresh

#### ClusteredNeuroVol (CRITICAL)
The primary data structure for parcellations. Located in `R/clustervol.R` and `R/all_class.R` (lines 816-884).

**Slots:**
- `mask`: LogicalNeuroVol - spatial domain of clusters (TRUE = included)
- `clusters`: integer vector - cluster ID for each masked voxel
- `label_map`: named list - cluster ID → human-readable label mapping
- `cluster_map`: environment - cluster ID → voxel indices (O(1) lookup)
- Inherits `data` (sparseVector) from SparseNeuroVol

**Constructor:**
```r
ClusteredNeuroVol(mask, clusters, label_map=NULL, label="")
```

**Key Methods:**
- `num_clusters(x)` - number of clusters
- `centroids(x, type="center_of_mass"|"medoid")` - cluster centers
- `split_clusters(x)` - returns deflist of ROIVol objects
- `as.array(x)` - convert to dense 3D array

**Companion 4D Class: ClusteredNeuroVec**
- Aggregates 4D data by cluster using reduction function
- Stores T × K matrix (timepoints × clusters) instead of full voxel data
- Constructor: `ClusteredNeuroVec(x, cvol, FUN=mean, weights=NULL)`

---

### Volume Class Hierarchy

```
NeuroObj (base: has @space slot)
  └── NeuroVol (3D volumes)
      ├── DenseNeuroVol (dense array storage)
      │   └── LogicalNeuroVol (binary masks)
      ├── SparseNeuroVol (sparse vector storage)
      │   └── ClusteredNeuroVol (partitioned clusters)
      └── IndexLookupVol (index mapping for sparse access)
```

**DenseNeuroVol**: Full 3D array, efficient for random access
**SparseNeuroVol**: Matrix::sparseVector, memory-efficient for sparse data
**LogicalNeuroVol**: Binary masks, extends DenseNeuroVol

---

### NeuroSpace (Coordinate System)
Core spatial reference for all neuroimaging objects.

**Slots:**
- `dim`: integer - grid dimensions
- `origin`: numeric - world coordinates of first voxel
- `spacing`: numeric - voxel physical sizes (mm)
- `axes`: AxisSet - anatomical orientation
- `trans`: 4×4 matrix - voxel→world transformation
- `inverse`: 4×4 matrix - world→voxel transformation

**Key Conversion Functions:**
- `index_to_grid(x, idx)` - linear index → (i,j,k)
- `grid_to_index(x, coords)` - (i,j,k) → linear index
- `coord_to_grid(x, coords)` - world → voxel
- `grid_to_coord(x, coords)` - voxel → world

---

### 4D Time Series (NeuroVec)

```
NeuroVec (base 4D class)
  ├── DenseNeuroVec (dense 4D array)
  ├── SparseNeuroVec (sparse: timepoints × masked_voxels matrix)
  ├── MappedNeuroVec (memory-mapped)
  ├── FileBackedNeuroVec (on-demand file loading)
  ├── BigNeuroVec (bigstatsr::FBM for large data)
  └── NeuroVecSeq (sequence of variable-length NeuroVecs)
```

**Key Methods:**
- `series(x, i, j, k)` - extract time series at voxel
- `series(x, coords_matrix)` - extract multiple time series
- `sub_vector(x, timepoints)` - temporal subset
- `x[[t]]` - extract single timepoint as NeuroVol

**5D Data: NeuroHyperVec**
- Dimensions: [features × trials × voxels]
- For multi-basis analyses, spectral decompositions

---

### ROI (Region of Interest) Classes

```
ROI (virtual base)
  └── ROICoords (just coordinates)
      └── ROIVol (coords + values)
          └── ROIVolWindow (centered ROI, e.g., searchlight)
  └── ROIVec (4D: coords + time series matrix)
      └── ROIVecWindow (centered 4D ROI)
```

**ROI Creation Functions:**
- `spherical_roi(bvol, centroid, radius)` → ROIVolWindow
- `spherical_roi_set(bvol, centroids, radius)` → list of ROIs
- `cuboid_roi(bvol, centroid, surround)` → ROIVolWindow
- `square_roi(bvol, centroid, surround, fixdim)` → 2D ROI

**Key Methods:**
- `coords(x, real=FALSE)` - extract coordinates
- `indices(x)` - convert to linear indices
- `values(x)` - get data values
- `series_roi(vec, coords)` → ROIVec

---

### I/O Functions
- `read_vol(filename)` - read 3D NIfTI/AFNI
- `read_vec(filename, mode="mmap"|"bigvec")` - read 4D
- `write_vol(vol, filename)` - write 3D
- `write_vec(vec, filename)` - write 4D

---

### Useful Processing Functions
- `conn_comp(vol)` - connected component labeling
- `resample(vol, outspace, interp)` - spatial resampling
- `gaussian_blur(vol, sigma)` - smoothing
- `searchlight(vol, radius, fun)` - searchlight analysis

---

## neuroatlas Library (Secondary Dependency)

### Package Overview
- **Version**: 0.1.0
- **Purpose**: Unified interface for brain atlases and parcellations
- **Integration**: Returns ClusteredNeuroVol from neuroim2

### Available Atlases

| Atlas | Function | Regions | Notes |
|-------|----------|---------|-------|
| Schaefer | `get_schaefer_atlas(parcels, networks)` | 100-1000 | 7 or 17 networks |
| Glasser | `get_glasser_atlas()` | 360 | Multi-modal cortical |
| ASEG | `get_aseg_atlas()` | ~33 | FreeSurfer subcortical |
| Olsen MTL | `get_olsen_mtl()` | 16 | Medial temporal lobe |

**Convenience:** `sy_200_7()`, `sy_400_17()`, etc.

### Atlas Object Structure
```r
atlas <- list(
  name = "Schaefer-300-7networks",
  atlas = ClusteredNeuroVol(...),  # The volume
  ids = integer_vector,            # Region IDs
  labels = character_vector,       # Short names
  orig_labels = character_vector,  # Full original names
  hemi = c("left", "right", ...),  # Hemisphere
  network = c("Vis", "SomMot",...),# Network assignment (Schaefer)
  cmap = data.frame(red, green, blue),  # Colors 0-255
  roi_metadata = tibble(...)       # Tidy metadata
)
class(atlas) <- c("schaefer", "volatlas", "atlas")
```

### Key Functions

**Access & Metadata:**
- `roi_metadata(atlas)` → tibble with id, label, hemi, network, colors
- `filter_atlas(atlas, condition)` → filtered atlas
- `get_roi(atlas, label|id, hemi)` → list of ROIVol

**Analysis:**
- `reduce_atlas(atlas, data, FUN)` → region-wise summary stats
- `map_atlas(atlas, vals, thresh)` → map values to regions
- `merge_atlases(atlas1, atlas2)` → combined atlas
- `dilate_atlas(atlas, mask, radius)` → expand parcels

**Templates (via TemplateFlow):**
- `get_template(space, suffix)` - standard templates
- `get_template_brainmask()` - brain mask

---

## Key Integration Points for neurothresh

### Working with Parcellations
```r
# Load atlas
atlas <- neuroatlas::get_schaefer_atlas(parcels="200", networks="7")

# Extract ClusteredNeuroVol
cvol <- atlas$atlas

# Get number of parcels
k <- neuroim2::num_clusters(cvol)

# Split into ROIs
rois <- neuroim2::split_clusters(cvol)

# Extract parcel indices for parcel k
indices_k <- cvol@cluster_map[[as.character(k)]]

# Convert parcel to mask
parcel_mask <- neuroim2::as.logical(rois[[k]])
```

### Reducing 4D Data by Parcels
```r
# Load 4D fMRI
vec <- neuroim2::read_vec("fmri.nii")

# Create clustered representation (parcel-wise means)
cvec <- neuroim2::ClusteredNeuroVec(vec, cvol, FUN=mean)

# Access T × K matrix
ts_matrix <- as.matrix(cvec)
```

### Spatial Operations
```r
# Get parcel centroids
cents <- neuroim2::centroids(cvol, type="center_of_mass")

# Create spherical ROI around centroid
roi <- neuroim2::spherical_roi(cvol, cents[1,], radius=6)

# Extract time series from ROI
roi_ts <- neuroim2::series(vec, neuroim2::coords(roi))
```

---

## Summary for neurothresh Implementation

**Primary Classes:**
1. `ClusteredNeuroVol` - parcellation representation
2. `ClusteredNeuroVec` - parcel-aggregated 4D data
3. `NeuroSpace` - coordinate system management
4. `ROIVol` - individual region extraction

**Key Operations:**
- Parcel extraction via `split_clusters()` or `cluster_map` environment
- Time series extraction via `series()` or `ClusteredNeuroVec`
- Coordinate mapping via `grid_to_index()` / `index_to_grid()`
- Atlas loading via neuroatlas functions

**Performance Considerations:**
- ClusteredNeuroVol is sparse - memory efficient for parcellations
- `cluster_map` environment provides O(1) lookup by cluster ID
- Use `SparseNeuroVec` for masked 4D data to save memory
- C++ implementations available for heavy computation (Rcpp)
