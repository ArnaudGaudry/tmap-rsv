# tmap-rsv

TMAP visualization of RSV anti-viral compound annotations with specificity analysis across plant extracts.

## Overview

This project generates interactive TMAP (Tree-Map) visualizations for chemical compound annotations from RSV (Respiratory Syncytial Virus) anti-viral research. The workflow includes:

- **TMAP Generation**: Creates interactive tree-maps of chemical compounds using MAP4 fingerprints and LSH Forest layout
- **Specificity Analysis**: Identifies chemical classes specific to active plant extracts
- **Lotus Integration**: Compares annotated compounds against the Lotus natural products database
- **Interactive Visualizations**: HTML-based TMAP browser with multiple categorical overlays

## Installation

### Prerequisites
- Python 3.9
- Conda package manager

### Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/tmap-rsv.git
   cd tmap-rsv
   ```

2. **Create the conda environment**
   ```bash
   conda env create -f environment.yml
   conda activate tmap_rsv
   ```

3. **Install additional dependencies**
   ```bash
   pip install marimo upsetplot
   ```

## Data Requirements

The following data files are required in the `data/` directory:

- `annotations_rsv_subset.csv` - RSV compound annotations with SMILES and InChIKey
- `anti_rsv_data.csv` - Anti-viral activity data for plant extracts
- `canopus_rsv_subset.csv` - CANOPUS chemical classification annotations
- `20260206_lotus_all_taxa.csv.gz` - Lotus natural products database (provided separately)

The CSV files should contain the following key columns:
- `annotations_rsv_subset.csv`: `sample`, `ik2d`, `smiles`
- `anti_rsv_data.csv`: `Plaque`, `Well`, `AOAO_N1`, `AOAO_N2`, `A549_N1`, `A549_N2`, `vgf_code`, `Family`
- `canopus_rsv_subset.csv`: `sample`, `np_class`

## Usage

### Running the Marimo Notebook

The main analysis is conducted in an interactive Marimo notebook:

```bash
cd src
marimo run marimo_notebook.py
```

Or to open in the browser for interactive exploration:

```bash
marimo edit marimo_notebook.py
```

### Workflow Steps

The notebook executes the following steps in sequence:

#### 1. Data Loading and Preparation
- Loads RSV annotations and removes duplicates
- Loads plant extract anti-viral activity data
- Categorizes extracts as active (in cells, organoids, or both) or inactive
- Merges annotations with metadata and activity classification

#### 2. TMAP Generation for Annotated Compounds
- Calculates MAP4 fingerprints for all annotated RSV compounds
- Builds LSH Forest index for efficient layout computation
- Generates 2D coordinates using tree-layout algorithm
- Outputs coordinates and chemical descriptors to `data/260210_rsv_annnotations_tmap/`

#### 3. TMAP Generation for Lotus Database
- Processes the Lotus natural products database
- Generates TMAP for all 2D unique structures
- Outputs to `data/260208_lotus_tmap/`

#### 4. Interactive Visualization Generation
- Creates Faerun-based interactive HTML visualizations
- **RSV Annotations TMAP** (`260210_rsv_annnotations_tmap/tmap.html`):
  - Color overlays for each active plant species (Litsea polyantha, Ampelocissus arachnoidea, Clausena wallichii, Vepris macrophylla)
  - Continuous color scales for molecular descriptors: HAC, C-fraction, ring atom fraction, largest ring size
- **Lotus TMAP** (`260223_lotus_tmap/tmap.html`):
  - Highlights annotated vs. unannotated compounds
  - Molecular descriptor visualizations

#### 5. Chemical Class Specificity Analysis
- Groups compounds by chemical class (CANOPUS)
- Calculates specificity scores: ratio of compound counts in target extract vs. average in other extracts
- Generates scatter plots highlighting high-specificity chemical classes (score > 10)
- Outputs PNG specificity plot: `{species_name} specificity.png`

#### 6. Upset Plots
- Creates upset plots showing annotation overlap across activity categories:
  - Active in both cell and organoid models
  - Active only in A549 cells
  - Active only in AOAO organoids
  - Inactive in both models
- Outputs: `upsetplot.png`

## Output Files

### Interactive HTML Visualizations
- `data/260210_rsv_annnotations_tmap/tmap.html` - Main RSV annotation TMAP with species and descriptor overlays
- `data/260208_lotus_tmap/tmap.html` - Lotus database TMAP with annotation overlay

### Static Plots (PNG)
- `src/{species_name} specificity.png` - Specificity analysis for each plant species
- `src/upsetplot.png` - Annotation overlap across activity categories

### Intermediate Data
- `data/260210_rsv_annnotations_tmap/attribute.csv` - Descriptor data for RSV compounds
- `data/260210_rsv_annnotations_tmap/coordinates.dat` - Pickled coordinate data
- `data/260210_rsv_annnotations_tmap/lsh_forest.dat` - LSH Forest index

## Configuration

Key parameters can be adjusted in the notebook:

- **MAP4 fingerprint dimensions**: `map4_dimensions=1024` (default)
- **LSH Forest depth**: `lsh_forest_depth=64` (default)
- **TMAP layout parameters**: Node size, k-neighbors, scaling type, etc.
- **Specificity threshold**: `threshold=5` (minimum compound frequency for inclusion)
- **Species selection**: Modify the `mapping` dictionary to analyze different extracts

## Active Extracts

The following plant extracts showed anti-viral activity and are highlighted in the visualizations:

| Code | Species | Family |
|------|---------|--------|
| VGF157_C10 | *Litsea polyantha* | Lauraceae |
| VGF154_H04 | *Ampelocissus arachnoidea* | Vitaceae |
| VGF157_D08 | *Clausena wallichii* | Rutaceae |
| VGF154_B02 | *Vepris macrophylla* | Rutaceae |

## Helper Scripts

Additional utility scripts are available:

- `get_annotations.py` - Function to retrieve compound annotations
- `get_class_annotations.py` - Extract chemical class-specific annotations
- `get_superclass_annotations.py` - Extract chemical superclass-specific annotations
- `plot_structural_tmap.py` - Advanced TMAP visualization utilities

## Requirements

Key Python packages (see `environment.yml` for complete list):
- `marimo` - Interactive notebook framework
- `tmap` - TMAP layout and LSH Forest algorithms
- `map4` - MAP4 molecular fingerprints
- `faerun` - Interactive visualization library
- `rdkit` - Chemical informatics
- `pandas` - Data manipulation
- `altair` - Declarative visualization
- `upsetplot` - Upset plot visualization
- `scipy` - Scientific computing

## Citation

If you use this work, please cite:

[Add citation information]

## License

See LICENSE file for details.