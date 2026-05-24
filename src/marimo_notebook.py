import marimo

__generated_with = "0.17.6"
app = marimo.App(width="medium")


@app.cell
def _():
    import os
    import pickle

    import marimo as mo
    import pandas as pd
    import polar as pl
    import tmap as tm
    import numpy as np
    import scipy.stats as ss
    import altair as alt
    import matplotlib.pyplot as plt

    from map4 import MAP4Calculator
    from rdkit.Chem import AllChem
    from rdkit import Chem
    from faerun import Faerun
    from tqdm import tqdm
    from pathlib import Path
    import matplotlib.colors as mcolors
    return (
        AllChem,
        Chem,
        Faerun,
        MAP4Calculator,
        Path,
        alt,
        mcolors,
        mo,
        np,
        os,
        pd,
        pickle,
        plt,
        ss,
        tm,
        tqdm,
    )


@app.cell
def _():
    from get_annotations import get_annotations
    return


@app.cell
def _(Counter):
    def Top_N_classes(N=0, input_data=[], input_labels=[]):
        """Keep only the top N classes for the input classes and replace the following by 'Other', default N = 0 return input"""
        if N == 0:
            return input_data, input_labels
        else:
            top_input = [i for i, _ in Counter(input_data).most_common(N)]
            output_labels = [(7, "Other")]
            input_map = [7] * len(input_data)
            value = 1
            for i, name in input_labels:
                if i in top_input:
                    v = value
                    if v == 7:
                        v = 0
                    output_labels.append((v, name))
                    input_map[i] = v
                    value += 1
            output_data = [input_map[val] for _, val in enumerate(input_data)]
            return output_data, output_labels
    return (Top_N_classes,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # TMAP Annotations

    Get TMAP from annotations
    """)
    return


@app.cell
def _(pd):
    #annotations = get_annotations()
    annotations = pd.read_csv('../data/annotations_rsv_subset.csv', sep=';')
    annotations_no_duplicates = annotations.drop_duplicates(subset='ik2d')
    annotations_no_duplicates = annotations_no_duplicates[annotations_no_duplicates['smiles'] != 'unknown']

    print(f"""
    Number of samples: {annotations['sample'].nunique()}
    Number of annotated features: {len(annotations)}
    Number of unique structures: {len(annotations_no_duplicates)}
    """)

    annotations.head()
    return annotations, annotations_no_duplicates


@app.cell
def _(AllChem, MAP4Calculator, Path, os, pd, pickle, tm, tqdm):
    def generate_tmap(
            df: pd.DataFrame,
            smiles_column: str,
            ik2d_column: str,
            output_dir: str,
            map4_dimensions: int = 1024,
            lsh_forest_depth: int = 64,
        ) -> None:
        """
        Generate TMAP visualization from chemical structures.

        Args:
            df: Input dataframe with chemical structures
            smiles_column: Name of column containing SMILES strings
            ik2d_column: Name of column containing InChIKey 2D
            output_dir: Output directory name (saved to ../data/{output_dir})
            map4_dimensions: MAP4 fingerprint dimensions
            lsh_forest_depth: LSH Forest depth parameter
        """

        output_path = Path("../data") / output_dir
        output_path.mkdir(parents=True, exist_ok=True)
        lf_store_path = os.path.join(output_path, "lsh_forest.dat")
        set_attributes_path = os.path.join(output_path, "attribute.csv")
        set_coordinates_path = os.path.join(output_path, "coordinates.dat")

        # Calculate fingerprints and descriptors
        fps = []
        descriptors = {"ik2d":[], "smiles":[], "hac": [], "c_frac": [], "ring_atom_frac": [], "largest_ring_size": []}

        for _, row in tqdm(df.iterrows()):
            mol = AllChem.MolFromSmiles(row[smiles_column])
            if mol is None:
                continue

            # Calculate descriptors
            atoms = mol.GetAtoms()
            size = mol.GetNumHeavyAtoms()
            n_c = sum(1 for atom in atoms if atom.GetSymbol().lower() == "c")
            n_ring_atoms = sum(1 for atom in atoms if atom.IsInRing())

            descriptors["ik2d"].append(row[ik2d_column])
            descriptors["smiles"].append(row[smiles_column])
            descriptors["hac"].append(size)
            descriptors["c_frac"].append(n_c / size if size > 0 else 0)
            descriptors["ring_atom_frac"].append(n_ring_atoms / size if size > 0 else 0)

            sssr = AllChem.GetSymmSSSR(mol)
            descriptors["largest_ring_size"].append(max([len(s) for s in sssr]) if sssr else 0)
            fps.append(mol)

        print("Parse smiles and get descriptors: done")
        # Build LSH Forest and TMAP
        map4 = MAP4Calculator(dimensions=map4_dimensions)
        lf = tm.LSHForest(map4_dimensions, lsh_forest_depth)

        fps = map4.calculate_many(fps)
        print("Get fingerprints: done")
        lf.batch_add(fps)
        lf.index()
        print("Get LSH forest: done")

        cfg = tm.LayoutConfiguration()
        cfg.k = 20
        cfg.sl_extra_scaling_steps = 10
        cfg.node_size = 1 / 25
        cfg.mmm_repeats = 2
        cfg.sl_repeats = 2
        cfg.sl_scaling_type = tm.RelativeToAvgLength

        x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
        print("Get coordinates: done")
        # To store coordinates
        x = list(x)
        y = list(y)
        s = list(s)
        t = list(t)

        pickle.dump(
            (x, y, s, t), open(set_coordinates_path, "wb+"), protocol=pickle.HIGHEST_PROTOCOL
        )

        lf.store(lf_store_path)

        descriptors_df = pd.DataFrame(descriptors)
        descriptors_df.to_csv(set_attributes_path)

        return None
    return (generate_tmap,)


@app.cell
def _(annotations_no_duplicates, generate_tmap):
    generate_tmap(df=annotations_no_duplicates, smiles_column ="smiles",
                  ik2d_column = 'ik2d', output_dir="260210_rsv_annnotations_tmap")
    return


@app.cell
def _(Path, os, pd, pickle):
    def load_tmap(output_dir: str):
        """
        Load TMAP.

        Args:
            output_dir: output directory name (saved to ../data/{output_dir})
        """
        output_path = Path("../data") / output_dir
        set_attributes_path = os.path.join(output_path, "attribute.csv")
        set_coordinates_path = os.path.join(output_path, "coordinates.dat")

        x, y, s, t = pickle.load(open(set_coordinates_path,
                                  "rb"))

        coordinates = {
            'x': x, 'y': y, 's': s, 't': t, 
        }
        descriptors = pd.read_csv(set_attributes_path)

        return coordinates, descriptors

    return (load_tmap,)


@app.cell
def _(load_tmap):
    coordinates_annot, descriptors_annot = load_tmap('260210_rsv_annnotations_tmap')
    return coordinates_annot, descriptors_annot


@app.cell
def _(descriptors_annot):
    descriptors_annot.dropna()
    return


@app.cell
def _(annotations, np, pd):
    extract_metadata = pd.read_csv('../data/anti_rsv_data.csv', sep=';')
    extract_metadata['Well'] = extract_metadata['Well'].apply(lambda x: x[0] + '0' + x[1:] if len(x) < 3 else x)
    extract_metadata['vgf_code'] = extract_metadata['Plaque'] + '_' + extract_metadata['Well']

    extract_metadata['category'] = np.where(
        ((extract_metadata.AOAO_N1 <= 20) & (extract_metadata.AOAO_N2 <= 20) &
        (extract_metadata.A549_N1 <= 20) & (extract_metadata.A549_N2 <= 20)), 'active_both', 'inactive')

    extract_metadata['category'] = np.where(
        ((extract_metadata.AOAO_N1 <= 20) & (extract_metadata.AOAO_N2 <= 20) &
        ((extract_metadata.A549_N1 > 20) | (extract_metadata.A549_N2 > 20))), 'active_organoids', extract_metadata['category'])

    extract_metadata['category'] = np.where(
        (((extract_metadata.AOAO_N1 > 20) | (extract_metadata.AOAO_N2 > 20)) &
        (extract_metadata.A549_N1 <= 20) & (extract_metadata.A549_N2 <= 20)), 'active_cells', extract_metadata['category'])

    mapping = {
        "VGF157_C10": "Litsea polyantha",
        "VGF154_H04": "Ampelocissus arachnoidea",
        "VGF157_D08": "Clausena wallichii",
        "VGF154_B02": "Vepris macrophylla"
    }

    extract_metadata["selection"] = extract_metadata["vgf_code"].map(mapping).fillna("Other")

    annotations_with_metadata = annotations.merge(extract_metadata, left_on='sample', right_on='vgf_code', how='left')
    annotations_with_metadata = annotations_with_metadata.groupby('ik2d').agg({
        'Family': lambda x: x.unique().tolist(),
        'category': lambda x: x.unique().tolist(),
        'npc_class': lambda x: x.unique().tolist(),
        'npc_superclass': lambda x: x.unique().tolist(),
        'npc_pathway': lambda x: x.unique().tolist(),
        'vgf_code': lambda x: x.unique().tolist(),
        'selection': lambda x: x.unique().tolist()}).reset_index()

    species_list = [
        "Litsea polyantha",
        "Ampelocissus arachnoidea",
        "Clausena wallichii",
        "Vepris macrophylla"
    ]

    for species in species_list:
        annotations_with_metadata[species] = (
            annotations_with_metadata["selection"]
            .apply(lambda x: "Yes" if species in x else "No")
        )

    annotations_with_metadata
    return annotations_with_metadata, mapping


@app.cell
def _(
    Faerun,
    Path,
    Top_N_classes,
    annotations_with_metadata,
    coordinates_annot,
    descriptors_annot,
    mcolors,
    np,
    os,
    ss,
):
    output_dir = '260210_rsv_annnotations_tmap'
    output_path = Path("../data") / output_dir
    tmap_filename = set_coordinates_path = os.path.join(output_path, "tmap")

    descriptors_annot["labels"] = (
            descriptors_annot["smiles"]
            + '__'
            + descriptors_annot["ik2d"]
            + "</a>"
    )    

    descriptors_annot["c_frac_ranked"] = ss.rankdata(
        np.array(descriptors_annot['c_frac']) / max(descriptors_annot['c_frac'])) / len(descriptors_annot['c_frac'])

    descriptors_annot_merged = descriptors_annot.merge(annotations_with_metadata[
        ['ik2d', "Litsea polyantha", "Ampelocissus arachnoidea", "Clausena wallichii", "Vepris macrophylla"]
        ], on= 'ik2d', how='left').fillna('No')

    dic_categories = {
        "Litsea polyantha": {'Ncat': 0},
        "Ampelocissus arachnoidea": {'Ncat': 0},
        "Clausena wallichii": {'Ncat': 0},
        "Vepris macrophylla": {'Ncat': 0}
    }

    for dic in dic_categories:
        labels, data = Faerun.create_categories(descriptors_annot_merged[str(dic)])
        dic_categories[dic]['data'], dic_categories[dic]['labels'] = Top_N_classes(dic_categories[dic]['Ncat'], data, labels)

    cmap = mcolors.ListedColormap(['#e6e6e6', "#023047"])
    cmap_litsea = mcolors.ListedColormap(['#e6e6e6', "#4394f7"])
    cmap_ampelocissus = mcolors.ListedColormap(['#e6e6e6', "#f0861d"])
    cmap_clausena = mcolors.ListedColormap(['#e6e6e6', "#d316f9"])
    cmap_vepris = mcolors.ListedColormap(['#e6e6e6', "#575763"])

    f = Faerun(view="front", coords=False,  clear_color='#ffffff')
    f.add_scatter(
        "annotations",
        {
            "x": coordinates_annot['x'],
            "y": coordinates_annot['y'],
            "c": [
                dic_categories['Litsea polyantha']['data'],
                dic_categories['Ampelocissus arachnoidea']['data'],
                dic_categories['Clausena wallichii']['data'],
                dic_categories['Vepris macrophylla']['data'],
                descriptors_annot_merged['hac'],
                descriptors_annot_merged['c_frac_ranked'],
                descriptors_annot_merged['ring_atom_frac'],
                descriptors_annot_merged['largest_ring_size'],
            ],
            "labels": descriptors_annot_merged["labels"],
        },
        shader="smoothCircle",
        point_scale=7.0,
        max_point_size=10,
        legend_labels=[
            dic_categories['Litsea polyantha']['labels'],
            dic_categories['Ampelocissus arachnoidea']['labels'],
            dic_categories['Clausena wallichii']['labels'],
            dic_categories['Vepris macrophylla']['labels']
            ],
        categorical=[True, True, True, True, False, False, False, False],
        colormap=[cmap_litsea, cmap_ampelocissus, cmap_clausena, cmap_vepris, "rainbow", "rainbow", "rainbow", "rainbow"],
        series_title=["Annotated in Litsea polyantha", "Annotated in Ampelocissus arachnoidea",
                      "Annotated in Clausena wallichii", "Annotated in Vepris macrophylla",
                      "HAC", "C Frac", "Ring Atom Frac", "Largest Ring Size"],
        has_legend=True
    )
    f.add_tree("annotations_tree", {"from": coordinates_annot['s'], "to": coordinates_annot['t']},
              point_helper="annotations", color='#e6e6e6')
    f.plot(output_dir, template="smiles")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # TMAP Lotus

    Get TMAP from Lotus compounds
    """)
    return


@app.cell(hide_code=True)
def _(AllChem, MAP4Calculator, Path, os, pd, pickle, tm, tqdm):
    def generate_lotus_tmap(
            df: pd.DataFrame,
            smiles_column: str,
            ik2d_column: str,
            output_dir: str,
            map4_dimensions: int = 1024,
            lsh_forest_depth: int = 64,
        ) -> None:
        """
        Generate TMAP visualization from chemical structures.

        Args:
            df: Input dataframe with chemical structures
            smiles_column: Name of column containing SMILES strings
            ik2d_column: Name of column containing InChIKey 2D
            output_dir: Output directory name (saved to ../data/{output_dir})
            map4_dimensions: MAP4 fingerprint dimensions
            lsh_forest_depth: LSH Forest depth parameter
        """

        output_path = Path("../data") / output_dir
        output_path.mkdir(parents=True, exist_ok=True)
        lf_store_path = os.path.join(output_path, "lsh_forest.dat")
        set_attributes_path = os.path.join(output_path, "attribute.csv")
        set_coordinates_path = os.path.join(output_path, "coordinates.dat")

        # Calculate fingerprints and descriptors
        fps = []
        descriptors = {"ik2d":[], "smiles":[], "hac": [], "c_frac": [], "ring_atom_frac": [], "largest_ring_size": []}

        for _, row in tqdm(df.iterrows()):
            mol = AllChem.MolFromSmiles(row[smiles_column])
            if mol is None:
                continue

            # Calculate descriptors
            atoms = mol.GetAtoms()
            size = mol.GetNumHeavyAtoms()
            n_c = sum(1 for atom in atoms if atom.GetSymbol().lower() == "c")
            n_ring_atoms = sum(1 for atom in atoms if atom.IsInRing())

            descriptors["ik2d"].append(row[ik2d_column])
            descriptors["smiles"].append(row[smiles_column])
            descriptors["hac"].append(size)
            descriptors["c_frac"].append(n_c / size if size > 0 else 0)
            descriptors["ring_atom_frac"].append(n_ring_atoms / size if size > 0 else 0)

            sssr = AllChem.GetSymmSSSR(mol)
            descriptors["largest_ring_size"].append(max([len(s) for s in sssr]) if sssr else 0)
            fps.append(mol)

        print("Parse smiles and get descriptors: done")
        # Build LSH Forest and TMAP
        map4 = MAP4Calculator(dimensions=map4_dimensions)
        lf = tm.LSHForest(map4_dimensions, lsh_forest_depth)

        fps = map4.calculate_many(fps)
        print("Get fingerprints: done")
        lf.batch_add(fps)
        lf.index()
        print("Get LSH forest: done")

        cfg = tm.LayoutConfiguration()
        cfg.k = 20
        cfg.sl_extra_scaling_steps = 10
        cfg.node_size = 1 / 25
        cfg.mmm_repeats = 2
        cfg.sl_repeats = 2
        cfg.sl_scaling_type = tm.RelativeToAvgLength

        x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
        print("Get coordinates: done")
        # To store coordinates
        x = list(x)
        y = list(y)
        s = list(s)
        t = list(t)

        pickle.dump(
            (x, y, s, t), open(set_coordinates_path, "wb+"), protocol=pickle.HIGHEST_PROTOCOL
        )

        lf.store(lf_store_path)

        descriptors_df = pd.DataFrame(descriptors)
        descriptors_df.to_csv(set_attributes_path)

        return None
    return (generate_lotus_tmap,)


@app.cell
def _(pd):
    #annotations = get_annotations()
    lotus = pd.read_csv('../data/20260206_lotus_all_taxa.csv.gz', compression='gzip')
    lotus['ik2d'] = lotus['compound_inchikey'].str[:14]

    print(f"""
    Number of unique 3D structures: {lotus.compound_inchikey.nunique()}
    Number of unique 2D structures: {lotus.ik2d.nunique()}
    Number of unique taxa: {lotus.taxon_name.nunique()}
    """)

    lotus.head()
    return (lotus,)


@app.cell
def _(AllChem, Chem, lotus, pd, tqdm):
    lotus_2d_for_tmap = {}
    lotus_2d_for_tmap["ik2d"] = []
    lotus_2d_for_tmap["smiles"] = []

    lotus_2d = lotus.drop_duplicates(subset='ik2d')

    for _, row in tqdm(lotus_2d.iterrows()):
        mol = AllChem.MolFromSmiles(row["compound_smiles"])
        if mol is None:
            continue
        lotus_2d_for_tmap["ik2d"].append(row["ik2d"])
        lotus_2d_for_tmap["smiles"].append(Chem.MolToSmiles(mol, isomericSmiles=False))

    lotus_2d_for_tmap = pd.DataFrame(lotus_2d_for_tmap)
    lotus_2d_for_tmap
    return (lotus_2d_for_tmap,)


@app.cell
def _(generate_lotus_tmap, lotus_2d_for_tmap):
    generate_lotus_tmap(df=lotus_2d_for_tmap, smiles_column ="smiles",
                  ik2d_column = 'ik2d', output_dir="260208_lotus_tmap")
    return


@app.cell
def _(load_tmap):
    coordinates_annot_lotus, descriptors_annot_lotus = load_tmap('260208_lotus_tmap')
    return coordinates_annot_lotus, descriptors_annot_lotus


@app.cell
def _(annotations_no_duplicates):
    annotations_no_duplicates
    return


@app.cell
def _(annotations_no_duplicates, descriptors_annot_lotus, np):
    descriptors_annot_lotus['annotated'] = np.where(descriptors_annot_lotus.ik2d.isin(annotations_no_duplicates.ik2d), 'Yes', 'No')
    descriptors_annot_lotus
    return


@app.cell
def _(
    Faerun,
    Top_N_classes,
    coordinates_annot_lotus,
    descriptors_annot_lotus,
    mcolors,
    np,
    ss,
):
    output_dir_lotus = '260223_lotus_tmap'

    descriptors_annot_lotus["labels"] = (
            descriptors_annot_lotus["smiles"]
            + '__'
            + descriptors_annot_lotus["ik2d"]
            + "</a>"

    )    

    descriptors_annot_lotus["c_frac_ranked"] = ss.rankdata(
        np.array(descriptors_annot_lotus['c_frac']) / max(descriptors_annot_lotus['c_frac'])) / len(descriptors_annot_lotus['c_frac'])

    dic_categories2 = {
        "annotated": {'Ncat': 0}
    }

    for dic2 in dic_categories2:
        labels2, data2 = Faerun.create_categories(descriptors_annot_lotus[str(dic2)])
        dic_categories2[dic2]['data'], dic_categories2[dic2]['labels'] = Top_N_classes(dic_categories2[dic2]['Ncat'], data2, labels2)

    cmap2 = mcolors.ListedColormap(['#e6e6e6', "#023047"])

    f2 = Faerun(view="front", coords=False,  clear_color='#ffffff')
    f2.add_scatter(
        "annotations",
        {
            "x": coordinates_annot_lotus['x'],
            "y": coordinates_annot_lotus['y'],
            "c": [
                dic_categories2['annotated']['data'],
                descriptors_annot_lotus['hac'].astype(float).to_list(),
                descriptors_annot_lotus['c_frac_ranked'].astype(float).to_list(),
                descriptors_annot_lotus['ring_atom_frac'].astype(float).to_list(),
                descriptors_annot_lotus['largest_ring_size'].astype(float).to_list(),
            ],
            "labels": descriptors_annot_lotus["labels"],
        },
        shader="smoothCircle",
        point_scale=7.0,
        max_point_size=5,
        legend_labels=[
            dic_categories2['annotated']['labels']
            ],
        categorical=[True, False, False, False, False],
        colormap=[cmap2, "rainbow", "rainbow", "rainbow", "rainbow"],
        series_title=["Annotated?", "HAC", "C Frac", "Ring Atom Frac", "Largest Ring Size"],
        has_legend=True
    )
    f2.add_tree("annotations_tree", {"from": coordinates_annot_lotus['s'], "to": coordinates_annot_lotus['t']},
              point_helper="annotations", color='#e6e6e6')
    f2.plot(output_dir_lotus, template="smiles")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Chemical Class Annotations

    Plot chemical class specificity in selected active extracts

    - "VGF157_C10": "Litsea polyantha",
    - "VGF154_H04": "Ampelocissus arachnoidea",
    - "VGF157_D08": "Clausena wallichii",
    - "VGF154_B02": "Vepris macrophylla"
    """)
    return


@app.cell
def _():
    sample_to_plot = 'VGF157_D08'
    threshold = 5 #minimal frequency of the chemical class in the sample_to_plot to be retained
    return sample_to_plot, threshold


@app.cell
def _(pd):
    canopus_annotation = pd.read_csv('../data/canopus_rsv_subset.csv', sep=';')
    canopus_annotation['np_class'] = canopus_annotation['np_class'].str.removeprefix('npc_')
    canopus_annotation['np_class'] = canopus_annotation['np_class'].replace({"_":" "}, regex=True)
    return (canopus_annotation,)


@app.cell
def _(canopus_annotation):
    canopus_annotation
    return


@app.cell
def _(canopus_annotation, pd, sample_to_plot, threshold):
    grouped_sample_class = canopus_annotation.groupby(['sample', 'np_class'] ).count().reset_index()
    sample_df = grouped_sample_class[grouped_sample_class['sample'] ==sample_to_plot]

    non_sample_df = grouped_sample_class[grouped_sample_class['sample'] !=sample_to_plot]
    all_samples = non_sample_df['sample'].unique()
    all_classes = non_sample_df['np_class'].unique()

    # 2. Build a MultiIndex from the product of these two lists
    full_index = pd.MultiIndex.from_product(
        [all_samples, all_classes], 
        names=['sample', 'np_class']
    )

    # 3. Set the current index to match the MultiIndex, then reindex
    # Any missing pairs will be created with NaN values
    df_filled = (non_sample_df.set_index(['sample', 'np_class'])
                 .reindex(full_index, fill_value=0)
                 .reset_index())

    average_count_non_sample = df_filled.groupby('np_class')['ion'].mean()

    plotting_df = sample_df.merge(average_count_non_sample, on='np_class', how='left')
    plotting_df['specificity_score'] = plotting_df['ion_x']/plotting_df['ion_y']
    plotting_df = plotting_df[plotting_df['ion_x']>=threshold]

    plotting_df
    return (plotting_df,)


@app.cell
def _(alt, mapping, plotting_df, sample_to_plot):
    # Shared axis titles
    x_axis_title = f"Count in {mapping[sample_to_plot]}"
    y_axis_title = "Average count in other extracts"

    # 1. The background: Gray dots for specificity_score <= 5
    # We explicitly remove the legend for this layer
    background = alt.Chart(plotting_df).mark_circle(color='black').encode(
        x=alt.X('ion_x', scale=alt.Scale(type='log', domain=[5, plotting_df['ion_x'].max()]), title=x_axis_title),
        y=alt.Y('ion_y', scale=alt.Scale(type='log'), title=y_axis_title),
        size=alt.Size('specificity_score', legend=None)
    ).transform_filter(
        alt.datum.specificity_score <= 10
    )

    # 2. The foreground: Colored dots for specificity_score > 5
    # This layer will generate the legend automatically
    foreground = alt.Chart(plotting_df).mark_circle().encode(
        x=alt.X('ion_x', scale=alt.Scale(type='log', domain=[5, plotting_df['ion_x'].max()])),
        y=alt.Y('ion_y', scale=alt.Scale(type='log')),
        size=alt.Size('specificity_score'),
        color=alt.Color('np_class:N', title="NP Class (Score > 10)")
    ).transform_filter(
        alt.datum.specificity_score > 10
    )

    # Combine them (Background first so colored dots sit on top)
    chart = (background + foreground).properties(
        width=300,
        height=200
    ).configure_title(
        # This ensures any markdown in titles is rendered correctly
        fontStyle='italic' 
    )

    chart
    return (chart,)


@app.cell
def _(chart, mapping, sample_to_plot):
    chart.save(f'{mapping[sample_to_plot]} specificity.png', ppi=600)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Upset plots

    Uncover the number of annotations that are specific and/or shared across active/inactive extracts
    """)
    return


@app.cell
def _():
    from upsetplot import UpSet
    from upsetplot import from_memberships
    from upsetplot import generate_counts, plot
    return from_memberships, plot


@app.cell
def _(annotations_with_metadata):
    annotations_with_metadata
    return


@app.cell
def _(annotations_with_metadata, from_memberships):
    categories = from_memberships(annotations_with_metadata.category)
    categories = categories.rename_axis(index={
        'active_both': 'Active both models',
        'active_cells': 'Active A549',
        'active_organoids': 'Active AOAO',
        'inactive': 'Inactive both models',
    })
    return (categories,)


@app.cell
def _(categories, plot, plt):
    aggregated = categories.groupby(level=list(range(categories.index.nlevels))).size()
    plot(aggregated, show_counts=True, sort_by='cardinality')
    plt.savefig('upsetplot.png', dpi=600)
    plt.show()
    return


if __name__ == "__main__":
    app.run()
