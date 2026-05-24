"""
Retrieve annotations for samples listed in `input_csv` and save to `output_csv`.
"""

import argparse
import os

import pandas as pd

from pandas import json_normalize
from SPARQLWrapper import SPARQLWrapper, JSON
from tqdm import tqdm

URI_REF = 'https://enpkg.commons-lab.org/kg/'

def get_annotations(input_csv=None,
                    output_csv=None):
    """Retrieve annotations for samples listed in `input_csv` and save to `output_csv`.

    Returns the resulting DataFrame.
    """
    # Get absolute paths relative to the script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    if input_csv is None:
        input_csv = os.path.join(project_root, 'data', 'anti_rsv_data.csv')
    if output_csv is None:
        output_csv = os.path.join(project_root, 'data', 'annotations_rsv_subset.csv')
    
    sparql = SPARQLWrapper('https://enpkg.commons-lab.org/graphdb/repositories/ENPKG')
    sparql.setReturnFormat(JSON)
    
    rsv_data = pd.read_csv(input_csv, sep=';')
    rsv_data['Well'] = rsv_data['Well'].apply(lambda x: x[0] + '0' + x[1:] if len(x) < 3 else x)
    rsv_data['vgf_code'] = rsv_data['Plaque'] + '_' + rsv_data['Well']

    vgf_rsv = rsv_data['vgf_code'].to_list()

    df_list = []
    for code in tqdm(vgf_rsv):
        query = f"""
        PREFIX enpkg: <https://enpkg.commons-lab.org/kg/>
        PREFIX enpkgmodule: <https://enpkg.commons-lab.org/module/>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>

        SELECT DISTINCT ?ik ?sample ?feature ?ik2d ?smiles ?npc_class ?npc_superclass ?npc_pathway 
            WHERE
            {{ 
                ?sample rdf:type enpkg:LabExtract .
                FILTER(regex(str(?sample), "{code}")) 
                ?sample enpkg:has_LCMS ?lcms .
                ?lcms rdf:type enpkg:LCMSAnalysis .
                ?lcms enpkg:has_lcms_feature_list ?feature_list .
                ?feature_list enpkg:has_lcms_feature ?feature .
                
                {{
                    ?feature enpkg:has_sirius_annotation ?annotation .
                    ?annotation enpkg:has_InChIkey2D ?ik2d .
                    ?annotation enpkg:has_cosmic_score ?cosmic .
                    ?annotation enpkg:has_zodiac_score ?zodiac .
                    FILTER((?cosmic > 0.5) && (?zodiac > 0.8))
                }}
                UNION
                {{
                    ?feature enpkg:has_isdb_annotation ?annotation .
                    ?annotation enpkg:has_taxo_score ?taxo .
                    FILTER(?taxo >= 6)
                }}

                ?annotation enpkg:has_InChIkey2D ?ik2d .
                ?ik2d enpkg:has_smiles ?smiles .
                ?ik2d enpkg:is_InChIkey2D_of ?ik .
                    ?ik enpkg:has_npc_class ?npc_class .
                    ?ik enpkg:has_npc_superclass ?npc_superclass .
                    ?ik enpkg:has_npc_pathway ?npc_pathway .
            }}
        """

        sparql.setQuery(query)
        results = sparql.queryAndConvert()
        df = json_normalize(results['results']["bindings"])
        if len(df) > 0:
            df = df.stack().str.replace(URI_REF, "", regex=False).unstack()
            df.drop(list(df.filter(regex = 'type')), axis = 1, inplace = True)
            df.columns = df.columns.str.replace('.value', '', regex=False)
            df_list.append(df)
        else:
            print(f'No annotation for {code}')

    if len(df_list) > 0:
        df_annotations = pd.concat(df_list, ignore_index=True)
    else:
        df_annotations = pd.DataFrame()

    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    df_annotations.to_csv(output_csv, index=False, sep=';')
    print(f'Wrote annotations to {output_csv}')
    return df_annotations


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Retrieve annotations and save CSV')
    parser.add_argument('--input', '-i', default=None,
                        help='Input CSV with sample list (default: data/anti_rsv_data.csv)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output CSV path (default: data/annotations_rsv_subset.csv)')
    args = parser.parse_args()
    get_annotations(args.input, args.output)
