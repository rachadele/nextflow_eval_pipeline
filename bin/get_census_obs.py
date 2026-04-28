#!/usr/bin/env python3

import cellxgene_census
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract CellXGene Census obs for a given collection.")
    parser.add_argument('--organism', type=str, default='homo_sapiens')
    parser.add_argument('--organ', type=str, default='brain')
    parser.add_argument('--census_version', type=str, default='2025-01-30')
    parser.add_argument('--ref_collections', type=str, nargs='+', default=[
        "Molecular and cellular evolution of the primate dorsolateral prefrontal cortex"
    ])
    parser.add_argument('--outfile', type=str, default='census_obs.tsv')

    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def main():
    args = parse_arguments()

    census = cellxgene_census.open_soma(census_version=args.census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()

    value_filter = (
        f"tissue_general == '{args.organ}' and "
        f"is_primary_data == True and "
        f"disease == 'normal'"
    )
    cellxgene_obs = cellxgene_census.get_obs(census, args.organism, value_filter=value_filter)

    cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None, "_y"))
    cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
    cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(args.ref_collections)]

    cellxgene_obs_filtered.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote {len(cellxgene_obs_filtered)} rows to {args.outfile}")

if __name__ == "__main__":
    main()
