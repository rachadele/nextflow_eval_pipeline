"""
Tests for the refactored Census reference-fetching functions in utils.py:
  subsample_cells, extract_ref, split_refs, get_census
"""
import sys
import os
import pytest
import pandas as pd
import numpy as np
import anndata as ad
from unittest.mock import patch, MagicMock, call

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "bin"))
from utils import subsample_cells, extract_ref, split_refs, get_census, CELL_COLUMNS

REF_KEYS = ["subclass", "class", "family"]


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def relabel_path(tmp_path):
    """Minimal relabel TSV (cell_type is the join key, matching census_map_human.tsv format)."""
    content = (
        "cell_type\tsubclass\tclass\tfamily\n"
        "typeA\tSubA\tClassA\tFamilyA\n"
        "typeB\tSubB\tClassB\tFamilyB\n"
        "typeC\tSubA\tClassA\tFamilyA\n"
    )
    p = tmp_path / "relabel.tsv"
    p.write_text(content)
    return str(p)


def make_obs(n_per_type, cell_types=("typeA", "typeB"), extra_cols=None):
    """Build a minimal obs DataFrame suitable for subsample_cells."""
    rows = []
    joinid = 0
    for ct in cell_types:
        for _ in range(n_per_type):
            row = {"cell_type": ct, "soma_joinid": joinid}
            if extra_cols:
                row.update(extra_cols)
            rows.append(row)
            joinid += 1
    return pd.DataFrame(rows)


def make_split_obs(dataset_ids, n_per_dataset=1200, collection_name="Good Collection"):
    """Build obs with dataset_id and collection_name columns for split_refs tests."""
    rows = []
    joinid = 0
    for did in dataset_ids:
        for i in range(n_per_dataset):
            rows.append({
                "cell_type": "typeA" if i % 2 == 0 else "typeB",
                "soma_joinid": joinid,
                "dataset_id": did,
                "collection_name": collection_name,
                "dataset_title": f"Title {did}",
            })
            joinid += 1
    return pd.DataFrame(rows)


def make_mock_adata(obs_df):
    """Return a minimal AnnData whose obs matches obs_df."""
    n = len(obs_df)
    adata = ad.AnnData(X=np.zeros((n, 5)))
    adata.obs = obs_df.reset_index(drop=True)
    adata.obs.index = adata.obs.index.astype(str)
    return adata


# ---------------------------------------------------------------------------
# subsample_cells
# ---------------------------------------------------------------------------

class TestSubsampleCells:

    def test_caps_at_subsample_limit(self, relabel_path):
        obs = make_obs(n_per_type=600)
        result = subsample_cells(obs, subsample=500, relabel_path=relabel_path,
                                 ref_keys=REF_KEYS, seed=42)
        # Two subclasses (SubA, SubB), each capped at 500
        assert len(result) == 1000

    def test_takes_all_when_below_limit(self, relabel_path):
        obs = make_obs(n_per_type=10)
        result = subsample_cells(obs, subsample=500, relabel_path=relabel_path,
                                 ref_keys=REF_KEYS, seed=42)
        assert len(result) == 20

    def test_deterministic_with_seed(self, relabel_path):
        obs = make_obs(n_per_type=600)
        r1 = subsample_cells(obs, subsample=100, relabel_path=relabel_path,
                             ref_keys=REF_KEYS, seed=42)
        r2 = subsample_cells(obs, subsample=100, relabel_path=relabel_path,
                             ref_keys=REF_KEYS, seed=42)
        assert r1 == r2

    def test_different_seeds_differ(self, relabel_path):
        obs = make_obs(n_per_type=600)
        r1 = subsample_cells(obs, subsample=100, relabel_path=relabel_path,
                             ref_keys=REF_KEYS, seed=0)
        r2 = subsample_cells(obs, subsample=100, relabel_path=relabel_path,
                             ref_keys=REF_KEYS, seed=99)
        assert r1 != r2

    def test_raises_when_join_key_missing(self, relabel_path):
        obs = pd.DataFrame({"soma_joinid": [1, 2], "wrong_col": ["x", "y"]})
        with pytest.raises(ValueError, match="cell_type not found in obs"):
            subsample_cells(obs, subsample=10, relabel_path=relabel_path,
                            ref_keys=REF_KEYS, seed=42)

    def test_subsamples_per_subclass_not_cell_type(self, relabel_path):
        # typeA and typeC both map to SubA — combined they form one bucket
        obs = make_obs(n_per_type=400, cell_types=("typeA", "typeC"))
        result = subsample_cells(obs, subsample=500, relabel_path=relabel_path,
                                 ref_keys=REF_KEYS, seed=42)
        # SubA bucket has 800 cells but we only subsample by ref_keys[0] (subclass),
        # so both typeA and typeC contribute to the same SubA bucket
        assert len(result) == 500

    def test_returns_soma_joinids(self, relabel_path):
        obs = make_obs(n_per_type=10)
        result = subsample_cells(obs, subsample=500, relabel_path=relabel_path,
                                 ref_keys=REF_KEYS, seed=42)
        assert all(jid in obs["soma_joinid"].values for jid in result)


# ---------------------------------------------------------------------------
# extract_ref
# ---------------------------------------------------------------------------

class TestExtractRef:

    @pytest.fixture
    def dataset_info(self):
        return pd.DataFrame({
            "dataset_id": ["ds1"],
            "dataset_title": ["My Dataset"],
            "collection_name": ["My Collection"],
        })

    @pytest.fixture
    def obs(self, relabel_path):
        return make_obs(n_per_type=10, extra_cols={"dataset_id": "ds1"})

    def _mock_adata(self, obs_df):
        """AnnData with dataset_id in obs, returned by mocked get_anndata."""
        return make_mock_adata(obs_df[["cell_type", "dataset_id"]].copy())

    def test_calls_get_anndata_with_subsampled_ids(self, obs, dataset_info, relabel_path):
        mock_adata = self._mock_adata(obs)
        with patch("utils.cellxgene_census.get_anndata", return_value=mock_adata) as mock_ga, \
             patch("utils.sc.pp.filter_genes"), \
             patch("utils.relabel", side_effect=lambda a, **kw: a):

            extract_ref(obs, census=MagicMock(), organism="Homo sapiens",
                        dataset_info=dataset_info, subsample=5,
                        relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)

            args, kwargs = mock_ga.call_args
            assert "obs_coords" in kwargs
            assert len(kwargs["obs_coords"]) <= 10   # 5 per type, 2 types
            assert kwargs["obs_embeddings"] == ["scvi"]
            assert set(kwargs["obs_column_names"]) == set(CELL_COLUMNS)

    def test_merges_dataset_info_onto_obs(self, obs, dataset_info, relabel_path):
        mock_adata = self._mock_adata(obs)
        with patch("utils.cellxgene_census.get_anndata", return_value=mock_adata), \
             patch("utils.sc.pp.filter_genes"), \
             patch("utils.relabel", side_effect=lambda a, **kw: a):

            result = extract_ref(obs, census=MagicMock(), organism="Homo sapiens",
                                 dataset_info=dataset_info, subsample=5,
                                 relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)

            assert "dataset_title" in result.obs.columns
            assert "collection_name" in result.obs.columns

    def test_maps_author_labels_when_provided(self, obs, dataset_info, relabel_path):
        from utils import get_original_celltypes
        census_version = "2025-01-30"
        original_celltypes = get_original_celltypes(
            columns_file=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}/original_celltype_columns.tsv",
            author_annotations_path=f"/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/{census_version}",
        )
        mock_adata = self._mock_adata(obs)
        with patch("utils.cellxgene_census.get_anndata", return_value=mock_adata), \
             patch("utils.sc.pp.filter_genes"), \
             patch("utils.relabel", side_effect=lambda a, **kw: a), \
             patch("utils.map_author_labels", return_value=mock_adata.obs) as mock_mal:

            extract_ref(obs, census=MagicMock(), organism="Homo sapiens",
                        dataset_info=dataset_info, subsample=5,
                        relabel_path=relabel_path, ref_keys=REF_KEYS,
                        original_celltypes=original_celltypes, seed=42)

            mock_mal.assert_called_once()

    def test_skips_author_labels_when_not_provided(self, obs, dataset_info, relabel_path):
        mock_adata = self._mock_adata(obs)
        with patch("utils.cellxgene_census.get_anndata", return_value=mock_adata), \
             patch("utils.sc.pp.filter_genes"), \
             patch("utils.relabel", side_effect=lambda a, **kw: a), \
             patch("utils.map_author_labels") as mock_mal:

            extract_ref(obs, census=MagicMock(), organism="Homo sapiens",
                        dataset_info=dataset_info, subsample=5,
                        relabel_path=relabel_path, ref_keys=REF_KEYS,
                        original_celltypes=None, seed=42)

            mock_mal.assert_not_called()


# ---------------------------------------------------------------------------
# split_refs
# ---------------------------------------------------------------------------

class TestSplitRefs:

    @pytest.fixture
    def dataset_info(self):
        return pd.DataFrame({
            "dataset_id": ["ds1", "ds2", "ds3"],
            "dataset_title": ["Title ds1", "Title ds2", "Title ds3"],
            "collection_name": ["Good", "Good", "Good"],
        })

    def _make_extract_ref_result(self, dataset_title):
        obs = pd.DataFrame({"dataset_title": [dataset_title] * 5,
                            "cell_type": ["typeA"] * 5})
        return make_mock_adata(obs)

    def test_skips_hvs_collection(self, dataset_info, relabel_path):
        obs = make_split_obs(["ds1"], n_per_dataset=1200,
                             collection_name="HVS: Human variation study")
        with patch("utils.extract_ref") as mock_er:
            result = split_refs(obs, split_column="dataset_id",
                                census=MagicMock(), organism="Homo sapiens",
                                dataset_info=dataset_info, subsample=10,
                                relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)
        mock_er.assert_not_called()
        assert result == {}

    def test_skips_splits_with_fewer_than_1000_cells(self, dataset_info, relabel_path):
        obs = make_split_obs(["ds1"], n_per_dataset=500)
        with patch("utils.extract_ref") as mock_er:
            result = split_refs(obs, split_column="dataset_id",
                                census=MagicMock(), organism="Homo sapiens",
                                dataset_info=dataset_info, subsample=10,
                                relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)
        mock_er.assert_not_called()
        assert result == {}

    def test_keys_by_dataset_title_when_split_by_dataset_id(self, dataset_info, relabel_path):
        obs = make_split_obs(["ds1", "ds2"], n_per_dataset=1200)

        def fake_extract_ref(split_obs, **kwargs):
            did = split_obs["dataset_id"].iloc[0]
            return make_mock_adata(pd.DataFrame({"dataset_title": [f"Title {did}"] * 3}))

        with patch("utils.extract_ref", side_effect=fake_extract_ref):
            result = split_refs(obs, split_column="dataset_id",
                                census=MagicMock(), organism="Homo sapiens",
                                dataset_info=dataset_info, subsample=10,
                                relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)

        assert "Title ds1" in result
        assert "Title ds2" in result

    def test_keys_by_split_value_when_split_by_tissue(self, dataset_info, relabel_path):
        obs = make_split_obs(["ds1"], n_per_dataset=1200)
        obs["tissue"] = "cortex"

        def fake_extract_ref(split_obs, **kwargs):
            return make_mock_adata(pd.DataFrame({"dataset_title": ["Title ds1"] * 3}))

        with patch("utils.extract_ref", side_effect=fake_extract_ref):
            result = split_refs(obs, split_column="tissue",
                                census=MagicMock(), organism="Homo sapiens",
                                dataset_info=dataset_info, subsample=10,
                                relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)

        assert "cortex" in result

    def test_extract_ref_receives_sliced_obs(self, dataset_info, relabel_path):
        """extract_ref should receive only the cells for its split, not the full obs."""
        obs = make_split_obs(["ds1", "ds2"], n_per_dataset=1200)
        captured = {}

        def fake_extract_ref(split_obs, **kwargs):
            did = split_obs["dataset_id"].iloc[0]
            captured[did] = len(split_obs)
            return make_mock_adata(pd.DataFrame({"dataset_title": [f"Title {did}"] * 3}))

        with patch("utils.extract_ref", side_effect=fake_extract_ref):
            split_refs(obs, split_column="dataset_id",
                       census=MagicMock(), organism="Homo sapiens",
                       dataset_info=dataset_info, subsample=10,
                       relabel_path=relabel_path, ref_keys=REF_KEYS, seed=42)

        assert captured["ds1"] == 1200
        assert captured["ds2"] == 1200


# ---------------------------------------------------------------------------
# get_census
# ---------------------------------------------------------------------------

class TestGetCensus:

    def test_uses_context_manager_and_returns_whole_cortex(self, relabel_path):
        mock_census = MagicMock()
        mock_census.__enter__ = MagicMock(return_value=mock_census)
        mock_census.__exit__ = MagicMock(return_value=False)

        dataset_info = pd.DataFrame({
            "dataset_id": ["ds1"],
            "dataset_title": ["Title ds1"],
            "collection_name": ["Good Collection"],
            "soma_joinid": [99],
        })
        mock_census["census_info"]["datasets"].read().concat().to_pandas.return_value = dataset_info

        filtered_obs = pd.DataFrame({
            "cell_type": ["typeA"] * 5,
            "soma_joinid": list(range(5)),
            "dataset_id": ["ds1"] * 5,
        })

        with patch("utils.cellxgene_census.open_soma", return_value=mock_census), \
             patch("utils.get_filtered_obs", return_value=filtered_obs), \
             patch("utils.split_refs", return_value={"ds1 ref": MagicMock()}) as mock_sr, \
             patch("utils.extract_ref", return_value=MagicMock()) as mock_er, \
             patch("utils.map_author_labels", side_effect=lambda o, _: o):

            result = get_census(
                census_version="2024-07-01",
                organism="homo_sapiens",
                subsample=5,
                ref_collections=["Good Collection"],
                relabel_path=relabel_path,
                ref_keys=REF_KEYS,
                seed=42,
            )

        assert "whole cortex" in result
        mock_sr.assert_called_once()
        mock_er.assert_called_once()

    def test_passes_remapped_organism_to_split_refs_and_extract_ref(self, relabel_path):
        mock_census = MagicMock()
        mock_census.__enter__ = MagicMock(return_value=mock_census)
        mock_census.__exit__ = MagicMock(return_value=False)

        dataset_info = pd.DataFrame({
            "dataset_id": ["ds1"],
            "dataset_title": ["Title ds1"],
            "collection_name": ["Good Collection"],
            "soma_joinid": [99],
        })
        mock_census["census_info"]["datasets"].read().concat().to_pandas.return_value = dataset_info

        filtered_obs = pd.DataFrame({
            "cell_type": ["typeA"] * 5,
            "soma_joinid": list(range(5)),
            "dataset_id": ["ds1"] * 5,
        })

        captured_organism = {}

        def fake_split_refs(obs, **kwargs):
            captured_organism["split"] = kwargs["organism"]
            return {}

        def fake_extract_ref(obs, **kwargs):
            captured_organism["extract"] = kwargs["organism"]
            return MagicMock()

        with patch("utils.cellxgene_census.open_soma", return_value=mock_census), \
             patch("utils.get_filtered_obs", return_value=filtered_obs), \
             patch("utils.split_refs", side_effect=fake_split_refs), \
             patch("utils.extract_ref", side_effect=fake_extract_ref):

            get_census(
                census_version="2024-07-01",
                organism="mus_musculus",
                ref_collections=["Good Collection"],
                relabel_path=relabel_path,
                ref_keys=REF_KEYS,
            )

        assert captured_organism["split"] == "Mus musculus"
        assert captured_organism["extract"] == "Mus musculus"
