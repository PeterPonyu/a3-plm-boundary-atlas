"""Fail-closed A3 source binding before the GPU/weight-dependent scoring step."""
from pathlib import Path
import hashlib
import importlib.util

import pytest
import torch

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("a3_candidate_rescore", HERE/"rescore_clinvar_candidate.py")
rescore = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rescore)


def test_manifest_and_source_are_one_to_one_and_reference_sequences_unique():
    grouped, refs, hashes, n6 = rescore.source_rows()
    assert set(grouped) == set(refs) == set(hashes) == set(n6["per_gene"])
    assert len(grouped) == 20 and sum(map(len, grouped.values())) == 4677
    assert all(rescore.source_wt(gene, refs[gene]) for gene in grouped)
    assert all(len(rescore.source_wt(gene, refs[gene])) == n6["per_gene"][gene]["wt_len"]
               for gene in grouped)


def test_changed_source_reference_is_rejected():
    grouped, refs, _, _ = rescore.source_rows()
    gene = next(iter(grouped))
    with pytest.raises(ValueError, match="Source sequence changed"):
        rescore.source_wt(gene, hashlib.sha256(b"incorrect reference").hexdigest())


def test_incomplete_checkpoint_is_rejected_before_model_import(tmp_path):
    checkpoint = tmp_path / "partial.pt"
    checkpoint.write_bytes(b"not a model")
    with pytest.raises(ValueError, match="Incomplete checkpoint"):
        rescore.require_complete_checkpoint(checkpoint)


def test_residue_index_and_windowed_llr_use_same_reference():
    class FakeAlphabet:
        def get_idx(self, aa):
            return 2 if aa == "C" else 1

    class FakeModel:
        def __call__(self, tokens):
            logits = torch.zeros((1, tokens.shape[1], 3), dtype=torch.float32)
            logits[:, :, 2] = 2.0
            return {"logits": logits}

    def convert(pairs):
        return ["wt"], [pairs[0][1]], torch.zeros((1, len(pairs[0][1])+2), dtype=torch.long)

    entries = [dict(example_id="toy1", ligand_id="A1C", label_value="1"),
               dict(example_id="toy2", ligand_id="A1025C", label_value="0")]
    rows = rescore.score_gene("TOY", entries, "A" * 1030, FakeModel(), FakeAlphabet(),
                             convert, "cpu")
    assert len(rows) == 2
    assert all(r["primary_mode"] == "windowed" for r in rows)
    assert all(r["primary_llr_fp32"] == pytest.approx(2.0) for r in rows)
    assert all(r["llr_whole_fp32"] == pytest.approx(r["llr_window_fp32"]) for r in rows)
    assert all(r["window_end_exclusive"]-r["window_start_zero_based"] <= rescore.CTX
               for r in rows)
    assert rows[0]["window_start_zero_based"] == 0
    assert rows[1]["window_start_zero_based"] > 0
