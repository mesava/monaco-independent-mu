import pytest

from indep_mu.beam_model.monaco_reference import parse_monaco_leaf_model_text


def test_monaco_reference_parser_preserves_raw_values() -> None:
    model = parse_monaco_leaf_model_text(
        """
        <LeafGap>0.100000</LeafGap>
        <LeafTransmission>0.005000</LeafTransmission>
        <LeafOffset>-0.050000</LeafOffset>
        """,
        source_name="LT6.txt",
    )

    assert model.source_name == "LT6.txt"
    assert model.require("LeafGap") == pytest.approx(0.1)
    assert model.require("LeafTransmission") == pytest.approx(0.005)
    assert model.require("LeafOffset") == pytest.approx(-0.05)


def test_duplicate_reference_parameter_is_rejected() -> None:
    with pytest.raises(ValueError, match="Duplicate"):
        parse_monaco_leaf_model_text(
            "<LeafGap>0.1</LeafGap><LeafGap>0.2</LeafGap>"
        )
