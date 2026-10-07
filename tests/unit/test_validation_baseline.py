from indep_mu.workflow.baseline import BaselineMismatch, _compare_subset


def _mismatches(expected, actual, *, atol=1e-6):
    result: list[BaselineMismatch] = []
    _compare_subset(
        expected,
        actual,
        path="",
        mismatches=result,
        float_atol=atol,
    )
    return result


def test_expected_mapping_is_a_subset_of_observation() -> None:
    expected = {
        "schema_version": 1,
        "case_id": "case001",
        "status": "METADATA_VERIFIED",
        "ct": {
            "slice_count": 187,
            "pixel_spacing_mm": [1.5625, 1.5625],
        },
    }
    actual = {
        "schema_version": 1,
        "ct": {
            "slice_count": 187,
            "pixel_spacing_mm": [1.5625, 1.5625],
            "extra_runtime_value": 123,
        },
        "extra_section": {"ignored": True},
    }

    assert _mismatches(expected, actual) == []


def test_float_comparison_uses_absolute_tolerance() -> None:
    assert _mismatches(
        {"mu": 694.867399},
        {"mu": 694.8673994},
        atol=1e-6,
    ) == []

    mismatches = _mismatches(
        {"mu": 694.867399},
        {"mu": 694.8685},
        atol=1e-6,
    )
    assert len(mismatches) == 1
    assert mismatches[0].path == "mu"


def test_list_order_and_length_are_regression_significant() -> None:
    expected = {
        "beams": [
            {"beam_number": 1, "mu": 100.0},
            {"beam_number": 2, "mu": 200.0},
        ]
    }

    mismatches = _mismatches(
        expected,
        {
            "beams": [
                {"beam_number": 2, "mu": 200.0},
                {"beam_number": 1, "mu": 100.0},
            ]
        },
    )
    assert mismatches


def test_missing_expected_key_is_reported_with_path() -> None:
    mismatches = _mismatches(
        {"rtplan": {"beam_count": 8}},
        {"rtplan": {}},
    )

    assert mismatches == [
        BaselineMismatch(
            "rtplan.beam_count",
            8,
            "<missing>",
        )
    ]
