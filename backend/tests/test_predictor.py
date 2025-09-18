import pytest
from backend.predictor import CYP3A4BasicPredictor

@pytest.fixture
def predictor():
    """Create predictor instance for testing"""
    return CYP3A4BasicPredictor()

def test_predict_with_nitrogen_compound(predictor):
    """Test prediction with nitrogen-containing compound"""
    result = predictor.predict("CCN")  # ethylamine

    assert "prediction" in result
    assert result["prediction"] in ["inhibitor", "non_inhibitor"]
    assert 0 <= result["probability"] <= 1
    assert 0 <= result["confidence"] <= 1
    assert "processing_time" in result

def test_predict_without_nitrogen_compound(predictor):
    """Test prediction with non-nitrogen compound"""
    result = predictor.predict("CCO")  # ethanol

    assert "prediction" in result
    assert result["prediction"] in ["inhibitor", "non_inhibitor"]
    assert 0 <= result["probability"] <= 1
    assert 0 <= result["confidence"] <= 1
    assert "warnings" in result

def test_predict_with_compound_id(predictor):
    """Test prediction with compound ID"""
    result = predictor.predict("CCO", compound_id="test_compound_1")

    assert result is not None
    assert "prediction" in result

def test_predict_invalid_smiles(predictor):
    """Test prediction with invalid SMILES"""
    with pytest.raises(ValueError):
        predictor.predict("")

def test_predictor_is_ready(predictor):
    """Test predictor readiness"""
    assert predictor.is_ready() == True

def test_batch_predict(predictor):
    """Test batch prediction functionality"""
    compounds = [
        {"smiles": "CCN", "compound_id": "comp1"},
        {"smiles": "CCO", "compound_id": "comp2"}
    ]

    results = predictor.batch_predict(compounds)

    assert len(results) == 2
    for result in results:
        assert "prediction" in result or "error" in result