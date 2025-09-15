"""
Basic API tests for ARR-MEDIC CYP3A4 Opensource
"""

import pytest
from fastapi.testclient import TestClient
from main import app

client = TestClient(app)


def test_health_endpoint():
    """Test health check endpoint"""
    response = client.get("/health")
    assert response.status_code == 200
    data = response.json()
    assert "status" in data
    assert data["status"] in ["healthy", "unhealthy"]
    assert "timestamp" in data
    assert "version" in data


def test_root_endpoint():
    """Test root endpoint"""
    response = client.get("/")
    assert response.status_code == 200
    data = response.json()
    assert "name" in data
    assert "version" in data
    assert data["name"] == "ARR-MEDIC CYP3A4 Opensource"


def test_predict_basic():
    """Test basic prediction endpoint"""
    payload = {
        "smiles": "CCO",  # Ethanol
        "compound_id": "ethanol_test"
    }
    response = client.post("/predict", json=payload)
    assert response.status_code == 200
    data = response.json()

    # Check required fields
    assert "prediction" in data
    assert "probability" in data
    assert "confidence" in data
    assert "smiles" in data
    assert data["smiles"] == "CCO"
    assert data["compound_id"] == "ethanol_test"


def test_predict_invalid_smiles():
    """Test prediction with invalid SMILES"""
    payload = {
        "smiles": "",  # Empty SMILES
        "compound_id": "test_invalid"
    }
    response = client.post("/predict", json=payload)
    assert response.status_code == 422  # Validation error


def test_predict_batch():
    """Test batch prediction endpoint"""
    payload = {
        "compounds": [
            {"smiles": "CCO", "compound_id": "ethanol"},
            {"smiles": "CC(=O)O", "compound_id": "acetic_acid"}
        ]
    }
    response = client.post("/predict/batch", json=payload)
    assert response.status_code == 200
    data = response.json()

    # Check batch response structure
    assert "results" in data
    assert "total_processed" in data
    assert "successful" in data
    assert "failed" in data
    assert len(data["results"]) <= 2  # Should have results for valid compounds


def test_docs_endpoint():
    """Test API documentation endpoint"""
    response = client.get("/docs")
    assert response.status_code == 200


def test_openapi_json():
    """Test OpenAPI JSON schema"""
    response = client.get("/openapi.json")
    assert response.status_code == 200
    data = response.json()
    assert "openapi" in data
    assert "info" in data