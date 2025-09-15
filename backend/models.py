"""
Pydantic Models for ARR-MEDIC CYP3A4 Opensource API

Data validation and serialization models.
"""

from pydantic import BaseModel, Field, validator
from typing import Dict, List, Any, Optional
from datetime import datetime
from enum import Enum

class StrictBase(BaseModel):
    class Config:
        extra = "forbid"  # disallow unknown fields

class PredictionType(str, Enum):
    INHIBITOR = "inhibitor"
    NON_INHIBITOR = "non_inhibitor"

class RiskLevel(str, Enum):
    LOW = "low"
    MODERATE = "moderate"
    HIGH = "high"

class PredictionRequest(StrictBase):
    """Request model for CYP3A4 prediction"""
    smiles: str = Field(..., description="SMILES string of the compound")
    compound_id: Optional[str] = Field(None, description="Optional compound identifier")
    
    @validator('smiles')
    def validate_smiles(cls, v):
        if not v or len(v.strip()) == 0:
            raise ValueError("SMILES string cannot be empty")
        if len(v) > 1000:
            raise ValueError("SMILES string too long (max 1000 characters)")
        return v.strip()

class PredictionResponse(StrictBase):
    """Response model for CYP3A4 prediction"""
    compound_id: Optional[str]
    smiles: str
    prediction: PredictionType
    probability: float = Field(..., ge=0.0, le=1.0, description="Prediction probability")
    confidence: float = Field(..., ge=0.0, le=1.0, description="Model confidence")
    risk_level: RiskLevel
    timestamp: datetime
    model_version: str = "basic-v1.0"
    processing_time: float = Field(..., description="Processing time in seconds")
    molecular_descriptors: Dict[str, float] = Field(default_factory=dict)
    warnings: List[str] = Field(default_factory=list)

class BatchPredictionRequest(StrictBase):
    """Request model for batch predictions"""
    compounds: List[Dict[str, str]] = Field(
        ..., 
        description="List of compounds with 'smiles' and optional 'compound_id'"
    )
    
    @validator('compounds')
    def validate_compounds(cls, v):
        if not v:
            raise ValueError("Compounds list cannot be empty")
        if len(v) > 100:
            raise ValueError("Batch size limited to 100 compounds")
        return v

class BatchPredictionResponse(StrictBase):
    """Response model for batch predictions"""
    results: List[PredictionResponse]
    total_processed: int
    successful: int
    failed: int
    processing_time: float

class HealthResponse(StrictBase):
    """Health check response model"""
    status: str = Field(..., description="Overall system status")
    timestamp: datetime
    version: str
    components: Dict[str, str] = Field(default_factory=dict)
    error: Optional[str] = None

class HistoryItem(StrictBase):
    """Historical prediction record"""
    id: Optional[int]
    compound_id: Optional[str] 
    smiles: str
    prediction: PredictionType
    probability: float
    confidence: float
    timestamp: datetime
    processing_time: float

class ApiStats(StrictBase):
    """API usage statistics"""
    version: str
    model_type: str
    accuracy: str
    uptime: str
    total_predictions: Optional[int] = 0
    features: List[str] = Field(default_factory=list)
    limitations: List[str] = Field(default_factory=list)