"""
ARR-MEDIC CYP3A4 Opensource - Main FastAPI Application

Basic CYP3A4 prediction service for academic research and learning.
MIT License - Educational purposes only.
"""

from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.middleware.gzip import GZipMiddleware
from contextlib import asynccontextmanager
import uvicorn
import logging
from typing import List, Dict, Any
import os
from datetime import datetime, timezone
from dotenv import load_dotenv

from .models import (
    PredictionRequest, PredictionResponse, HealthResponse,
    BatchPredictionRequest, BatchPredictionResponse
)
from .predictor import CYP3A4BasicPredictor
from .database import get_database, Database, close_database

# Load environment variables
load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Get CORS origins from environment
_origins = os.getenv("CORS_ALLOW_ORIGINS", "").strip()
ALLOWED_ORIGINS = [o.strip() for o in _origins.split(",") if o.strip()] or ["http://localhost:3000"]

# Get allowed hosts from environment
_allowed_hosts = os.getenv("FASTAPI_ALLOWED_HOSTS", "").strip()
ALLOWED_HOSTS = [h.strip() for h in _allowed_hosts.split(",") if h.strip()] or [
    "localhost",
    "127.0.0.1",
    "testserver",  # Allow TestClient requests
    "your-domain.com",
]

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan manager"""
    # Startup
    logger.info("Starting ARR-MEDIC CYP3A4 Opensource API...")
    db = await get_database()
    if not db:
        logger.error("Database initialization failed")
        raise RuntimeError("Database initialization failed")
    logger.info("Application startup complete")

    yield

    # Shutdown
    logger.info("Shutting down ARR-MEDIC CYP3A4 Opensource API...")
    await close_database()
    logger.info("Application shutdown complete")

# Initialize FastAPI app
app = FastAPI(
    title="ARR-MEDIC CYP3A4 Opensource",
    description="Basic CYP3A4 Drug Interaction Prediction API",
    version="1.0.0-opensource",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)

# Security middleware
app.add_middleware(GZipMiddleware, minimum_size=1024)
app.add_middleware(TrustedHostMiddleware, allowed_hosts=ALLOWED_HOSTS)

# CORS middleware with environment-based origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global exception handler (clean JSON)
from starlette.requests import Request

@app.exception_handler(Exception)
async def _unhandled(request: Request, exc: Exception):
    logger.exception("Unhandled error: %s", exc)
    return JSONResponse(status_code=500, content={"detail": "Internal Server Error"})

# Initialize predictor
predictor = CYP3A4BasicPredictor()

# Lifespan events handled in lifespan context manager above

@app.get("/", response_model=Dict[str, str])
async def root():
    """Root endpoint - API information"""
    return {
        "name": "ARR-MEDIC CYP3A4 Opensource",
        "version": "1.0.0",
        "description": "Basic CYP3A4 prediction for research and education",
        "license": "MIT",
        "docs": "/docs"
    }

@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint"""
    try:
        # Basic health checks
        db_status = "healthy"  # Simplified for opensource
        model_status = "loaded" if predictor.is_ready() else "loading"
        
        return HealthResponse(
            status="healthy",
            timestamp=datetime.now(timezone.utc),
            version="1.0.0-opensource",
            components={
                "database": db_status,
                "ml_model": model_status,
                "api": "healthy"
            }
        )
    except Exception as e:
        logger.error(f"Health check failed: {e}")
        return HealthResponse(
            status="unhealthy",
            timestamp=datetime.now(timezone.utc),
            version="1.0.0-opensource",
            components={
                "database": "unknown",
                "ml_model": "unknown", 
                "api": "error"
            },
            error=str(e)
        )

@app.post("/predict", response_model=PredictionResponse)
async def predict_cyp3a4(
    request: PredictionRequest,
    db: Database = Depends(get_database)
):
    """
    Predict CYP3A4 inhibition for given compound
    
    Basic prediction with ~70% accuracy for research purposes.
    """
    try:
        logger.info(f"Processing prediction request for compound: {request.compound_id}")
        
        # Validate input
        if not request.smiles:
            raise HTTPException(status_code=400, detail="SMILES string is required")
        
        # Make prediction
        result = predictor.predict(
            smiles=request.smiles,
            compound_id=request.compound_id
        )
        
        # Store in database (simplified)
        if db:
            await db.store_prediction(request, result)
        
        # Format response
        response = PredictionResponse(
            compound_id=request.compound_id,
            smiles=request.smiles,
            prediction=result["prediction"],
            confidence=result["confidence"],
            probability=result["probability"],
            risk_level=result["risk_level"],
            timestamp=datetime.now(timezone.utc),
            model_version="basic-v1.0",
            processing_time=result.get("processing_time", 0.0),
            molecular_descriptors=result.get("descriptors", {}),
            warnings=result.get("warnings", [])
        )
        
        logger.info(f"Prediction completed for {request.compound_id}: {result['prediction']}")
        return response

    except ValueError as e:
        logger.error(f"Invalid input: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")

@app.post("/predict/batch", response_model=BatchPredictionResponse)
async def predict_batch(
    payload: BatchPredictionRequest,
    db: Database = Depends(get_database)
):
    """Batch prediction endpoint for multiple compounds"""
    try:
        logger.info(f"Processing batch of {len(payload.compounds)} compounds")
        start_time = datetime.now(timezone.utc)
        results = []
        successful = 0
        failed = 0

        for compound in payload.compounds:
            try:
                # Validate compound data
                if "smiles" not in compound:
                    failed += 1
                    continue

                # Make prediction
                result = predictor.predict(
                    smiles=compound["smiles"],
                    compound_id=compound.get("compound_id")
                )

                # Store in database
                if db:
                    request_obj = type("Request", (), {
                        "smiles": compound["smiles"],
                        "compound_id": compound.get("compound_id")
                    })()
                    await db.store_prediction(request_obj, result)

                # Create response
                response = PredictionResponse(
                    compound_id=compound.get("compound_id"),
                    smiles=compound["smiles"],
                    prediction=result["prediction"],
                    probability=result["probability"],
                    confidence=result["confidence"],
                    risk_level=result["risk_level"],
                    timestamp=datetime.now(timezone.utc),
                    model_version="basic-v1.0",
                    processing_time=result.get("processing_time", 0.0),
                    molecular_descriptors=result.get("descriptors", {}),
                    warnings=result.get("warnings", [])
                )

                results.append(response)
                successful += 1

            except Exception as e:
                logger.error(f"Failed to process compound {compound.get('compound_id', 'unknown')}: {e}")
                failed += 1
                continue

        processing_time = (datetime.now(timezone.utc) - start_time).total_seconds()

        batch_response = BatchPredictionResponse(
            results=results,
            total_processed=len(payload.compounds),
            successful=successful,
            failed=failed,
            processing_time=processing_time
        )

        logger.info(f"Batch processing completed: {successful} successful, {failed} failed")
        return batch_response

    except Exception as e:
        logger.error(f"Batch prediction failed: {e}")
        raise HTTPException(status_code=500, detail="Batch prediction failed")

@app.get("/history")
async def get_prediction_history(
    limit: int = 50,
    db: Database = Depends(get_database)
):
    """Get recent prediction history"""
    try:
        if not db:
            return {"message": "Database not available", "predictions": []}
        
        history = await db.get_recent_predictions(limit=limit)
        return {
            "predictions": history,
            "total": len(history),
            "limit": limit
        }
    except Exception as e:
        logger.error(f"Failed to get history: {e}")
        raise HTTPException(status_code=500, detail="Failed to retrieve history")

@app.get("/stats")
async def get_api_stats():
    """Get basic API usage statistics"""
    return {
        "version": "1.0.0-opensource",
        "model_type": "basic",
        "accuracy": "~70%",
        "uptime": "available",
        "features": [
            "Basic CYP3A4 prediction",
            "Molecular descriptor analysis",
            "Simple risk assessment",
            "REST API access"
        ],
        "limitations": [
            "No clinical validation",
            "No EMR integration", 
            "No advanced ML models",
            "Research use only"
        ]
    }

if __name__ == "__main__":
    port = int(os.getenv("PORT", 8000))
    host = os.getenv("HOST", "0.0.0.0")
    
    logger.info("Starting ARR-MEDIC CYP3A4 Opensource API")
    uvicorn.run(app, host=host, port=port)