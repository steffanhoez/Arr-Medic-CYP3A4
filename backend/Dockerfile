FROM python:3.10-slim

WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt /app/requirements.txt

# Install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# rdkit may need additional system packages
# RUN apt-get update && apt-get install -y libxrender1 && rm -rf /var/lib/apt/lists/*

# Copy application code
COPY . /app/

# Environment variables
ENV HOST=0.0.0.0
ENV PORT=8000
ENV DB_PATH=/app/arr_opensource.db

# Expose port
EXPOSE 8000

# Run the application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]