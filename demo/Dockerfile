FROM continuumio/miniconda3:latest

WORKDIR /code

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Install RDKit and scientific packages via conda (most stable approach)
RUN conda install -c conda-forge rdkit pillow numpy pandas -y \
    && conda clean --all -y

# Copy requirements and install Python packages
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY . .

# Expose Gradio port
EXPOSE 7860

# Run the Gradio application
CMD ["python", "app.py"]