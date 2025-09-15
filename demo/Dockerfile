FROM continuumio/miniconda3:latest

WORKDIR /code

# Update conda and install scientific packages via conda (most stable)
RUN conda update -n base -c defaults conda && \
    conda install -c conda-forge \
    python=3.10 \
    rdkit \
    pillow \
    numpy \
    pandas \
    pydantic \
    python-dotenv \
    -y && \
    conda clean --all -y

# Install Gradio and joblib via pip (minimal pip usage for compatibility)
RUN pip install --no-cache-dir gradio~=3.50 joblib

# Copy application files
COPY . .

# Expose Gradio port
EXPOSE 7860

# Run the multilingual Gradio application
CMD ["python", "app.py"]