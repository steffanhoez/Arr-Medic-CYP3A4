import gradio as gr
import sys
import os

# Add the backend directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), 'backend')))

try:
    # 백엔드 폴더에서 예측기 클래스를 가져옵니다.
    from predictor import CYP3A4BasicPredictor
    # 예측기를 한 번만 초기화합니다.
    predictor = CYP3A4BasicPredictor()
    predictor_ready = True
except Exception as e:
    predictor_ready = False
    # 초기화 실패 시 오류 메시지를 저장합니다.
    predictor_error = f"Failed to load predictor. Is 'rdkit' installed correctly? Error: {e}"

def predict_cyp3a4(smiles_string: str):
    """
    Gradio 인터페이스를 위한 래퍼 함수입니다.
    """
    if not predictor_ready:
        return {"error": predictor_error}
    
    if not smiles_string or not smiles_string.strip():
        return {"error": "Please enter a valid SMILES string."}
        
    try:
        # 예측기 실행
        result = predictor.predict(smiles_string.strip())
        return result
    except Exception as e:
        return {"error": f"An error occurred during prediction: {str(e)}"}

# Gradio UI 정의
demo = gr.Interface(
    fn=predict_cyp3a4,
    inputs=gr.Textbox(label="SMILES String", placeholder="e.g., CCO for Ethanol"),
    outputs=gr.JSON(label="Prediction Result"),
    title="ARR-MEDIC CYP3A4 Inhibitor Prediction",
    description="""
    화학 화합물이 CYP3A4 효소를 억제하는지 예측하는 데모입니다.
    SMILES 문자열을 입력하고 결과를 확인하세요.
    **본 데모는 학술 및 연구 목적으로만 사용해야 합니다.**
    """,
    examples=[
        ["CCO"],  # Ethanol
        ["CC(=O)OC1=CC=CC=C1C(=O)O"], # Aspirin
        ["CCN"], # Ethylamine
    ],
    allow_flagging="never"
)

if __name__ == "__main__":
    demo.launch()
