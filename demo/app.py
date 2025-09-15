"""
ARR-MEDIC CYP3A4 Opensource - Hugging Face Spaces Demo
Interactive CYP3A4 inhibition prediction with molecular visualization
"""

import gradio as gr
import sys
import os
import traceback
from typing import Dict, Any, Optional

# Add backend path to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'backend'))

try:
    from predictor import CYP3A4BasicPredictor
    PREDICTOR_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import predictor: {e}")
    PREDICTOR_AVAILABLE = False

# Optional RDKit import for molecular visualization
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    import io
    import base64
    from PIL import Image
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available. Molecular visualization will be limited.")

# Initialize predictor
predictor = None
if PREDICTOR_AVAILABLE:
    try:
        predictor = CYP3A4BasicPredictor()
    except Exception as e:
        print(f"Failed to initialize predictor: {e}")

def create_molecule_image(smiles: str) -> Optional[Image.Image]:
    """Create molecular structure image from SMILES"""
    if not RDKIT_AVAILABLE or not smiles.strip():
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Generate 2D coordinates
        from rdkit.Chem import rdDepictor
        rdDepictor.Compute2DCoords(mol)

        # Create image
        img = Draw.MolToImage(mol, size=(300, 300))
        return img
    except Exception as e:
        print(f"Error creating molecule image: {e}")
        return None

def predict_cyp3a4(smiles: str, compound_id: str = "") -> tuple:
    """
    Main prediction function for Gradio interface

    Returns:
        tuple: (prediction_text, confidence_text, molecular_image, warning_text)
    """

    if not smiles.strip():
        return (
            "❌ **오류**: SMILES 문자열을 입력해주세요",
            "",
            None,
            ""
        )

    # Create molecular image
    mol_image = create_molecule_image(smiles)

    # Make prediction
    if not predictor:
        return (
            "⚠️ **예측 불가**: 예측기를 초기화할 수 없습니다",
            "",
            mol_image,
            "시스템 오류로 인해 예측을 수행할 수 없습니다."
        )

    try:
        # Run prediction
        result = predictor.predict(smiles, compound_id or f"compound_{hash(smiles)}")

        # Format prediction result
        prediction = result['prediction']
        probability = result['probability']
        confidence = result['confidence']
        risk_level = result['risk_level']

        # Create result text
        if prediction == "inhibitor":
            pred_text = f"🔴 **CYP3A4 억제제** (확률: {probability:.1%})"
            risk_color = "🔴" if risk_level == "high" else "🟡"
        else:
            pred_text = f"🟢 **비억제제** (확률: {(1-probability):.1%})"
            risk_color = "🟢"

        confidence_text = f"**신뢰도**: {confidence:.1%} | **위험도**: {risk_color} {risk_level.title()}"

        # Format molecular descriptors
        descriptors = result.get('descriptors', {})
        desc_lines = []
        if descriptors:
            desc_lines.append("### 📊 분자 특성")
            desc_lines.append(f"- **분자량**: {descriptors.get('molecular_weight', 'N/A'):.1f}")
            desc_lines.append(f"- **LogP**: {descriptors.get('logp', 'N/A'):.2f}")
            desc_lines.append(f"- **TPSA**: {descriptors.get('tpsa', 'N/A'):.1f}")
            desc_lines.append(f"- **H-결합 공여체**: {descriptors.get('hbd', 'N/A')}")
            desc_lines.append(f"- **H-결합 수용체**: {descriptors.get('hba', 'N/A')}")
            desc_lines.append(f"- **회전 결합수**: {descriptors.get('rotatable_bonds', 'N/A')}")

        descriptor_text = "\n".join(desc_lines) if desc_lines else ""

        # Format warnings
        warnings = result.get('warnings', [])
        warning_text = ""
        if warnings:
            warning_text = "### ⚠️ 주의사항\n" + "\n".join([f"- {w}" for w in warnings])

        full_result = f"{pred_text}\n\n{descriptor_text}"

        return (
            full_result,
            confidence_text,
            mol_image,
            warning_text
        )

    except Exception as e:
        error_text = f"❌ **예측 오류**: {str(e)}"
        return (
            error_text,
            "",
            mol_image,
            "예측 중 오류가 발생했습니다. SMILES 문자열을 확인해주세요."
        )

def get_example_compounds():
    """Return example compounds for quick testing"""
    return [
        ["CCO", "에탄올"],
        ["CCN", "에틸아민"],
        ["CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "이부프로펜"],
        ["CC1=CC=C(C=C1)C2=CC(=NN2C3=CC=C(C=C3)S(=O)(=O)N)C(F)(F)F", "셀레콕시브"],
        ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "카페인"],
    ]

# Create Gradio interface
def create_interface():
    with gr.Blocks(
        title="ARR-MEDIC CYP3A4 Demo",
        theme=gr.themes.Soft(),
        css="""
        .gradio-container {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }
        .warning-text {
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 4px;
            padding: 10px;
            margin: 10px 0;
        }
        """
    ) as demo:

        gr.Markdown("""
        # 🧬 ARR-MEDIC CYP3A4 Opensource Demo

        **CYP3A4 약물 상호작용 예측 시스템** - 연구 및 교육 목적

        🔴 **주의**: 이 도구는 연구 및 교육 목적으로만 사용하세요. 임상 진단이나 치료 결정에 사용하지 마세요.
        """)

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### 📥 입력")

                smiles_input = gr.Textbox(
                    label="SMILES 문자열",
                    placeholder="예: CCO (에탄올)",
                    info="화합물의 SMILES 표기법을 입력하세요",
                    lines=2
                )

                compound_id = gr.Textbox(
                    label="화합물 ID (선택사항)",
                    placeholder="예: compound_001",
                    info="식별을 위한 화합물 이름이나 ID"
                )

                predict_btn = gr.Button(
                    "🔍 CYP3A4 예측 실행",
                    variant="primary",
                    size="lg"
                )

                gr.Markdown("### 🧪 예시 화합물")

                examples = get_example_compounds()
                for smiles, name in examples:
                    gr.Button(
                        f"{name} ({smiles})",
                        size="sm",
                        variant="secondary"
                    ).click(
                        lambda s=smiles: s,
                        outputs=smiles_input
                    )

            with gr.Column(scale=2):
                gr.Markdown("### 📊 예측 결과")

                with gr.Row():
                    with gr.Column():
                        prediction_output = gr.Markdown(
                            label="예측 결과",
                            value="결과가 여기에 표시됩니다..."
                        )

                        confidence_output = gr.Markdown(
                            label="신뢰도 정보"
                        )

                    with gr.Column():
                        molecule_image = gr.Image(
                            label="분자 구조",
                            show_label=True,
                            height=300
                        )

                warning_output = gr.Markdown(
                    label="주의사항",
                    elem_classes=["warning-text"]
                )

        # Button click event
        predict_btn.click(
            predict_cyp3a4,
            inputs=[smiles_input, compound_id],
            outputs=[prediction_output, confidence_output, molecule_image, warning_output]
        )

        # Footer information
        gr.Markdown("""
        ---

        ### 📚 정보

        - **정확도**: ~70% (교육용 기본 모델)
        - **방법**: 규칙 기반 분자 기술자 분석
        - **라이선스**: MIT License
        - **GitHub**: [ARR-MEDIC CYP3A4 Opensource](https://github.com/your-org/arr-medic-cyp3a4-opensource)

        ### 🔬 CYP3A4란?

        CYP3A4는 인간의 간에서 가장 중요한 약물 대사 효소입니다. 많은 의약품들이 이 효소에 의해 대사되므로,
        CYP3A4를 억제하는 화합물은 다른 약물과 상호작용을 일으킬 수 있습니다.

        **⚠️ 면책조항**: 이 예측은 연구 및 교육 목적으로만 제공됩니다.
        임상적 결정이나 규제 승인에 사용하지 마세요.
        """)

    return demo

# Launch demo
if __name__ == "__main__":
    demo = create_interface()
    demo.launch(
        server_name="0.0.0.0",
        server_port=7860,
        share=False,
        show_error=True
    )