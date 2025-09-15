"""
ARR-MEDIC CYP3A4 Opensource - Multilingual Demo
Korean-English language toggle with enhanced UI
"""

import gradio as gr
import sys
import os
from typing import Dict, Any, Optional

# Add backend path
sys.path.append('.')

try:
    from predictor import CYP3A4BasicPredictor
    PREDICTOR_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import predictor: {e}")
    PREDICTOR_AVAILABLE = False

# Optional RDKit import
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from PIL import Image
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("Warning: RDKit not available.")

# Initialize predictor
predictor = None
if PREDICTOR_AVAILABLE:
    try:
        predictor = CYP3A4BasicPredictor()
    except Exception as e:
        print(f"Failed to initialize predictor: {e}")

# Language definitions
LANGUAGES = {
    'en': {
        'title': '🧬 ARR-MEDIC CYP3A4 Interactive Demo',
        'disclaimer': '🔴 **Not for Clinical Use** - Research and educational purposes only',
        'smiles_label': 'SMILES String',
        'smiles_placeholder': 'e.g., CCO (ethanol)',
        'smiles_info': 'Enter the SMILES notation of your compound',
        'compound_label': 'Compound ID (optional)',
        'compound_placeholder': 'e.g., compound_001',
        'predict_button': '🔍 Run CYP3A4 Prediction',
        'results_title': '📊 Prediction Results',
        'examples_title': '🧪 Example Compounds',
        'error_no_smiles': '❌ **Error**: Please enter a SMILES string',
        'waiting_results': 'Results will appear here...',
        'about_title': '📚 About CYP3A4',
        'about_text': '''CYP3A4 is the most important drug-metabolizing enzyme:
- **Metabolizes 50%+ of drugs** - Critical for drug safety
- **Drug interactions** - Can cause dangerous combinations
- **Educational tool** - Learn molecular property relationships''',
        'ethanol_btn': 'Ethanol (CCO)',
        'caffeine_btn': 'Caffeine',
        'ibuprofen_btn': 'Ibuprofen',
    },
    'ko': {
        'title': '🧬 ARR-MEDIC CYP3A4 인터랙티브 데모',
        'disclaimer': '🔴 **임상 사용 금지** - 연구 및 교육 목적으로만 사용',
        'smiles_label': 'SMILES 문자열',
        'smiles_placeholder': '예: CCO (에탄올)',
        'smiles_info': '화합물의 SMILES 표기법을 입력하세요',
        'compound_label': '화합물 ID (선택사항)',
        'compound_placeholder': '예: compound_001',
        'predict_button': '🔍 CYP3A4 예측 실행',
        'results_title': '📊 예측 결과',
        'examples_title': '🧪 예시 화합물',
        'error_no_smiles': '❌ **오류**: SMILES 문자열을 입력해주세요',
        'waiting_results': '결과가 여기에 표시됩니다...',
        'about_title': '📚 CYP3A4에 대하여',
        'about_text': '''CYP3A4는 가장 중요한 약물 대사 효소입니다:
- **50% 이상 약물 대사** - 약물 안전성에 핵심
- **약물 상호작용** - 위험한 조합을 일으킬 수 있음
- **교육 도구** - 분자 특성 관계 학습''',
        'ethanol_btn': '에탄올 (CCO)',
        'caffeine_btn': '카페인',
        'ibuprofen_btn': '이부프로펜',
    }
}

def get_text(lang: str, key: str) -> str:
    return LANGUAGES.get(lang, LANGUAGES['en']).get(key, key)

def create_molecule_image(smiles: str) -> Optional[Image.Image]:
    if not RDKIT_AVAILABLE or not smiles.strip():
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        from rdkit.Chem import rdDepictor
        rdDepictor.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(300, 300))
        return img
    except Exception as e:
        print(f"Error creating molecule image: {e}")
        return None

def predict_cyp3a4(smiles: str, compound_id: str, language: str) -> tuple:
    if not smiles.strip():
        return (
            get_text(language, 'error_no_smiles'),
            "",
            None,
            ""
        )

    mol_image = create_molecule_image(smiles)

    if not predictor:
        error_msg = "⚠️ **예측 불가**: 예측기 초기화 실패" if language == 'ko' else "⚠️ **Prediction Error**: Predictor initialization failed"
        return (error_msg, "", mol_image, "")

    try:
        result = predictor.predict(smiles, compound_id or f"compound_{hash(smiles)}")

        prediction = result['prediction']
        probability = result['probability']
        confidence = result['confidence']
        risk_level = result['risk_level']

        if language == 'ko':
            if prediction == "inhibitor":
                pred_text = f"🔴 **CYP3A4 억제제** (확률: {probability:.1%})"
            else:
                pred_text = f"🟢 **비억제제** (확률: {(1-probability):.1%})"

            confidence_text = f"**신뢰도**: {confidence:.1%} | **위험도**: {risk_level}"

            descriptors = result.get('descriptors', {})
            if descriptors:
                desc_text = f"""

### 📊 분자 특성
- **분자량**: {descriptors.get('molecular_weight', 'N/A'):.1f} Da
- **LogP**: {descriptors.get('logp', 'N/A'):.2f}
- **TPSA**: {descriptors.get('tpsa', 'N/A'):.1f} Ų
- **H-결합 공여체**: {descriptors.get('hbd', 'N/A')}개
- **H-결합 수용체**: {descriptors.get('hba', 'N/A')}개
"""
            else:
                desc_text = ""

        else:
            if prediction == "inhibitor":
                pred_text = f"🔴 **CYP3A4 Inhibitor** (Probability: {probability:.1%})"
            else:
                pred_text = f"🟢 **Non-Inhibitor** (Probability: {(1-probability):.1%})"

            confidence_text = f"**Confidence**: {confidence:.1%} | **Risk Level**: {risk_level}"

            descriptors = result.get('descriptors', {})
            if descriptors:
                desc_text = f"""

### 📊 Molecular Properties
- **Molecular Weight**: {descriptors.get('molecular_weight', 'N/A'):.1f} Da
- **LogP**: {descriptors.get('logp', 'N/A'):.2f}
- **TPSA**: {descriptors.get('tpsa', 'N/A'):.1f} Ų
- **H-Bond Donors**: {descriptors.get('hbd', 'N/A')}
- **H-Bond Acceptors**: {descriptors.get('hba', 'N/A')}
"""
            else:
                desc_text = ""

        full_result = f"{pred_text}{desc_text}"

        warnings = result.get('warnings', [])
        warning_text = ""
        if warnings:
            if language == 'ko':
                warning_text = "### ⚠️ 주의사항\n" + "\n".join([f"- {w}" for w in warnings])
            else:
                warning_text = "### ⚠️ Important Notes\n" + "\n".join([f"- {w}" for w in warnings])

        return (full_result, confidence_text, mol_image, warning_text)

    except Exception as e:
        error_text = f"❌ **예측 오류**: {str(e)}" if language == 'ko' else f"❌ **Prediction Error**: {str(e)}"
        return (error_text, "", mol_image, "")

def update_interface(language):
    return [
        get_text(language, 'title'),
        get_text(language, 'disclaimer'),
        gr.update(
            label=get_text(language, 'smiles_label'),
            placeholder=get_text(language, 'smiles_placeholder'),
            info=get_text(language, 'smiles_info')
        ),
        gr.update(
            label=get_text(language, 'compound_label'),
            placeholder=get_text(language, 'compound_placeholder')
        ),
        gr.update(value=get_text(language, 'predict_button')),
        get_text(language, 'results_title'),
        get_text(language, 'examples_title'),
        get_text(language, 'about_title'),
        get_text(language, 'about_text'),
        get_text(language, 'waiting_results'),
        gr.update(value=get_text(language, 'ethanol_btn')),
        gr.update(value=get_text(language, 'caffeine_btn')),
        gr.update(value=get_text(language, 'ibuprofen_btn')),
    ]

# Create interface
with gr.Blocks(title="ARR-MEDIC CYP3A4 Demo", theme=gr.themes.Soft()) as demo:

    # Language selector
    with gr.Row():
        language_selector = gr.Radio(
            choices=[("English", "en"), ("한국어", "ko")],
            value="en",
            label="🌐 Language / 언어",
            container=False
        )

    # Main content
    title_md = gr.Markdown("🧬 ARR-MEDIC CYP3A4 Interactive Demo")
    disclaimer_md = gr.Markdown("🔴 **Not for Clinical Use** - Research and educational purposes only")

    with gr.Row():
        with gr.Column():
            gr.Markdown("### 📥 Input")

            smiles_input = gr.Textbox(
                label="SMILES String",
                placeholder="e.g., CCO (ethanol)",
                info="Enter the SMILES notation of your compound"
            )

            compound_input = gr.Textbox(
                label="Compound ID (optional)",
                placeholder="e.g., compound_001"
            )

            predict_btn = gr.Button(
                "🔍 Run CYP3A4 Prediction",
                variant="primary",
                size="lg"
            )

            examples_title = gr.Markdown("### 🧪 Example Compounds")

            with gr.Row():
                ethanol_btn = gr.Button("Ethanol (CCO)", size="sm")
                caffeine_btn = gr.Button("Caffeine", size="sm")
                ibuprofen_btn = gr.Button("Ibuprofen", size="sm")

        with gr.Column():
            results_title = gr.Markdown("### 📊 Prediction Results")

            with gr.Row():
                with gr.Column():
                    prediction_output = gr.Markdown("Results will appear here...")
                    confidence_output = gr.Markdown("")

                with gr.Column():
                    molecule_image = gr.Image(label="Molecular Structure", height=300)

            warning_output = gr.Markdown("")

    # Information section
    with gr.Row():
        with gr.Column():
            about_title = gr.Markdown("### 📚 About CYP3A4")
            about_text = gr.Markdown("""CYP3A4 is the most important drug-metabolizing enzyme:
- **Metabolizes 50%+ of drugs** - Critical for drug safety
- **Drug interactions** - Can cause dangerous combinations
- **Educational tool** - Learn molecular property relationships""")

    # Event handlers
    predict_btn.click(
        predict_cyp3a4,
        inputs=[smiles_input, compound_input, language_selector],
        outputs=[prediction_output, confidence_output, molecule_image, warning_output]
    )

    # Example buttons
    ethanol_btn.click(lambda: "CCO", outputs=smiles_input)
    caffeine_btn.click(lambda: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", outputs=smiles_input)
    ibuprofen_btn.click(lambda: "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", outputs=smiles_input)

    # Language change
    language_selector.change(
        update_interface,
        inputs=[language_selector],
        outputs=[
            title_md, disclaimer_md, smiles_input, compound_input, predict_btn,
            results_title, examples_title, about_title, about_text, prediction_output,
            ethanol_btn, caffeine_btn, ibuprofen_btn
        ]
    )

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860, share=False)