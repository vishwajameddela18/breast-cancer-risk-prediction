import streamlit as st
import pandas as pd
import joblib
import requests

# ----- Load Dataset and Trained Model -----
df = pd.read_excel("mmc4.xlsx")
df = df.dropna(subset=["Chromosome", "Positionb", "EAFc", "Overall Breast Cancerd"])
model = joblib.load("model.pkl")

# ----- Configure Streamlit Page -----
st.set_page_config(page_title="Breast Cancer SNP Risk Predictor", layout="wide")
st.title("🧬 Breast Cancer SNP Risk Prediction Tool")
st.caption(
    "This tool analyzes a single genetic variant (SNP) to estimate breast cancer risk based on known associations. "
    "It combines genomic data with predictive modeling to support clinical decision-making."
)

# ----- Sample Selection Interface -----
st.write(f"Available Samples: 0 to {len(df)-1}")
sample_idx = st.number_input(
    "🔎 Select a sample index from the dataset:",
    min_value=0,
    max_value=len(df) - 1,
    value=0,
    step=1
)
sample = df.iloc[sample_idx]

# ----- Model Prediction -----
X = sample[["Chromosome", "Positionb", "EAFc", "ER-positivee", "ER-negativef"]].values.reshape(1, -1)
pred_score = model.predict(X)[0]

def label_risk(score):
    if score >= 0.01:
        return "🟥 High"
    elif score > 0:
        return "🟧 Moderate"
    else:
        return "🟩 Low"
risk_label = label_risk(pred_score)

# ----- Map Chromosome/Position to rsID -----
@st.cache_data(show_spinner=False)
def map_to_rsid(chrom, pos):
    url = f"https://grch37.rest.ensembl.org/overlap/region/human/{int(chrom)}:{int(pos)}-{int(pos)}?feature=variation"
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=10)
        if r.status_code == 200 and r.json():
            return r.json()[0]["id"]
    except:
        return None
    return None

rsid = map_to_rsid(sample["Chromosome"], sample["Positionb"])

# ----- Fetch Variant Annotation from Ensembl -----
def get_variant_info(rsid):
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?"
    headers = {"Content-Type": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=10)
        if r.status_code != 200:
            return None
        data = r.json()
        tc = data[0].get('transcript_consequences', [{}])[0]
        return {
            "rsID": rsid,
            "Gene": tc.get('gene_symbol', 'Unknown'),
            "Impact": tc.get('impact', 'Unknown'),
            "Effect": tc.get('consequence_terms', ['Unknown'])[0],
            "Clinical Significance": data[0].get('clinical_significance', ['Not reported'])[0]
        }
    except:
        return None

variant = get_variant_info(rsid) if rsid else None

# ----- Display Prediction and Variant Info -----
col1, col2 = st.columns(2)

with col1:
    st.subheader("📈 Model Risk Prediction")
    st.metric(label="Predicted Risk Score", value=f"{pred_score:.6f}")
    st.metric(label="Risk Level", value=risk_label)
    st.caption("The predicted score reflects the statistical association of this variant with breast cancer risk.")

with col2:
    st.subheader("🧬 Genomic Variant Details")
    if rsid:
        st.write(f"**rsID**: {rsid}")
        if variant:
            st.write(f"**Gene**: {variant['Gene']}")
            st.write(f"**Effect**: {variant['Effect']} (e.g., intron variant, missense)")
            st.write(f"**Impact**: {variant['Impact']} (predicted functional impact on gene)")
            st.write(f"**Clinical Significance**: {variant['Clinical Significance']}")
        else:
            st.warning("No detailed variant annotation available.")
    else:
        st.info("No rsID found for this chromosomal position.")

# ----- Recommendation Based on Clinical Context -----
st.divider()
st.subheader("🏥 Clinical Recommendation")

def get_recommendation(impact, clin_sig):
    clin_sig = (clin_sig or "").lower()

    if impact == "HIGH" or clin_sig in ["pathogenic", "likely_pathogenic"]:
        return (
            "🔴 This variant is classified as high impact or clinically pathogenic. "
            "Referral to a clinical geneticist is recommended. Consider confirmatory testing, "
            "comprehensive family history assessment, and risk-reducing strategies based on guidelines."
        )
    elif impact == "MODERATE":
        return (
            "🟠 This variant is predicted to have moderate functional impact. "
            "Periodic surveillance and stratified risk assessment based on personal and family history are advised."
        )
    else:
        return (
            "🟢 No strong evidence of clinical pathogenicity is currently associated with this variant. "
            "Routine follow-up is appropriate. Recommend re-evaluation if new evidence emerges or if patient has elevated familial risk."
        )

if variant:
    recommendation = get_recommendation(variant["Impact"], variant["Clinical Significance"])
    st.success(recommendation)
else:
    st.info("No recommendation available due to missing variant details.")
