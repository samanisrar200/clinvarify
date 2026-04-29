import streamlit as st
import pandas as pd
import gzip
import os
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

# Page config
st.set_page_config(page_title="ClinVarify Dashboard", layout="wide")

st.title("🧬 ClinVarify Variant Dashboard")
st.markdown("Interactive viewer for VCFs (Zenodo Downloads or Pipeline Outputs).")

# Sidebar for file selection - prioritize Downloads
st.sidebar.header("Select VCF")
downloads_dir = Path.home() / "Downloads"
pipeline_dir = Path("results/variants/haplotypecaller")

# Collect files from both locations
vcf_files = []
if downloads_dir.exists():
    vcf_files.extend(list(downloads_dir.glob("*.vcf.gz")))
if pipeline_dir.exists():
    vcf_files.extend(list(pipeline_dir.glob("*.vcf.gz")))

if vcf_files:
    # Display with paths
    file_options = {f.name: f for f in vcf_files}
    selected_name = st.sidebar.selectbox("VCF File", list(file_options.keys()))
    selected_vcf = file_options[selected_name]
    
    st.sidebar.info(f"Loading from: {selected_vcf.parent.name}")
else:
    st.sidebar.warning("No .vcf.gz in ~/Downloads or results/. Download from Zenodo first.")

if 'selected_vcf' in locals() and selected_vcf:
    # Load VCF (cached)
    @st.cache_data
    def load_vcf(vcf_path):
        variants = []
        try:
            with gzip.open(vcf_path, 'rt') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 8:
                        variants.append({
                            'CHROM': fields[0],
                            'POS': int(fields[1]),
                            'ID': fields[2],
                            'REF': fields[3],
                            'ALT': fields[4],
                            'QUAL': float(fields[5]) if fields[5] != '.' else None,
                            'FILTER': fields[6],
                            'INFO': fields[7]
                        })
            return pd.DataFrame(variants)
        except Exception as e:
            st.error(f"Error loading {vcf_path.name}: {e}")
            return pd.DataFrame()

    df = load_vcf(selected_vcf)

    if not df.empty:
        st.success(f"✅ Loaded {len(df):,} variants from **{selected_vcf.name}** ({selected_vcf.parent.name}).")

        # Metrics row
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Variants", len(df))
        with col2:
            st.metric("PASS Variants", len(df[df['FILTER'] == 'PASS']))
        with col3:
            st.metric("Mean QUAL", f"{df['QUAL'].mean():.1f}" if pd.notna(df['QUAL']).any() else "N/A")
        with col4:
            total_info = len(df['INFO'].str.contains('CLINVAR=|ANN=', na=False))  # ClinVar/SnpEff placeholder
            st.metric("Annotated", total_info)

        # Tabs for views
        tab1, tab2, tab3 = st.tabs(["📊 Table", "📈 QUAL Hist", "🧬 Chrom Bars"])

        with tab1:
            search_gene = st.text_input("Filter INFO (e.g., CLNSIG, GENE)")
            filtered_df = df[df['INFO'].str.contains(search_gene, case=False, na=False)] if search_gene else df
            st.dataframe(filtered_df, use_container_width=True, height=500,
                         column_config={"POS": st.column_config.NumberColumn("Position", format="%.0f")})

        with tab2:
            if pd.notna(df['QUAL']).any():
                fig_hist = px.histogram(df.dropna(subset=['QUAL']), x='QUAL', nbins=50,
                                        title="Quality Score Distribution",
                                        labels={'QUAL': 'QUAL Score'})
                st.plotly_chart(fig_hist, use_container_width=True)

        with tab3:
            chrom_counts = df['CHROM'].value_counts().reset_index()
            chrom_counts.columns = ['CHROM', 'Count']
            fig_bar = px.bar(chrom_counts.head(20), x='CHROM', y='Count',
                             title="Top Chromosomes by Variant Count")
            st.plotly_chart(fig_bar, use_container_width=True)

        # Export
        csv_data = filtered_df.to_csv(index=False).encode('utf-8')
        st.download_button("💾 Download CSV", csv_data, f"{selected_vcf.stem}_clinvarify.csv",
                           "text/csv")
    else:
        st.warning("No variants parsed. Check VCF format.")
