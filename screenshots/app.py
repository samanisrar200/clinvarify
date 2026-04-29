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
st.markdown("Interactive viewer for HaplotypeCaller VCFs with ClinVar-ready annotations.")

# Sidebar for file selection
st.sidebar.header("Select VCF")
vcf_dir = Path("results/variants/haplotypecaller")
vcf_files = [f for f in vcf_dir.glob("*.vcf.gz") if f.is_file()]
selected_vcf = st.sidebar.selectbox("VCF File", vcf_files, format_func=lambda x: x.name)

if selected_vcf:
    # Load VCF
    @st.cache_data
    def load_vcf(vcf_path):
        variants = []
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

    df = load_vcf(selected_vcf)

    if not df.empty:
        st.success(f"Loaded {len(df):,} variants from {selected_vcf.name}.")

        # Metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Variants", len(df))
        with col2:
            st.metric("PASS Filter", len(df[df['FILTER'] == 'PASS']))
        with col3:
            st.metric("Avg QUAL", f"{df['QUAL'].mean():.1f}" if df['QUAL'].mean() else "N/A")
        with col4:
            st.metric("Unique Genes", df['INFO'].str.contains('GENEINFO=', na=False).sum())  # Placeholder

        # Tabs
        tab1, tab2, tab3 = st.tabs(["Variant Table", "Quality Distribution", "Chromosome View"])

        with tab1:
            st.subheader("Variants")
            # Search/filter
            search = st.text_input("Search by gene/symbol (INFO field)")
            filtered_df = df[df['INFO'].str.contains(search, case=False, na=False)] if search else df
            st.dataframe(filtered_df, use_container_width=True, height=600)

        with tab2:
            st.subheader("QUAL Scores")
            if df['QUAL'].notna().any():
                fig = px.histogram(df, x='QUAL', nbins=50, title="Variant Quality Distribution")
                st.plotly_chart(fig, use_container_width=True)

        with tab3:
            st.subheader("Variants by Chromosome")
            chrom_df = df['CHROM'].value_counts().reset_index()
            chrom_df.columns = ['CHROM', 'Count']
            fig = px.bar(chrom_df, x='CHROM', y='Count', title="Variants per Chromosome")
            st.plotly_chart(fig, use_container_width=True)

        # Download
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button("Download filtered CSV", csv, f"{selected_vcf.stem}_variants.csv", "text/csv")
    else:
        st.warning("No variants found in VCF.")
else:
    st.info("No VCF files detected. Place .vcf.gz files in results/variants/haplotypecaller/.")
