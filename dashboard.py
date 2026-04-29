import streamlit as st
import pandas as pd
import gzip
import re
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go

# ─────────────────────────────────────────────
# Page config
# ─────────────────────────────────────────────
st.set_page_config(
    page_title="ClinVarify Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─────────────────────────────────────────────
# Custom CSS — dark clinical theme
# ─────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600&family=Syne:wght@400;700;800&display=swap');

html, body, [class*="css"] {
    font-family: 'Syne', sans-serif;
}
code, pre, .stDataFrame {
    font-family: 'JetBrains Mono', monospace !important;
}

/* Metric cards */
[data-testid="metric-container"] {
    background: linear-gradient(135deg, #0f1923 0%, #1a2942 100%);
    border: 1px solid #1e3a5f;
    border-radius: 12px;
    padding: 16px 20px;
}
[data-testid="metric-container"] label {
    color: #7eb8d4 !important;
    font-size: 0.75rem !important;
    text-transform: uppercase;
    letter-spacing: 0.1em;
}
[data-testid="metric-container"] [data-testid="metric-value"] {
    color: #e8f4f8 !important;
    font-size: 1.8rem !important;
    font-weight: 800;
}
[data-testid="metric-container"] [data-testid="metric-delta"] {
    font-size: 0.7rem !important;
}

/* Sidebar */
section[data-testid="stSidebar"] {
    background: #08111c;
    border-right: 1px solid #1e3a5f;
}
section[data-testid="stSidebar"] * {
    color: #b0cfe0 !important;
}

/* Tabs */
button[data-baseweb="tab"] {
    font-family: 'JetBrains Mono', monospace !important;
    font-size: 0.75rem;
    letter-spacing: 0.05em;
}

/* Download button */
div.stDownloadButton > button {
    background: linear-gradient(90deg, #0d6b9e, #0a8a7a) !important;
    color: white !important;
    border: none !important;
    border-radius: 8px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.8rem;
    letter-spacing: 0.05em;
}
div.stDownloadButton > button:hover {
    opacity: 0.88;
}

/* ACMG badge styling */
.acmg-badge {
    display: inline-block;
    padding: 2px 8px;
    border-radius: 4px;
    font-size: 0.7rem;
    font-weight: 700;
    font-family: 'JetBrains Mono', monospace;
    letter-spacing: 0.05em;
}
</style>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────
# Header
# ─────────────────────────────────────────────
st.markdown("""
<div style="
    background: linear-gradient(135deg, #07111d 0%, #0d2137 100%);
    border: 1px solid #1e3a5f;
    border-radius: 16px;
    padding: 28px 32px;
    margin-bottom: 24px;
">
    <h1 style="
        font-family: 'Syne', sans-serif;
        font-weight: 800;
        font-size: 2.2rem;
        color: #e8f4f8;
        margin: 0 0 4px 0;
        letter-spacing: -0.02em;
    ">🧬 ClinVarify</h1>
    <p style="
        color: #5a8fa8;
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.78rem;
        margin: 0;
        letter-spacing: 0.06em;
    ">GERMLINE VARIANT CALLING · ACMG CLASSIFICATION · INTERACTIVE DASHBOARD</p>
</div>
""", unsafe_allow_html=True)

# ─────────────────────────────────────────────
# ACMG colour palette
# ─────────────────────────────────────────────
ACMG_COLORS = {
    "Pathogenic":        "#e63946",
    "Likely Pathogenic": "#f4a261",
    "VUS":               "#e9c46a",
    "Likely Benign":     "#57cc99",
    "Benign":            "#48cae4",
}

# ─────────────────────────────────────────────
# Helper: parse AF from INFO field
# ─────────────────────────────────────────────
def parse_af(info: str) -> float | None:
    """Extract allele frequency from INFO string (AF=0.5 or AF=0.5,0.3)."""
    m = re.search(r'(?:^|;)AF=([^;]+)', info)
    if m:
        try:
            return float(m.group(1).split(',')[0])
        except ValueError:
            pass
    return None

# ─────────────────────────────────────────────
# Helper: classify variant per ACMG rules
# ─────────────────────────────────────────────
KNOWN_PATHOGENIC_IDS = {
    # A small representative set of known pathogenic rsIDs (extend as needed)
    "rs80357906", "rs28934578", "rs80357382", "rs786201005",
    "rs121912666", "rs121912667", "rs121912668",
}

def classify_acmg(row: pd.Series) -> str:
    qual = row["QUAL"]
    af   = row["AF"]
    rsid = str(row["ID"]).lower()

    # Pathogenic: known rsID AND quality support
    if rsid in KNOWN_PATHOGENIC_IDS and (pd.isna(qual) or qual > 800):
        return "Pathogenic"
    # Likely Pathogenic: very high quality + common allele
    if not pd.isna(qual) and qual > 1500 and not pd.isna(af) and af >= 0.4:
        return "Likely Pathogenic"
    # Benign: low quality OR very rare allele
    if (not pd.isna(qual) and qual < 500) or (not pd.isna(af) and af < 0.1):
        return "Benign"
    # Likely Benign: moderate quality
    if not pd.isna(qual) and 500 <= qual <= 900:
        return "Likely Benign"
    # Everything else → VUS
    return "VUS"

# ─────────────────────────────────────────────
# VCF loader (cached by path string)
# ─────────────────────────────────────────────
@st.cache_data(show_spinner="Parsing VCF …")
def load_vcf(vcf_path: str) -> pd.DataFrame:
    variants = []
    try:
        opener = gzip.open if vcf_path.endswith(".gz") else open
        mode   = "rt"
        with opener(vcf_path, mode) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                # Safe QUAL parse
                try:
                    qual = float(fields[5])
                except (ValueError, IndexError):
                    qual = None

                info = fields[7]
                af   = parse_af(info)
                is_annotated = bool(
                    re.search(r'(CLINVAR=|ANN=|CSQ=|CLNSIG=)', info)
                )

                variants.append({
                    "CHROM":      fields[0],
                    "POS":        int(fields[1]),
                    "ID":         fields[2],
                    "REF":        fields[3],
                    "ALT":        fields[4],
                    "QUAL":       qual,
                    "FILTER":     fields[6],
                    "AF":         af,
                    "ANNOTATED":  is_annotated,
                    "INFO":       info,
                })
    except Exception as e:
        st.error(f"Error loading VCF: {e}")
        return pd.DataFrame()

    df = pd.DataFrame(variants)
    if not df.empty:
        df["ACMG"] = df.apply(classify_acmg, axis=1)
    return df

# ─────────────────────────────────────────────
# Sidebar — file selection
# ─────────────────────────────────────────────
st.sidebar.markdown("## 📂 Select VCF File")

downloads_dir   = Path.home() / "Downloads"
pipeline_dir    = Path("results/variants/haplotypecaller")

vcf_files: list[Path] = []
for d in [downloads_dir, pipeline_dir]:
    if d.exists():
        vcf_files.extend(sorted(d.glob("*.vcf.gz")))
        vcf_files.extend(sorted(d.glob("*.vcf")))

selected_vcf: Path | None = None

if vcf_files:
    file_options = {f.name: f for f in vcf_files}
    selected_name = st.sidebar.selectbox("VCF File", list(file_options.keys()))
    selected_vcf  = file_options[selected_name]
    st.sidebar.info(f"📁 Source: `{selected_vcf.parent.name}/`")
else:
    st.sidebar.warning(
        "No `.vcf` or `.vcf.gz` found in:\n"
        "- `~/Downloads`\n"
        "- `results/variants/haplotypecaller`\n\n"
        "Download data from Zenodo first."
    )

st.sidebar.markdown("---")
st.sidebar.markdown("### 🎚️ Filters")
min_qual_filter = st.sidebar.slider("Minimum QUAL", 0, 5000, 0, 50)
acmg_filter     = st.sidebar.multiselect(
    "ACMG Class",
    options=list(ACMG_COLORS.keys()),
    default=list(ACMG_COLORS.keys()),
)
pass_only       = st.sidebar.checkbox("PASS variants only", value=False)

st.sidebar.markdown("---")
st.sidebar.markdown(
    "<small style='color:#3d6b82;'>ClinVarify · Research Use Only · MIT</small>",
    unsafe_allow_html=True,
)

# ─────────────────────────────────────────────
# Main content
# ─────────────────────────────────────────────
if selected_vcf is None:
    st.info("👈 Select a VCF file from the sidebar to begin.")
    st.stop()

df = load_vcf(str(selected_vcf))

if df.empty:
    st.warning("No variants parsed. Check VCF format.")
    st.stop()

# Apply sidebar filters
mask = (
    (df["QUAL"].fillna(0) >= min_qual_filter)
    & (df["ACMG"].isin(acmg_filter))
)
if pass_only:
    mask &= df["FILTER"].str.upper() == "PASS"

filtered_df = df[mask].copy()

st.success(
    f"✅ Loaded **{len(df):,}** variants from `{selected_vcf.name}` "
    f"({selected_vcf.parent.name}) · Showing **{len(filtered_df):,}** after filters."
)

# ─────────────────────────────────────────────
# Metrics row
# ─────────────────────────────────────────────
acmg_counts   = filtered_df["ACMG"].value_counts()
n_pass        = int((filtered_df["FILTER"].str.upper() == "PASS").sum())
n_annotated   = int(filtered_df["ANNOTATED"].sum())
n_pathogenic  = int(acmg_counts.get("Pathogenic", 0))
n_vus         = int(acmg_counts.get("VUS", 0))
n_benign      = int(
    acmg_counts.get("Benign", 0) + acmg_counts.get("Likely Benign", 0)
)
mean_qual     = filtered_df["QUAL"].mean()

c1, c2, c3, c4, c5, c6 = st.columns(6)
c1.metric("Total Variants",     f"{len(filtered_df):,}")
c2.metric("PASS",               f"{n_pass:,}",
          delta=f"{n_pass/max(len(filtered_df),1)*100:.0f}%")
c3.metric("🔴 Pathogenic",      f"{n_pathogenic:,}",
          delta=f"{n_pathogenic/max(len(filtered_df),1)*100:.1f}%")
c4.metric("🟡 VUS",             f"{n_vus:,}",
          delta=f"{n_vus/max(len(filtered_df),1)*100:.1f}%")
c5.metric("🟢 Benign + LB",     f"{n_benign:,}")
c6.metric("Mean QUAL",
          f"{mean_qual:.0f}" if pd.notna(mean_qual) else "N/A")

st.markdown("<br>", unsafe_allow_html=True)

# ─────────────────────────────────────────────
# Tabs
# ─────────────────────────────────────────────
tab_table, tab_pie, tab_3d, tab_hist, tab_chrom = st.tabs([
    "📋 Variants Table",
    "🥧 ACMG Distribution",
    "🌐 3-D Quality Landscape",
    "📈 QUAL Histogram",
    "🧬 Chromosomes",
])

# ── Tab 1: Data table ──────────────────────────────────────────
with tab_table:
    col_search, col_acmg_sel = st.columns([3, 2])
    with col_search:
        search_term = st.text_input(
            "🔍 Filter INFO / rsID (e.g. CLNSIG, GENE, rs80357906)",
            placeholder="Leave blank to show all",
        )
    with col_acmg_sel:
        table_acmg = st.multiselect(
            "ACMG class (table only)",
            options=list(ACMG_COLORS.keys()),
            default=list(ACMG_COLORS.keys()),
            key="table_acmg",
        )

    tdf = filtered_df[filtered_df["ACMG"].isin(table_acmg)].copy()
    if search_term:
        mask_search = (
            tdf["INFO"].str.contains(search_term, case=False, na=False)
            | tdf["ID"].str.contains(search_term, case=False, na=False)
        )
        tdf = tdf[mask_search]

    display_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "AF", "ACMG", "ANNOTATED"]
    st.dataframe(
        tdf[display_cols],
        use_container_width=True,
        height=480,
        column_config={
            "POS":  st.column_config.NumberColumn("Position", format="%d"),
            "QUAL": st.column_config.NumberColumn("QUAL", format="%.1f"),
            "AF":   st.column_config.NumberColumn("AF", format="%.3f"),
            "ACMG": st.column_config.TextColumn("ACMG Class"),
            "ANNOTATED": st.column_config.CheckboxColumn("Annotated"),
        },
    )

    st.download_button(
        "💾 Download filtered table as CSV",
        data=tdf[display_cols].to_csv(index=False).encode("utf-8"),
        file_name=f"{selected_vcf.stem}_clinvarify.csv",
        mime="text/csv",
    )

# ── Tab 2: ACMG pie chart ──────────────────────────────────────
with tab_pie:
    st.markdown("#### ACMG Classification Distribution")
    pie_data = (
        filtered_df["ACMG"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "ACMG", "ACMG": "Count"})
    )
    # Ensure column names are correct regardless of pandas version
    pie_data.columns = ["ACMG", "Count"]

    fig_pie = go.Figure(
        go.Pie(
            labels=pie_data["ACMG"],
            values=pie_data["Count"],
            marker_colors=[ACMG_COLORS.get(c, "#888") for c in pie_data["ACMG"]],
            hole=0.42,
            textinfo="label+percent",
            hovertemplate="<b>%{label}</b><br>Count: %{value:,}<br>%{percent}<extra></extra>",
        )
    )
    fig_pie.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font_color="#b0cfe0",
        legend=dict(font=dict(size=12)),
        margin=dict(t=20, b=20, l=20, r=20),
        height=420,
        annotations=[dict(
            text=f"<b>{len(filtered_df):,}</b><br>variants",
            x=0.5, y=0.5, font_size=15, showarrow=False,
            font_color="#e8f4f8",
        )],
    )
    st.plotly_chart(fig_pie, use_container_width=True)

# ── Tab 3: 3D quality landscape ───────────────────────────────
with tab_3d:
    st.markdown("#### 3-D Quality Landscape — Position × QUAL × Allele Frequency")

    plot_3d_df = filtered_df.dropna(subset=["QUAL", "AF"]).copy()

    if plot_3d_df.empty:
        st.info(
            "No variants have both QUAL and AF populated in INFO. "
            "Try a VCF annotated with GATK HaplotypeCaller (which writes AF to INFO)."
        )
    else:
        fig_3d = go.Figure(
            go.Scatter3d(
                x=plot_3d_df["POS"],
                y=plot_3d_df["QUAL"],
                z=plot_3d_df["AF"],
                mode="markers",
                marker=dict(
                    size=4,
                    color=[ACMG_COLORS.get(c, "#888") for c in plot_3d_df["ACMG"]],
                    opacity=0.82,
                    line=dict(width=0),
                ),
                customdata=plot_3d_df[["ID", "CHROM", "REF", "ALT", "ACMG"]].values,
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Chr: %{customdata[1]}  Pos: %{x:,}<br>"
                    "QUAL: %{y:.1f}  AF: %{z:.3f}<br>"
                    "REF: %{customdata[2]} → ALT: %{customdata[3]}<br>"
                    "ACMG: <b>%{customdata[4]}</b><extra></extra>"
                ),
            )
        )
        fig_3d.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            scene=dict(
                bgcolor="rgba(7,17,29,1)",
                xaxis=dict(title="Position", color="#5a8fa8", gridcolor="#1e3a5f"),
                yaxis=dict(title="QUAL Score", color="#5a8fa8", gridcolor="#1e3a5f"),
                zaxis=dict(title="Allele Freq (AF)", color="#5a8fa8", gridcolor="#1e3a5f"),
            ),
            font_color="#b0cfe0",
            margin=dict(t=10, b=10, l=10, r=10),
            height=560,
        )
        st.plotly_chart(fig_3d, use_container_width=True)

        # Colour legend
        leg_cols = st.columns(len(ACMG_COLORS))
        for col, (label, color) in zip(leg_cols, ACMG_COLORS.items()):
            n = int((filtered_df["ACMG"] == label).sum())
            col.markdown(
                f'<div style="text-align:center">'
                f'<span style="background:{color};border-radius:50%;'
                f'display:inline-block;width:12px;height:12px;margin-right:4px;"></span>'
                f'<span style="font-size:0.75rem;color:#b0cfe0;">{label}<br>'
                f'<b style="color:#e8f4f8;">{n:,}</b></span></div>',
                unsafe_allow_html=True,
            )

# ── Tab 4: QUAL histogram ─────────────────────────────────────
with tab_hist:
    st.markdown("#### Quality Score Distribution by ACMG Class")
    qual_df = filtered_df.dropna(subset=["QUAL"])
    if qual_df.empty:
        st.info("No QUAL values available.")
    else:
        fig_hist = px.histogram(
            qual_df,
            x="QUAL",
            color="ACMG",
            nbins=60,
            color_discrete_map=ACMG_COLORS,
            labels={"QUAL": "QUAL Score", "count": "Variant Count"},
            barmode="overlay",
            opacity=0.75,
        )
        fig_hist.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(7,17,29,1)",
            font_color="#b0cfe0",
            xaxis=dict(gridcolor="#1e3a5f"),
            yaxis=dict(gridcolor="#1e3a5f"),
            legend_title_text="ACMG",
            margin=dict(t=20, b=20),
            height=420,
        )
        st.plotly_chart(fig_hist, use_container_width=True)

# ── Tab 5: Chromosome bar chart ───────────────────────────────
with tab_chrom:
    st.markdown("#### Variant Count per Chromosome")
    chrom_df = (
        filtered_df.groupby(["CHROM", "ACMG"])
        .size()
        .reset_index(name="Count")
    )
    # Natural chromosome sort
    chrom_order = sorted(
        filtered_df["CHROM"].unique(),
        key=lambda c: (
            int(c.replace("chr", "")) if c.replace("chr", "").isdigit() else 100 + ord(c[-1])
        ),
    )
    fig_bar = px.bar(
        chrom_df,
        x="CHROM",
        y="Count",
        color="ACMG",
        color_discrete_map=ACMG_COLORS,
        category_orders={"CHROM": chrom_order},
        labels={"CHROM": "Chromosome", "Count": "Variant Count"},
        barmode="stack",
    )
    fig_bar.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(7,17,29,1)",
        font_color="#b0cfe0",
        xaxis=dict(gridcolor="#1e3a5f"),
        yaxis=dict(gridcolor="#1e3a5f"),
        legend_title_text="ACMG",
        margin=dict(t=20, b=20),
        height=420,
    )
    st.plotly_chart(fig_bar, use_container_width=True)
