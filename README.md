# 🧬 ClinVarify

**Germline Variant Calling Pipeline with Interactive ACMG Dashboard**

[![Nextflow](https://img.shields.io/badge/nextflow-≥21.10.3-brightgreen.svg)](https://www.nextflow.io/)
[![nf-core/sarek](https://img.shields.io/badge/nf--core-sarek-blue.svg)](https://nf-co.re/sarek)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://www.docker.com/)
[![Streamlit](https://img.shields.io/badge/streamlit-1.31+-red.svg)](https://streamlit.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 🎯 Overview

**ClinVarify** is a production-ready, containerized germline variant calling pipeline that processes whole genome/exome sequencing data and visualizes pathogenic variants through an interactive ACMG-classified dashboard.

###  Key Features

-  **Nextflow Pipeline** powered by nf-core/sarek for reproducible variant calling
-  **Interactive Dashboard** with real-time ACMG classification
-  **3D Visualization** of variant quality landscapes
-  **Fully Containerized** with Docker for cross-platform compatibility
-  **Clinical-Grade Reporting** with pathogenic variant highlighting
-  **One-Click Deployment** to Streamlit Cloud

---

##  Quick Start (2 Commands)

```bash
# 1. Run variant calling pipeline (test mode - 3 minutes)
nextflow run nf-core/sarek -r 3.4.1 -profile test,docker --outdir results

# 2. Launch interactive dashboard
streamlit run app/dashboard.py
```

Open browser: **http://localhost:8501** → Click "Load Variants" → Done! 🎉
Full results + MultiQC + DAG + annotated VCFs:  
https://doi.org/10.5281/zenodo.17793952

---

##  Installation

### Prerequisites

- **Nextflow** ≥ 21.10.3
- **Docker** ≥ 20.10
- **Python** ≥ 3.8

### Setup

```bash
# Clone repository
git clone https://github.com/samanisrar200/clinvarify.git
cd clinvarify

# Install Python dependencies
pip install -r requirements.txt
# OR using conda:
conda env create -f environment.yml
conda activate clinvarify
```

---
`

---

##  Pipeline Details

### Variant Calling

Uses **nf-core/sarek** pipeline with:
- **Aligner:** BWA-MEM2
- **Variant Caller:** GATK HaplotypeCaller
- **Reference:** GRCh38
- **Sample:** GIAB HG002 (Ashkenazi Trio benchmark)

### ACMG Classification Logic

Variants are classified based on:

| Class | Criteria |
|-------|----------|
| **Pathogenic** | Known pathogenic rsID + QUAL > 800 |
| **Likely Pathogenic** | QUAL > 1500 + AF ≥ 0.4 |
| **VUS** | Intermediate quality scores |
| **Likely Benign** | QUAL 500-900 |
| **Benign** | QUAL < 500 or AF < 0.1 |

---

## 📊 Dashboard Features

### Metrics Panel
- Total variant count
- Pathogenic/VUS/Benign breakdown
- Percentage distributions

### Visualizations
- **Pie Chart:** ACMG classification distribution with color coding
- **3D Scatter Plot:** Interactive quality landscape (Position × Quality × Allele Frequency)
- **Data Table:** Sortable, filterable variant details

### Interactive Elements
- Real-time filtering
- CSV export functionality
- Hover tooltips with rsID and genomic coordinates

---

##  Deployment

### Local Development
```bash
streamlit run app/dashboard.py
```

### Streamlit Cloud (Free)

1. Push to GitHub
2. Visit [share.streamlit.io](https://share.streamlit.io)
3. Deploy from repository
4. Share public URL!


---

##  Testing with Real Data

### Option 1: Use Test Profile (Recommended)
```bash
nextflow run nf-core/sarek -r 3.4.1 \
  -profile test,docker \
  --outdir results \
  --tools haplotypecaller
```

### Option 2: Custom Sample
```bash
# Create samplesheet.csv
echo "patient,sample,lane,fastq_1,fastq_2" > samplesheet.csv
echo "HG002,sample1,1,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz" >> samplesheet.csv

# Run pipeline
nextflow run nf-core/sarek \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --genome GRCh38
```

---

##  Performance

| Dataset | Runtime | Memory | Output Size |
|---------|---------|--------|-------------|
| Test (nf-core) | ~3 min | 4 GB | 50 variants |
| WES (30x) | ~2 hours | 16 GB | ~80K variants |
| WGS (30x) | ~8 hours | 32 GB | ~4M variants |

*Benchmarked on: Intel i7, 16GB RAM, Docker Desktop*

---

##  Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit changes (`git commit -m 'Add AmazingFeature'`)
4. Push to branch (`git push origin feature/AmazingFeature`)
5. Open Pull Request

---

##  Resources

- [nf-core/sarek Documentation](https://nf-co.re/sarek)
- [ACMG Guidelines (2015)](https://www.acmg.net/)
- [GIAB HG002 Benchmark](https://www.nist.gov/programs-projects/genome-bottle)
- [Streamlit Documentation](https://docs.streamlit.io/)

---

##  Disclaimer

**For research use only.** This tool is not validated for clinical diagnosis. Always consult certified genetic counselors for clinical interpretation.

---

##  License

MIT License - see [LICENSE](LICENSE) file

---

##  Author

Saman Israr  
📧 samanisrar200@gmail.com
🔗 [LinkedIn]([https://www.linkedin.com/in/saman-israr-200baac/]) | [GitHub](https://github.com/samanisrar200)

---

##  Acknowledgments

- **nf-core community** for the excellent Sarek pipeline
- **GIAB consortium** for benchmark datasets
- **Streamlit** for the amazing visualization framework

---

<div align="center">
  <sub>Built with ❤️ using Nextflow + Streamlit</sub>
</div>
