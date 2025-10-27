import os
import pandas as pd
import streamlit as st
from datetime import datetime

# Set the page icon to the nf-xen logo
st.set_page_config(
    page_title="nf-xen: GenSheet",
    layout="wide",
    page_icon="assets/nfxen.png"
)

# Default values
DEFAULTS = {
    'region': 'all_regions',
    'min_counts': 10,
    'max_counts': 500,
    'min_genes': 5,
    'max_genes': 110,
    'min_cell_area': 20,
    'max_cell_area': 500,
    'min_area_ratio': 0.1,
    'max_area_ratio': 0.9,
    'max_nucleus_area': 110,
    'min_cells_per_gene': 5,
    'qc': 'YES',
    'clust': 'YES',
    'clust_res': 'NA',
    'celldex_ref': 'hpca__2024-02-26',
    'celldex_labs': 'label.main'
}

st.title("nf-xen sample sheet generator")

uploaded_files = st.file_uploader(
    "Upload .h5ad files", 
    type=['h5ad'], 
    accept_multiple_files=True
)

if uploaded_files:
    data = []
    for file in uploaded_files:
        full_path = os.path.abspath(file.name)
        sample_name = os.path.splitext(os.path.basename(file.name))[0]
        entry = {
            'sampleName': sample_name,
            'h5ad': full_path
        }
        entry.update(DEFAULTS)
        data.append(entry)

    df = pd.DataFrame(data)
    # Remove longest common suffix from sample names (only if >1 sample)
    def longest_common_suffix(strs):
        if not strs or len(strs) < 2:
            return ""
        rev_strs = [s[::-1] for s in strs]
        shortest = min(rev_strs, key=len)
        for i, ch in enumerate(shortest):
            for other in rev_strs:
                if other[i] != ch:
                    return shortest[:i][::-1]
        return shortest[::-1]

    if len(df['sampleName']) > 1:
        suffix = longest_common_suffix(df['sampleName'].tolist())
        if suffix:
            df['sampleName'] = df['sampleName'].str[:-len(suffix)]

    st.write("### Edit Sample Sheet Below:")
    st.write("You can edit the values below. The columns `qc` and `clust` are dropdowns with options `YES` and `NO`")
    st.write("The `h5ad` column is the path to the uploaded file. You can edit it if needed, but it should point to the correct file path.")
    st.write("The sample names are derived from the uploaded file names, and you can edit them as needed.")
    st.write("The longest common suffix has been removed from the sample names to make them more readable.")
    st.write("Don't forget to scroll horizontally to see all columns as they may not be visible.")
    st.write("Github: [nf-xen](https://github.com/addityea/nf-xen): Aditya Singh")
    st.write("#### Editor")
    st.write("Functions similar to a spreadsheet editor such as Excel or Google Sheets.")
    st.write("For bulk editing, edit one cell and drag the fill handle (small square at the bottom-right corner of the cell) to fill other cells with the same value.")
    edited_df = st.data_editor(
        df,
        column_config={
            "qc": st.column_config.SelectboxColumn("QC", options=["YES", "NO"], default="YES", required=True),
            "clust": st.column_config.SelectboxColumn("Clust", options=["YES", "NO"], default="YES", required=True)
        },
        hide_index=True
    )

    csv = edited_df.to_csv(index=False).encode('utf-8')
    now_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    st.download_button(
        label="📥 Download Sample Sheet CSV",
        data=csv,
        file_name=f'nf-gen_sample_sheet_{now_str}.csv',
        mime='text/csv'
    )
    # Additional Streamlit UI for workflow parameters
    st.sidebar.image("../assets/nfxen.png", width=200)
    # Add a line
    st.sidebar.markdown("---")
    st.sidebar.header("Workflow Parameters")

    # Profile Selection
    available_profiles = ["local", "docker", "arm", "apptainer", "conda", "uppmax", "test_full", "offline"]
    selected_profiles = st.sidebar.multiselect("Select Profiles", available_profiles, default=["docker"])

    # Optional Params Inputs
    cpus_col, retry_col = st.sidebar.columns(2)
    cpus = cpus_col.number_input("CPUs", min_value=1, value=(os.cpu_count() -1), max_value = os.cpu_count())
    retry = retry_col.number_input("Max Retries", min_value=0, value=3)
    mcol1, mcol2 = st.sidebar.columns(2)
    mem = mcol1.number_input("Memory", min_value=1, value=20, max_value=1000)
    mem_unit = mcol2.selectbox("Memory Unit", ["GB", "MB", "TB"], index=0)
    memory = f"{mem}.{mem_unit}" if mem_unit != "GB" else f"{mem}.GB"
    #memory = st.sidebar.text_input("Memory (e.g., 20.GB)", value="20.GB")
    st.sidebar.write("Max Time")
    col1, col2, col3, col4 = st.sidebar.columns(4)
    days = col1.number_input("Days", min_value=0, value=10)
    hours = col2.number_input("Hours", min_value=0, max_value=24, value=0)
    mins = col3.number_input("Minutes", min_value=0, max_value=60, value=0)
    secs = col4.number_input("Seconds", min_value=0, max_value=60, value=0)
    time_param = ""
    if days == 0 and hours == 0 and mins == 0 and secs == 0:
        st.sidebar.error("At least one time unit must be greater than 0.")
        st.stop()
    if days == 0:
        time_param = f"{hours:02d}:{mins:02d}:{secs:02d}"
    else:
        time_param = f"{days}-{hours:02d}:{mins:02d}:{secs:02d}"
    account = st.sidebar.text_input("Account (Only for UPPMAX)", value="")

    # Clustering Method Dropdown (Louvain or Leiden)
    clust_method = st.sidebar.selectbox("Clustering Method", ["leiden", "louvain"], index=0)
    clust_res = st.sidebar.text_input("Clustering Resolution(s)", value="0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5")

    # Output Directory
    outdir = st.sidebar.text_input("Output Directory", value="results")

    # Generate the --profile parameter (comma-separated)
    profile_param = ",".join(selected_profiles)

    # Sample Sheet Path
    sample_sheet_path = 'sample_sheet.csv'

    # Nextflow Command Preview
    st.write("### Generated Nextflow Run Command")
    st.write("Make sure to run this command from within the `nf-xen` directory.")
    st.write("Update the `--sampleSheet` parameter if it's not in the nf-xen directory.")
    nextflow_cmd = f"nextflow run main.nf -profile {profile_param} --sampleSheet {sample_sheet_path} --outdir {outdir} --cpus {cpus} --memory {memory} --retry {retry} --time {time_param} --clust_method {clust_method} --clust_res \"{clust_res}\""
    if account:
        nextflow_cmd += f" --account {account}"

    st.code(nextflow_cmd, language='bash')

else:
    st.info("Upload one or more .h5ad files to start.")
st.markdown("---")
# Add a footer with emojis and art
st.markdown("""
    <style>
        footer {
            text-align: center;
            padding: 10px;
            background-color: #f1f1f1;
            border-top: 1px solid #e0e0e0;
        }
        footer p {
            margin: 0;
        }
        .footer-icons {
            font-size: 1.5em;
            margin-top: 8px;
        }
        .footer-art {
            font-family: monospace;
            font-size: 0.95em;
            color: #888;
            margin-top: 10px;
            white-space: pre;
        }
    </style>
    <footer>
        <p>Made by Aditya Singh</p>
        <p>Follow me on GitHub: <a href="https://www.github.com/addityea" target="_blank">addityea</a></p>
        <p>Powered by Python 🐍, Streamlit 🎈 and Nextflow 🧬</p>
        <div class="footer-icons">🚀✨🧑‍💻📊🧬</div>
    </footer>
""", unsafe_allow_html=True)