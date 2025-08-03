import os
import pandas as pd
import streamlit as st
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# Enable wide layout
st.set_page_config(
    page_title="nf-xen: GenSheet",
    layout="wide",
    page_icon="fa-solid fa-file-csv"
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
    # Remove longest common suffix from sample names
    # Remove longest common suffix from sample names
    def longest_common_suffix(strs):
        if not strs:
            return ""
        rev_strs = [s[::-1] for s in strs]
        shortest = min(rev_strs, key=len)
        for i, ch in enumerate(shortest):
            for other in rev_strs:
                if other[i] != ch:
                    return shortest[:i][::-1]
        return shortest[::-1]

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
    # Bulk Edit Section
    st.write("### Bulk Edit Column Values")
    st.write("You can bulk edit a column by selecting it and entering a new value. This will update all rows in that column.")
    with st.form(key='bulk_edit_form', clear_on_submit=False):
        col_to_edit = st.selectbox("Select Column to Edit", df.columns)
        new_value = st.text_input("Enter New Value for All Rows", "")
        submit_button = st.form_submit_button(label='Apply to Column')

        if submit_button and col_to_edit and new_value != "":
            df[col_to_edit] = new_value
            st.success(f"Updated all values in '{col_to_edit}' column to '{new_value}'")
    # Configure AgGrid
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_default_column(editable=True, resizable=True, cellStyle={'textOverflow': 'ellipsis', 'whiteSpace': 'nowrap', 'overflow': 'hidden'})

    # Dropdowns for qc and clust
    gb.configure_column(
        "qc", 
        cellEditor='agSelectCellEditor', 
        cellEditorParams={'values': ["YES", "NO"]},
        editable=True
    )

    gb.configure_column(
        "clust", 
        cellEditor='agSelectCellEditor', 
        cellEditorParams={'values': ["YES", "NO"]},
        editable=True
    )

    # Limit width for certain columns (like h5ad) to trigger ellipsis
    gb.configure_column("h5ad", width=400)  # Adjust width as needed

    # Grid options
    gb.configure_grid_options(domLayout='normal')  # No autoHeight to keep single-line rows
    gb.configure_side_bar()

    grid_options = gb.build()

    grid_response = AgGrid(
        df,
        gridOptions=grid_options,
        editable=True,
        fit_columns_on_grid_load=False,  # Let user scroll horizontally
        theme='balham',
        update_mode=GridUpdateMode.MODEL_CHANGED,
        height=600,
        width='100%'
    )

    edited_df = grid_response['data']

    csv = edited_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download Sample Sheet CSV",
        data=csv,
        file_name='sample_sheet.csv',
        mime='text/csv'
    )
    # Additional Streamlit UI for workflow parameters
    st.sidebar.header("Workflow Parameters")

    # Profile Selection
    available_profiles = ["local", "docker", "arm", "apptainer", "conda", "uppmax", "test_full"]
    selected_profiles = st.sidebar.multiselect("Select Profiles", available_profiles, default=["docker"])

    # Optional Params Inputs
    cpus = st.sidebar.number_input("CPUs", min_value=1, value=14)
    memory = st.sidebar.text_input("Memory (e.g., 20.GB)", value="20.GB")
    retry = st.sidebar.number_input("Max Retries", min_value=0, value=3)
    time_param = st.sidebar.text_input("Max Time (e.g., 10-00:00:00)", value="10-00:00:00")
    account = st.sidebar.text_input("Account (Only for UPPMAX)", value="")

    # Clustering Method Dropdown (Louvain or Leiden)
    clust_method = st.sidebar.selectbox("Clustering Method", ["louvain", "leiden"], index=0)
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
    nextflow_cmd = f"nextflow run nf-xen -profile {profile_param} --sampleSheet {sample_sheet_path} --outdir {outdir} --cpus {cpus} --memory {memory} --retry {retry} --time {time_param} --clust_method {clust_method} --clust_res \"{clust_res}\""
    if account:
        nextflow_cmd += f" --account {account}"

    st.code(nextflow_cmd, language='bash')

else:
    st.info("Upload one or more .h5ad files to start.")
