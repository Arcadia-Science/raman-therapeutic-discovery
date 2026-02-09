"""
Utility functions for parsing Raman spectroscopy data files.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import ramanspy as rp

ANNOTATION_FILE_DIR = "../peak-ID/"
POSSIBLE_WASATCH_COLUMNS = ['Wavelength', 'Wavenumber', 'Processed', 'Raw', 'Dark']

def collect_spectra_from_dir(data_dir, instrument, include_files=None, recursive=False) -> tuple:
    '''
    Collect data from files in the directory.

    Inputs:
        data_dir: str, directory with files
        instrument: str, type of instrument ('horiba', 'wasatch', 'openraman', 'instantraman')
        include_files: list of str, filenames to include. If none, all files in the directory will be included
        recursive: bool, if True, search for files recursively in subdirectories

    Returns:
        data_list: list of pd.DataFrame, list of dataframes with collected data
        sample_df: pd.DataFrame, sample table with columns 'ID', 'Name', 'Group ID'
    '''
    acceptable_ext = ['.txt', '.tsv', '.csv']
    if include_files:
        files = include_files
    else:
        if recursive:
            files = []
            data_path = Path(data_dir)
            for file_path in data_path.rglob('*'):
                if file_path.is_file() and file_path.suffix.lower() in acceptable_ext:
                    # Store relative path from data_dir
                    relative_path = file_path.relative_to(data_path)
                    files.append(str(relative_path))
        else:
            data_path = Path(data_dir)
            files = [f.name for f in data_path.iterdir() if f.is_file() and f.suffix.lower() in acceptable_ext]
    
    # remove sample_table.csv from files
    if 'sample_table.csv' in files:
        files.remove('sample_table.csv')
    
    # sort files in alphabetical order, first by parent directory, then by filename
    files.sort(key=lambda x: (Path(x).parent.as_posix(), Path(x).name))

    # Create sample_table.csv if it doesn't exist
    data_path = Path(data_dir)
    sample_table_path = data_path / 'sample_table.csv'
    if not sample_table_path.exists():
        print(f"Warning: Creating new sample_table.csv in {data_dir}")
        sample_data = []
        for i, file in enumerate(files):
            file_path = Path(file)
            sample_data.append({
                'ID': file,
                'Name': file_path.with_suffix(''),
                'Group ID': i
            })
        sample_df = pd.DataFrame(sample_data)
        sample_df.to_csv(sample_table_path, index=False)
    else:
        print(f"Loading existing sample_table.csv: {sample_table_path}")
        # Load existing sample table
        try:
            sample_df = pd.read_csv(sample_table_path, index_col=0, dtype={'ID': str, 'Name': str, 'Group ID': int})
        except Exception as e:
            raise ValueError(f"Error reading sample_table.csv: {e}. Ensure it has columns 'ID', 'Name', 'Group ID'.")
    data_list = []
    data_path = Path(data_dir)
    for file in files:
        file_path = data_path / file
        if instrument.lower() == 'horiba':
            df = parse_horiba_output(file_path)
        elif instrument.lower() == 'wasatch':
            df = parse_wasatch_output(file_path)
        elif instrument.lower() == 'openraman':
            df = parse_openraman_output(file_path)
        elif instrument.lower() == 'instantraman':
            df = parse_instantraman_output(file_path)
        else:
            raise ValueError('Unsupported instrument type: ', instrument)
        if df is not None:
            data_list.append((file, df))
    
    if not data_list:
        raise ValueError(f"No valid data files found in {data_dir}. Check the directory. Set `recursive=True` to search subdirectories.")
    print(f"Collected {len(data_list)} files from {data_dir}")
    return data_list, sample_df



def plot_spectra_in_dir(data_list,
                        data_dir: str | Path,
                        x_to_plot: str = "Wavenumber",
                        col_to_plot: str = "Intensity",
                        rp_pipeline=None,
                        include: list | None = None,
                        colors=None,
                        subtract_id: str | dict | None = None,
                        y_offset = 0,
                        **plotting_kwargs):
    '''
    Plot normalized data using RamanSPy for preprocessing.

    Inputs:
        data_list: list of pd.DataFrame, list of spectrum dataframes with normalized data
        col_to_plot: str, column in dataframe to plot
        data_dir: pathlib.Path, directory with files
        rp_pipeline: `ramanspy.preprocessing.Pipeline`, preprocessing pipeline to apply (optional)
        include: list (optional), sample names to include in the plot. Must exactly match the 'Name' column in `sample_table.csv`
        colors: list of str, colors to use for plotting. If None, uses default Plotly colors
        x_to_plot: str, column in dataframe to use as x-axis (default is "Wavenumber")
        subtract_id: str or dict, spectrum to subtract from other spectra (optional)
            - If str: filename or sample ID of spectrum to subtract from all other spectra
            - If dict: keys are filenames/sample IDs, values are filenames/sample IDs of spectra to subtract from each key
        y_offset: float, offset to apply to plot for clarity (optional). The default is 0
        plotting_kwargs: additional plotly keyword arguments for plotting

    Returns:
        fig: plotly Figure object
        spectra_names: list of str, names of processed spectra
        spectral_container: `ramanspy.SpectralContainer` containing all plotted spectra
    '''


    if colors is None:
        colors = px.colors.qualitative.D3

    fig = go.Figure()
    processed_spectra = []  # Store all processed spectra for SpectralContainer
    spectra_names = []  # Store names of processed spectra

    data_path = Path(data_dir)
    sample_table_path = data_path / 'sample_table.csv'
    if sample_table_path.exists():
        sample_table = pd.read_csv(sample_table_path, index_col=0, dtype={'ID': str, 'Name': str, 'Group ID': int})

        if include is not None:
            # filter sample_table to only include samples in include
            sample_table = sample_table[sample_table['Name'].isin(include)]

            # filter data_list to only include samples mentioned in sample_table
            data_list = [f for f in data_list if f[0] in sample_table.index.values]

        # Factorize 'Group ID' for color indexing
        sample_table["Group ID factorized"] = pd.factorize(sample_table['Group ID'])[0]

        # sort data_list based on sample_table
        sorted_data_list = []
        for index, _ in sample_table.iterrows():
            # Find the corresponding file in data_list
            matching_files = [f for f in data_list if f[0] == index]
            if len(matching_files) > 1:
                print(f"Warning: Multiple files match sample '{index}'. Using the first match.")
            elif len(matching_files) == 0:
                print(f"Warning: No files match sample '{index}'. Skipping.")
            if matching_files:
                sorted_data_list.append(matching_files[0])
            
   
    else:
        print('No sample table found!')
        sample_table = None

        if include is not None:
            # filter files to only include samples mentioned in include
            data_list = [f for f in data_list if Path(f[0]).stem in include]
        sorted_data_list = data_list

    # Process subtract spectrum if provided
    # Can be either a single string (subtract from all) or dict (specific subtractions)
    subtract_spectrum = None
    subtract_spectra_dict = {}

    if subtract_id is not None:
        # Helper function to find spectrum data by ID
        def find_spectrum_data(target_id):
            for file, df in sorted_data_list:
                if isinstance(sample_table, pd.DataFrame) and file in sample_table.index:
                    sample_name = sample_table.loc[file]['Name']
                    if target_id == sample_name:
                        return (file, df)
                if target_id == file:
                    return (file, df)
            return None

        if isinstance(subtract_id, str):
            # Original behavior: single subtract spectrum for all
            subtract_data = find_spectrum_data(subtract_id)

            if subtract_data is None:
                print(f"Warning: Could not find subtract spectrum '{subtract_id}'. Proceeding without subtraction.")
            else:
                # Process subtract spectrum
                file, df = subtract_data
                if col_to_plot not in df.columns:
                    raise ValueError(f"Column '{col_to_plot}' not found in subtract spectrum data.")
                if x_to_plot not in df.columns:
                    raise ValueError(f"Column '{x_to_plot}' not found in subtract spectrum data.")
                subtract_spectrum = rp.Spectrum(df[col_to_plot].values, df[x_to_plot].values)

        elif isinstance(subtract_id, dict):
            # New behavior: dict mapping sample IDs to their respective subtract spectra
            for key, value in subtract_id.items():
                subtract_data = find_spectrum_data(value)

                if subtract_data is None:
                    print(f"Warning: Could not find subtract spectrum '{value}' for key '{key}'. Skipping this subtraction.")
                else:
                    file, df = subtract_data
                    if col_to_plot not in df.columns:
                        raise ValueError(f"Column '{col_to_plot}' not found in subtract spectrum '{value}'.")
                    if x_to_plot not in df.columns:
                        raise ValueError(f"Column '{x_to_plot}' not found in subtract spectrum '{value}'.")
                    subtract_spectra_dict[key] = rp.Spectrum(df[col_to_plot].values, df[x_to_plot].values)

        else:
            raise TypeError("subtract_id must be either a str or dict")

    for i, (file, df) in enumerate(sorted_data_list):
        if col_to_plot not in df.columns:
            raise ValueError(f"Column '{col_to_plot}' not found in '{file}' spectrum data.")
        if x_to_plot not in df.columns:
            raise ValueError(f"Column '{x_to_plot}' not found in '{file}' spectrum data.")

        sample_id = Path(file).stem
        # Create RamanSPy Spectrum object
        spectrum = rp.Spectrum(df[col_to_plot].values, df[x_to_plot].values)

        # Get sample name if available
        if isinstance(sample_table, pd.DataFrame) and file in sample_table.index:
            sample_name = sample_table.loc[file]['Name']
        else:
            sample_name = sample_id

        # Subtract reference spectrum if provided (before applying pipeline)
        spectrum_to_subtract = None

        if subtract_spectrum is not None:
            # String case: subtract from all spectra (except the subtract spectrum itself)
            is_subtract_spectrum = (subtract_id == file or subtract_id == sample_id or subtract_id == sample_name)

            # Skip subtraction if this is the subtract spectrum itself
            if not is_subtract_spectrum:
                spectrum_to_subtract = subtract_spectrum

        elif subtract_spectra_dict:
            # Dict case: check if this spectrum has a specific subtract spectrum
            # Try to find the spectrum in the dict using file, sample_id, or sample_name as key
            for key_to_try in [file, sample_id, sample_name]:
                if key_to_try in subtract_spectra_dict:
                    spectrum_to_subtract = subtract_spectra_dict[key_to_try]
                    break

        # Perform subtraction if a spectrum was found
        if spectrum_to_subtract is not None:
            y_data_subtract = spectrum_to_subtract.spectral_data
            # Ensure spectral axes match for subtraction
            if len(spectrum.spectral_data) == len(y_data_subtract):
                spectrum.spectral_data = spectrum.spectral_data - y_data_subtract
        
        # Apply preprocessing pipeline if provided (after subtraction)
        if rp_pipeline is not None:
            spectrum = rp_pipeline.apply(spectrum)
            
        # Extract processed data for plotting
        x_data = spectrum.spectral_axis
        y_data = spectrum.spectral_data
        
        if isinstance(sample_table, pd.DataFrame):
            data_name = sample_table.loc[file]['Name']
            group_label = sample_table.loc[file]['Group ID factorized']
        else:
            data_name = sample_id
            group_label = i
            
        fig.add_trace(
            go.Scatter(
                x=x_data,
                y=y_data + y_offset * group_label,  # Offset each spectrum vertically for visibility
                name=data_name,
                line=dict(
                    color=colors[group_label]
                ),
                text=data_name.replace(' ', '<br>'),
                hoverinfo='text+x+y'
            )
        )
        
        # Add processed spectrum to collection
        processed_spectra.append(spectrum)
        spectra_names.append(data_name)

    
    fig.update_layout(
        showlegend=True,
        title_text=data_dir.parent.name,
        template='plotly_white',
        **plotting_kwargs
    )

    if len(processed_spectra) == 0:
        raise ValueError("No valid spectra found to plot. Check your data directory and filtering conditions.")
    # Create SpectralContainer with all processed spectra
    spectral_container = rp.SpectralContainer.from_stack(processed_spectra)

    return fig, spectra_names, spectral_container

def parse_wasatch_output(file_path):
    """
    Parse Wasatch spectrometer output file.
    
    The file format includes:
    - Metadata section with key,value pairs
    - Empty line separator
    - Data section with header: Wavelength,Wavenumber,Processed,Raw,Dark
    - Data rows with comma-separated values
    
    Parameters:
        file_path (str): Path to the Wasatch CSV file
        
    Returns:
        pd.DataFrame: DataFrame with columns ['Wavenumber', 'Intensity'] where
                     Intensity comes from the 'Processed' column
    """
    metadata = {}
    data_start = 0

    
    # Read the file line by line to extract metadata and identify where data starts
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # Parse metadata
    for i, line in enumerate(lines):
        if line.strip() == "":  # End of metadata block
            data_start = i + 1
            break
        if "," in line:
            key, _, value = line.partition(",")
            metadata[key.strip()] = value.strip()

    # Load the data portion into a DataFrame
    df = pd.read_csv(file_path, skiprows=data_start)

    # Keep only columns present in POSSIBLE_WASATCH_COLUMNS
    present_cols = [col for col in POSSIBLE_WASATCH_COLUMNS if col in df.columns]
    df = df[present_cols]

    # Rename Processed column to Intensity for compatibility
    df.rename(columns={'Processed': 'Intensity'}, inplace=True)

    # Remove columns that are all NaN (e.g. when "Dark" column is blank)
    df.dropna(axis='columns', how='all', inplace=True)
    # Remove rows with NaN values
    df.dropna(axis='rows', how="any", inplace=True)

    if not _is_like_spectrum(df):
        print(f"Warning: {file_path} does not resemble a spectrum. Expected at least two numeric columns. Ignoring this file.")
        return None
    return df


def parse_horiba_output(file_path):
    """
    Parse Horiba spectrometer output file.
    
    Files may have a header with metadata (denoted by #), followed by 2-column data (no column labels)
    
    Parameters:
        file_path (str): Path to the Horiba file
        
    Returns:
        pd.DataFrame: DataFrame with columns ['Wavenumber', 'Intensity']
    """
    def determine_sep(file_path):
        path = Path(file_path)
        ext = path.suffix.lower()
        if ext == '.csv':
            sep = ','
        elif ext == '.tsv':
            sep = '\t'
        elif ext == '.txt':
            sep = '\t'
        else:
            raise ValueError(f"Unsupported file extension: {ext}")
        return sep

    sep = determine_sep(file_path)
    data_lines = []
    metadata = {}
    
    with open(file_path, 'r', encoding='latin-1') as file:
        for line in file:
            if line.startswith("#"):
                # Split metadata by "=" and strip any extra spaces or tabs
                key_value = line[1:].split("=", 1)
                if len(key_value) == 2:
                    metadata[key_value[0].strip()] = key_value[1].strip()
            else:
                # Process non-metadata lines as data
                split_line = line.strip().split(sep)
                # if header line, skip
                if split_line[0][0].isalpha():
                    continue
                if len(split_line) == 2:
                    data_lines.append([float(split_line[0]), float(split_line[1])])

    # Convert the data into a DataFrame
    df = pd.DataFrame(data_lines, columns=["Wavenumber", "Intensity"])
    if not _is_like_spectrum(df):
        print(f"Warning: {file_path} does not resemble a spectrum. Expected at least two numeric columns. Ignoring this file.")
        return None
    return df


def parse_openraman_output(file_path):
    """
    Parse OpenRaman output file.
    
    Parameters:
        file_path (str): Path to the OpenRaman CSV file
        
    Returns:
        pd.DataFrame: DataFrame with columns ['Wavenumber', 'Intensity']
    """
    df = pd.read_csv(file_path)
    df.rename(columns={'Intensity (a.u.)': 'Intensity'}, inplace=True)
    if not _is_like_spectrum(df):
        print(f"Warning: {file_path} does not resemble a spectrum. Expected at least two numeric columns. Ignoring this file.")
        return None
    return df

def parse_instantraman_output(file_path):
    """
    Parse InstantRaman output file. Assuming no header and two columns: ['Pixel', 'Intensity'] but Pixel is actually Wavenumber.
    
    Parameters:
        file_path (str): Path to the InstantRaman CSV file
        
    Returns:
        pd.DataFrame: DataFrame with columns ['Pixel', 'Intensity']
    """
    df = pd.read_csv(file_path)
    if 'Wavenumber' not in df.columns and 'Pixel' in df.columns:
        # If Wavenumber column does not exist, rename 'Pixel' to 'Wavenumber'
        print("Warning: 'Pixel' column found, renaming to 'Wavenumber'.")
        df.rename(columns={'Pixel': 'Wavenumber'}, inplace=True)

    if not _is_like_spectrum(df):
        print(f"Warning: {file_path} does not resemble a spectrum. Expected at least two numeric columns. Ignoring this file.")
        return None
    return df

def _is_like_spectrum(df):
    """
    Heuristic to determine if a DataFrame resembles a spectrum.
    Checks to make sure there are at least 2 numeric columns.
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    return len(numeric_cols) >= 2
    

def parse_proteins_json(file_path):
    """
    Parse proteins.json file containing peak identification data.
    
    Parameters:
        file_path (str): Path to the proteins.json file
        
    Returns:
        dict: Dictionary with 'single' and 'range' keys containing peak data
    """
    with open(file_path, 'r') as f:
        data = json.load(f)
    return data


def add_peak_annotations(fig,
                         annotation_files = "proteins",
                         y_range=None):
    """
    Add peak annotations to a plotly figure based on peak annotation data.
    
    Parameters:
        fig (plotly.graph_objects.Figure): Plotly figure to annotate
        annotation_files (str or list): Name(s) of annotation file(s) to parse. 
                                       If string, will parse "../peak-ID/{name}.json"
                                       If list, will parse all files in the list
        y_range (tuple, optional): (y_min, y_max) for vertical lines and shaded areas
        
    Returns:
        plotly.graph_objects.Figure: Updated figure with peak annotations
    """

    # Handle single string or list of annotation files
    if isinstance(annotation_files, str):
        annotation_files = [annotation_files]
    
    # Collect all annotations from specified files
    all_protein_annotations = {}
    for annotation_file in annotation_files:
        annotation_dir = Path(ANNOTATION_FILE_DIR)
        file_path = annotation_dir / (annotation_file + ".json")
        try:
            file_annotations = parse_proteins_json(file_path)
            # Merge annotations, handling potential overlaps
            for protein_type, peaks in file_annotations.items():
                if protein_type in all_protein_annotations:
                    # If protein type already exists, merge the peaks
                    all_protein_annotations[protein_type].update(peaks)
                else:
                    all_protein_annotations[protein_type] = peaks
        except FileNotFoundError:
            print(f"Warning: Annotation file not found: {file_path}")
        except Exception as e:
            print(f"Warning: Error parsing annotation file {file_path}: {e}")
    
    if not all_protein_annotations:
        print("Warning: No valid annotation files found or loaded")
        return fig
    # Get x-axis range from existing data
    x_min = float('inf')
    x_max = float('-inf')
    for trace in fig.data:
        if hasattr(trace, 'x') and trace.x is not None:
            x_min = min(x_min, min(trace.x))
            x_max = max(x_max, max(trace.x))
    
    if y_range is None:
        # Get y-axis range from existing data
        y_min = float('inf')
        y_max = float('-inf')
        for trace in fig.data:
            if hasattr(trace, 'y') and trace.y is not None:
                y_min = min(y_min, min(trace.y))
                y_max = max(y_max, max(trace.y))
        y_range = (y_min, y_max)
    
    # Process each protein type and its peaks
    for protein_type, peaks in all_protein_annotations.items():
        for freq_key in peaks.keys():
            # Determine if this is a single frequency or range based on presence of "-"
            if '-' in freq_key:
                # Handle range - add shaded area
                start, end = freq_key.split('-')
                start_val = float(start)
                end_val = float(end)

                # Only add annotation if any part of the range overlaps with x-axis range
                if start_val <= x_max and end_val >= x_min:
                    # Clamp the range to the x-axis bounds
                    range_start = max(start_val, x_min)
                    range_end = min(end_val, x_max)

                    fig.add_vrect(
                        x0=range_start,
                        x1=range_end,
                        fillcolor="lightblue",
                        opacity=0.2,
                        layer="below",
                        line_width=0,
                        annotation_text=f"{protein_type}: {freq_key}",
                        annotation_position="top inside",
                        annotation_font_size=8
                    )
            else:
                # Handle single frequency - add vertical line
                freq_val = float(freq_key)
                # Only add annotation if frequency is within the x-axis range
                if x_min <= freq_val <= x_max:
                    fig.add_vline(
                        x=freq_val,
                        line=dict(color='red', width=1, dash='dash'),
                        annotation_text=f"{protein_type}: {freq_key}",
                        annotation_position="top",
                        annotation_font_size=8
                    )

    return fig
