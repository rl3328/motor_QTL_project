#!/usr/bin/env python3
# Usage:
# # Combine all TWAS files in a directory
# python combine_results.py /path/to/your/folder --type twas

# # Combine all MR files in a directory
# python combine_results.py /path/to/your/folder --type mr

# # Combine both TWAS and MR files (auto-detect)
# python combine_results.py /path/to/your/folder --type both

# # Combine with custom output name
# python combine_results.py /path/to/your/folder --type twas --output my_combined_results

# # List all available file patterns
# python combine_results.py /path/to/your/folder --list-patterns

"""
Script to combine TWAS and MR (Mendelian Randomization) result files
Handles *.twas.tsv.gz and *.mr_result.tsv.gz files with proper error handling and validation
Skips empty files but preserves all logging messages
"""

import os
import glob
import gzip
import pandas as pd
import argparse
import re
import sys
import logging
from datetime import datetime
from pathlib import Path
from collections import Counter

def setup_logging(directory, output_name=None, file_type='results'):
    """Setup logging to both console and file"""
    # Create log filename
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if output_name:
        # Extract just the basename if a full path was provided
        base_name = os.path.basename(output_name)
        # Remove any extensions
        base_name = base_name.replace('.tsv.gz', '').replace('.gz', '').replace('.tsv', '')
        log_filename = f"{base_name}_combine_{file_type}_{timestamp}.log"
    else:
        log_filename = f"combine_{file_type}_{timestamp}.log"
    
    log_path = os.path.join(directory, log_filename)
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logger = logging.getLogger(__name__)
    logger.info(f"=== File Combination Log ===")
    logger.info(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Log file: {log_filename}")
    logger.info(f"Working directory: {directory}")
    
    return logger, log_path

def extract_region_from_filename(filename):
    """Extract region information from filename (format: chrN_start_end)"""
    # Pattern to match chr followed by chromosome identifier and coordinates
    pattern = r'(chr\w+_\d+_\d+)'
    match = re.search(pattern, filename)
    if match:
        return match.group(1)
    
    # Alternative pattern for just chromosome
    pattern2 = r'(chr\w+)'
    match2 = re.search(pattern2, filename)
    if match2:
        return match2.group(1)
    
    return None

def extract_study_context_from_filename(filename, file_type='twas'):
    """Extract study and context information from filename"""
    # Remove the appropriate extension
    if file_type == 'twas':
        base = filename.replace('.twas.tsv.gz', '')
    elif file_type == 'mr':
        base = filename.replace('.mr_result.tsv.gz', '')
    else:
        base = filename
    
    # Split by dots and underscores to find components
    parts = re.split(r'[._]', base)
    
    study_context = []
    for part in parts:
        if part not in ['chr', 'twas', 'mr', 'result', 'tsv', 'gz'] and not re.match(r'^\d+$', part):
            study_context.append(part)
    
    return '_'.join(study_context[:3]) if study_context else 'unknown'

def find_result_files(directory, file_type='twas'):
    """Find all result files based on type"""
    if file_type == 'twas':
        pattern = "*.twas.tsv.gz"
    elif file_type == 'mr':
        pattern = "*.mr_result.tsv.gz"
    elif file_type == 'both':
        twas_files = glob.glob(os.path.join(directory, "*.twas.tsv.gz"))
        mr_files = glob.glob(os.path.join(directory, "*.mr_result.tsv.gz"))
        return sorted(twas_files), sorted(mr_files)
    else:
        raise ValueError(f"Unknown file type: {file_type}")
    
    files = glob.glob(os.path.join(directory, pattern))
    return sorted(files)

def is_file_empty(filepath, logger=None):
    """Check if a gzipped file is empty (header only or completely empty)"""
    try:
        with gzip.open(filepath, 'rt') as f:
            # Check if file has any content
            header = f.readline().strip()
            if not header:
                return True, "File is completely empty"
            
            # Check if file has only header (no data rows)
            data_line = f.readline().strip()
            if not data_line:
                return True, "File contains only header, no data rows"
            
            return False, None
            
    except Exception as e:
        error_msg = f"Error checking if file is empty: {str(e)}"
        if logger:
            logger.error(error_msg)
        return True, error_msg

def validate_file_structure(filepath, file_type='twas', expected_columns=None, logger=None):
    """Validate file structure and return column info"""
    try:
        # First check if file is empty
        is_empty, empty_reason = is_file_empty(filepath, logger)
        if is_empty:
            return False, f"Empty file: {empty_reason}"
        
        with gzip.open(filepath, 'rt') as f:
            # Read first few lines to check structure
            header = f.readline().strip()
            if not header:
                return False, "Empty file or no header"
            
            columns = header.split('\t')
            
            # Check for essential columns based on file type
            if file_type == 'twas':
                essential_cols = ['chr', 'molecular_id', 'twas_z', 'twas_pval']
            elif file_type == 'mr':
                essential_cols = ['gene_name', 'context', 'gwas_study']
            else:
                essential_cols = []
            
            missing_essential = [col for col in essential_cols if col not in columns]
            
            if missing_essential:
                error_msg = f"Missing essential columns: {missing_essential}"
                if logger:
                    logger.warning(f"File validation failed: {error_msg}")
                return False, error_msg
            
            # Try to read one data line to ensure file is readable
            try:
                data_line = f.readline().strip()
                if data_line:
                    data_parts = data_line.split('\t')
                    if len(data_parts) != len(columns):
                        error_msg = f"Column count mismatch: header has {len(columns)}, data has {len(data_parts)}"
                        if logger:
                            logger.warning(f"File validation failed: {error_msg}")
                        return False, error_msg
            except:
                pass  # File might be header-only, which we now catch as empty
            
            return True, columns
            
    except Exception as e:
        error_msg = f"Error reading file: {str(e)}"
        if logger:
            logger.error(f"File validation error: {error_msg}")
        return False, error_msg

def read_result_file(filepath, file_type='twas', logger=None):
    """Read a result file with comprehensive error handling"""
    filename = os.path.basename(filepath)
    
    # First check if file is empty
    is_empty, empty_reason = is_file_empty(filepath, logger)
    if is_empty:
        skip_msg = f"SKIPPED (empty): {empty_reason}"
        if logger:
            logger.info(f"{filename} - {skip_msg}")
        print(f"    {skip_msg}")
        return None, f"Empty file: {empty_reason}"
    
    # Then validate the file structure
    is_valid, result = validate_file_structure(filepath, file_type=file_type, logger=logger)
    if not is_valid:
        if logger:
            logger.warning(f"{filename} - {result}")
        print(f"    WARNING: {filename} - {result}")
        return None, result
    
    try:
        with gzip.open(filepath, 'rt') as f:
            # Read with comprehensive missing value recognition
            df = pd.read_csv(f, sep='\t', 
                           na_values=['', ' ', 'NA', 'na', 'N/A', 'n/a', 'NaN', 'nan', 
                                     'NULL', 'null', 'None', 'none', '.', '-', '--', '---',
                                     'missing', 'Missing', 'MISSING', 'Inf', '-Inf'],
                           keep_default_na=True,
                           dtype=str)  # Read as strings first to handle mixed types
        
        # Convert numeric columns based on file type
        if file_type == 'twas':
            numeric_cols = ['TSS', 'start', 'end', 'rsq_cv', 'pval_cv', 'twas_z', 'twas_pval']
        elif file_type == 'mr':
            numeric_cols = ['nsnp', 'b', 'se', 'pval', 'lo_ci', 'up_ci', 'or', 'or_lci95', 'or_uci95']
        else:
            numeric_cols = []
        
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Basic data validation
        if len(df) == 0:
            skip_msg = "SKIPPED (empty): No data rows after parsing"
            if logger:
                logger.info(f"{filename} - {skip_msg}")
            print(f"    {skip_msg}")
            return None, "No data rows after parsing"
        
        # Check for completely empty rows
        empty_rows = df.isnull().all(axis=1).sum()
        if empty_rows > 0:
            msg = f"Removing {empty_rows} completely empty rows"
            if logger:
                logger.info(f"{filename} - {msg}")
            print(f"    INFO: {msg}")
            df = df.dropna(how='all')
            
            # Check if all rows were empty
            if len(df) == 0:
                skip_msg = "SKIPPED (empty): All rows were empty after cleanup"
                if logger:
                    logger.info(f"{filename} - {skip_msg}")
                print(f"    {skip_msg}")
                return None, "All rows were empty after cleanup"
        
        return df, None
        
    except pd.errors.EmptyDataError:
        error_msg = "SKIPPED (empty): EmptyDataError - no data to parse"
        if logger:
            logger.info(f"{filename} - {error_msg}")
        print(f"    {error_msg}")
        return None, error_msg
    except pd.errors.ParserError as e:
        error_msg = f"Parser error: {str(e)}"
        if logger:
            logger.error(f"{filename} - {error_msg}")
        return None, error_msg
    except Exception as e:
        error_msg = f"Unexpected error: {str(e)}"
        if logger:
            logger.error(f"{filename} - {error_msg}")
        return None, error_msg

def get_all_columns(files, file_type='twas', max_files_to_check=10, logger=None):
    """Scan files to get complete column list, skipping empty files"""
    msg = f"Scanning up to {max_files_to_check} files for column information..."
    if logger:
        logger.info(msg)
    print(msg)
    
    all_columns = set()
    column_frequency = Counter()
    
    files_checked = 0
    files_skipped_empty = 0
    
    for filepath in files[:max_files_to_check]:
        filename = os.path.basename(filepath)
        
        # Check if file is empty first
        is_empty, empty_reason = is_file_empty(filepath, logger)
        if is_empty:
            files_skipped_empty += 1
            if logger:
                logger.info(f"Skipping empty file during column scan: {filename} - {empty_reason}")
            continue
            
        try:
            with gzip.open(filepath, 'rt') as f:
                header = f.readline().strip().split('\t')
                all_columns.update(header)
                column_frequency.update(header)
                files_checked += 1
        except Exception as e:
            warning_msg = f"Could not read header from {filename}: {e}"
            if logger:
                logger.warning(warning_msg)
            print(f"    Warning: {warning_msg}")
    
    result_msg = f"Checked {files_checked} files, skipped {files_skipped_empty} empty files, found {len(all_columns)} unique columns"
    if logger:
        logger.info(result_msg)
    print(f"    {result_msg}")
    
    # Show column frequency for debugging
    if len(all_columns) > 20:
        columns_msg = "Most common columns:"
        if logger:
            logger.info(columns_msg)
        print(f"    {columns_msg}")
        for col, freq in column_frequency.most_common(10):
            col_msg = f"{col}: {freq}/{files_checked} files"
            if logger:
                logger.info(f"      {col_msg}")
            print(f"      {col_msg}")
    
    return sorted(all_columns)

def get_standard_column_order(file_type='twas'):
    """Define standard column order based on file type"""
    if file_type == 'twas':
        return [
            'chr', 'molecular_id', 'TSS', 'start', 'end', 'context', 
            'gwas_study', 'method', 'is_imputable', 'is_selected_method',
            'rsq_cv', 'pval_cv', 'twas_z', 'twas_pval', 'type', 'block'
        ]
    elif file_type == 'mr':
        return [
            'exposure', 'outcome', 'method', 'nsnp', 'b', 'se', 'pval',
            'lo_ci', 'up_ci', 'or', 'or_lci95', 'or_uci95',
            'Q', 'Q_df', 'Q_pval', 'egger_intercept', 'se_intercept', 'pval_intercept'
        ]
    else:
        return []

def add_missing_columns(df, all_columns):
    """Add any missing columns with NA values"""
    for col in all_columns:
        if col not in df.columns:
            df[col] = pd.NA
    return df

def reorder_columns(df, standard_order):
    """Reorder columns with standard order first, then metadata"""
    final_order = []
    
    # Add standard columns first (if they exist)
    for col in standard_order:
        if col in df.columns:
            final_order.append(col)
    
    # Add any extra columns (excluding metadata)
    for col in sorted(df.columns):
        if col not in standard_order and col not in ['region', 'source_file', 'study_context']:
            final_order.append(col)
    
    # Add metadata columns at the end
    for meta_col in ['region', 'study_context', 'source_file']:
        if meta_col in df.columns:
            final_order.append(meta_col)
    
    return df[final_order]

def combine_result_files(directory, file_type='twas', output_name=None):
    """
    Combine all result files in a directory with comprehensive error handling
    Skips empty files but preserves all logging messages
    
    Args:
        directory: Path to directory containing result files
        file_type: Type of files to combine ('twas', 'mr', or 'both')
        output_name: Optional custom output filename
    """
    
    # Setup logging
    logger, log_path = setup_logging(directory, output_name, file_type)
    
    logger.info(f"=== Combining {file_type.upper()} files from {directory} ===")
    print(f"=== Combining {file_type.upper()} files from {directory} ===")
    
    # Find all result files
    if file_type == 'both':
        twas_files, mr_files = find_result_files(directory, file_type)
        files_to_process = [
            ('twas', twas_files),
            ('mr', mr_files)
        ]
    else:
        files = find_result_files(directory, file_type)
        files_to_process = [(file_type, files)]
    
    # Process each file type
    for current_type, files in files_to_process:
        if not files:
            warning_msg = f"No {current_type.upper()} files found in directory"
            logger.warning(warning_msg)
            print(warning_msg)
            continue
        
        process_file_type(directory, current_type, files, output_name, logger, log_path)

def process_file_type(directory, file_type, files, output_name, logger, log_path):
    """Process files of a specific type"""
    
    info_msg = f"\n=== Processing {file_type.upper()} files ==="
    logger.info(info_msg)
    print(info_msg)
    
    info_msg = f"Found {len(files)} {file_type.upper()} files to process"
    logger.info(info_msg)
    print(info_msg)
    
    # Get all possible columns (this will now skip empty files)
    all_columns = get_all_columns(files, file_type=file_type, logger=logger)
    standard_order = get_standard_column_order(file_type)
    
    # Track processing statistics
    combined_data = []
    regions_processed = []
    study_contexts = set()
    total_results = 0
    failed_files = []
    empty_files = []
    
    for i, filepath in enumerate(files, 1):
        filename = os.path.basename(filepath)
        region = extract_region_from_filename(filename)
        study_context = extract_study_context_from_filename(filename, file_type)
        
        progress_msg = f"[{i}/{len(files)}] Processing: {filename}"
        logger.info(progress_msg)
        print(progress_msg)
        
        if region:
            regions_processed.append(region)
        study_contexts.add(study_context)
        
        # Read the file
        df, error = read_result_file(filepath, file_type=file_type, logger=logger)
        
        if df is not None:
            # Add missing columns to ensure consistency
            df = add_missing_columns(df, all_columns)
            
            # Add metadata
            if region:
                df['region'] = region
            df['study_context'] = study_context
            df['source_file'] = filename
            
            # Reorder columns
            df = reorder_columns(df, standard_order)
            
            combined_data.append(df)
            result_count = len(df)
            total_results += result_count
            
            success_msg = f"SUCCESS: {result_count:,} results"
            logger.info(f"    {success_msg}")
            print(f"    {success_msg}")
            
        else:
            # Check if it was empty vs other error
            if "empty" in error.lower():
                empty_files.append((filename, error))
                empty_msg = f"EMPTY FILE: {error}"
                logger.info(f"    {empty_msg}")
                print(f"    {empty_msg}")
            else:
                failed_files.append((filename, error))
                fail_msg = f"FAILED: {error}"
                logger.error(f"    {fail_msg}")
                print(f"    {fail_msg}")
    
    # Report empty files separately from failures
    if empty_files:
        empty_header = f"=== EMPTY FILES SKIPPED ({len(empty_files)} files) ==="
        logger.info(empty_header)
        print(f"\n{empty_header}")
        for filename, error in empty_files:
            empty_detail = f"{filename}: {error}"
            logger.info(f"  {empty_detail}")
            print(f"  {empty_detail}")
    
    # Report actual failures
    if failed_files:
        failure_header = f"=== PROCESSING FAILURES ({len(failed_files)} files) ==="
        logger.warning(failure_header)
        print(f"\n{failure_header}")
        for filename, error in failed_files:
            failure_detail = f"{filename}: {error}"
            logger.warning(f"  {failure_detail}")
            print(f"  {failure_detail}")
    
    if not combined_data:
        error_msg = "No data successfully loaded - cannot create combined file!"
        logger.error(error_msg)
        print(error_msg)
        return
    
    # Combining phase
    combine_header = "=== COMBINING DATA ==="
    logger.info(combine_header)
    print(f"\n{combine_header}")
    
    files_with_data = len(combined_data)
    files_empty = len(empty_files)
    files_failed = len(failed_files)
    
    success_rate = f"Files with data: {files_with_data}/{len(files)}"
    empty_rate = f"Empty files skipped: {files_empty}/{len(files)}"
    fail_rate = f"Failed files: {files_failed}/{len(files)}"
    total_msg = f"Total results: {total_results:,}"
    regions_msg = f"Unique regions: {len(set(regions_processed))}"
    contexts_msg = f"Study contexts: {len(study_contexts)}"
    
    for msg in [success_rate, empty_rate, fail_rate, total_msg, regions_msg, contexts_msg]:
        logger.info(msg)
        print(msg)
    
    # Combine all dataframes
    try:
        concat_msg = "Concatenating dataframes..."
        logger.info(concat_msg)
        print(concat_msg)
        
        combined_df = pd.concat(combined_data, ignore_index=True, sort=False)
        
        shape_msg = f"Combined shape: {combined_df.shape[0]:,} rows × {combined_df.shape[1]} columns"
        logger.info(shape_msg)
        print(shape_msg)
        
    except Exception as e:
        error_msg = f"Error during concatenation: {e}"
        logger.error(error_msg)
        print(error_msg)
        return
    
    # Data quality checks
    quality_header = "=== DATA QUALITY CHECKS ==="
    logger.info(quality_header)
    print(f"\n{quality_header}")
    
    # Check for missing values in key columns based on file type
    if file_type == 'twas':
        key_cols = ['chr', 'molecular_id', 'twas_z', 'twas_pval']
        pval_col = 'twas_pval'
    elif file_type == 'mr':
        key_cols = ['exposure', 'outcome', 'method', 'b', 'se', 'pval']
        pval_col = 'pval'
    else:
        key_cols = []
        pval_col = None
    
    for col in key_cols:
        if col in combined_df.columns:
            missing_count = combined_df[col].isnull().sum()
            missing_pct = (missing_count / len(combined_df)) * 100
            missing_msg = f"Missing values in {col}: {missing_count:,} ({missing_pct:.1f}%)"
            logger.info(missing_msg)
            print(missing_msg)
    
    # Check for duplicates (excluding metadata columns)
    data_cols = [col for col in combined_df.columns 
                 if col not in ['region', 'source_file', 'study_context']]
    duplicates = combined_df.duplicated(subset=data_cols).sum()
    if duplicates > 0:
        dup_msg = f"Duplicate results: {duplicates:,}"
        logger.warning(dup_msg)
        print(dup_msg)
        logger.info("    (Consider removing duplicates if not expected)")
        print("    (Consider removing duplicates if not expected)")
    
    # Summary statistics
    if pval_col and pval_col in combined_df.columns:
        valid_pvals = combined_df[pval_col].dropna()
        if len(valid_pvals) > 0:
            sig_count = (valid_pvals < 0.05).sum()
            sig_msg = f"Significant results (p < 0.05): {sig_count:,} ({sig_count/len(valid_pvals)*100:.1f}%)"
            logger.info(sig_msg)
            print(sig_msg)
    
    # Create output filename
    if output_name:
        output_file = f"{output_name}.combined_{file_type}.tsv.gz"
    else:
        output_file = f"combined_{file_type}_results.tsv.gz"
    
    output_path = os.path.join(directory, output_file)
    
    # Save combined data
    save_header = "=== SAVING RESULTS ==="
    logger.info(save_header)
    print(f"\n{save_header}")
    
    output_msg = f"Output file: {output_file}"
    logger.info(output_msg)
    print(output_msg)
    
    try:
        with gzip.open(output_path, 'wt') as f:
            combined_df.to_csv(f, sep='\t', index=False, na_rep='NA')
        
        save_success = f"Successfully saved {len(combined_df):,} results to {output_file}"
        logger.info(save_success)
        print(save_success)
        
        # File size info
        file_size = os.path.getsize(output_path) / (1024 * 1024)  # MB
        size_msg = f"Output file size: {file_size:.1f} MB"
        logger.info(size_msg)
        print(size_msg)
        
    except Exception as e:
        error_msg = f"Error saving file: {e}"
        logger.error(error_msg)
        print(error_msg)
        return
    
    # Final summary
    summary_header = "=== SUMMARY ==="
    logger.info(summary_header)
    print(f"\n{summary_header}")
    
    # Summary with empty file information
    summary_stats = [
        f"Total files processed: {len(files)}",
        f"Files with data: {files_with_data}",
        f"Empty files skipped: {files_empty}",
        f"Failed files: {files_failed}",
        f"Total results combined: {total_results:,}"
    ]
    
    for stat in summary_stats:
        logger.info(stat)
        print(stat)
    
    if 'region' in combined_df.columns:
        region_counts = combined_df['region'].value_counts()
        region_summary = f"Regions with data: {len(region_counts)}"
        logger.info(region_summary)
        print(region_summary)
        
        if len(region_counts) <= 10:
            for region, count in region_counts.items():
                region_detail = f"  {region}: {count:,} results"
                logger.info(region_detail)
                print(region_detail)
        else:
            logger.info("Top 5 regions:")
            print("Top 5 regions:")
            for region, count in region_counts.head(5).items():
                region_detail = f"  {region}: {count:,} results"
                logger.info(region_detail)
                print(region_detail)
    
    if 'study_context' in combined_df.columns:
        context_counts = combined_df['study_context'].value_counts()
        context_summary = f"Study contexts: {len(context_counts)}"
        logger.info(context_summary)
        print(context_summary)
        
        for context, count in context_counts.items():
            context_detail = f"  {context}: {count:,} results"
            logger.info(context_detail)
            print(context_detail)
    
    # Log completion
    completion_msg = f"Process completed successfully at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    logger.info(completion_msg)
    logger.info(f"Full log saved to: {os.path.basename(log_path)}")
    print(f"\nFull log saved to: {os.path.basename(log_path)}")

def list_file_patterns(directory):
    """List all TWAS and MR file patterns found in directory"""
    twas_files = find_result_files(directory, 'twas')
    mr_files = find_result_files(directory, 'mr')
    
    if not twas_files and not mr_files:
        print(f"No result files (*.twas.tsv.gz or *.mr_result.tsv.gz) found in '{directory}'")
        return
    
    # Display TWAS files
    if twas_files:
        print(f"\n=== TWAS Files ===")
        print(f"Found {len(twas_files)} TWAS files:")
        print(f"{'Filename':<50} {'Region':<20} {'Study Context':<20} {'Status'}")
        print("-" * 95)
        
        empty_count = 0
        for filepath in twas_files[:20]:  # Show first 20 files
            filename = os.path.basename(filepath)
            region = extract_region_from_filename(filename) or "unknown"
            study_context = extract_study_context_from_filename(filename, 'twas')
            
            # Check if file is empty
            is_empty, _ = is_file_empty(filepath)
            status = "EMPTY" if is_empty else "OK"
            if is_empty:
                empty_count += 1
                
            print(f"{filename:<50} {region:<20} {study_context:<20} {status}")
        
        if len(twas_files) > 20:
            print(f"... and {len(twas_files) - 20} more files")
        
        print(f"\n  Total TWAS files: {len(twas_files)}")
        print(f"  Empty TWAS files: {empty_count}")
    
    # Display MR files
    if mr_files:
        print(f"\n=== MR Files ===")
        print(f"Found {len(mr_files)} MR files:")
        print(f"{'Filename':<50} {'Region':<20} {'Study Context':<20} {'Status'}")
        print("-" * 95)
        
        empty_count = 0
        for filepath in mr_files[:20]:  # Show first 20 files
            filename = os.path.basename(filepath)
            region = extract_region_from_filename(filename) or "unknown"
            study_context = extract_study_context_from_filename(filename, 'mr')
            
            # Check if file is empty
            is_empty, _ = is_file_empty(filepath)
            status = "EMPTY" if is_empty else "OK"
            if is_empty:
                empty_count += 1
                
            print(f"{filename:<50} {region:<20} {study_context:<20} {status}")
        
        if len(mr_files) > 20:
            print(f"... and {len(mr_files) - 20} more files")
        
        print(f"\n  Total MR files: {len(mr_files)}")
        print(f"  Empty MR files: {empty_count}")
    
    # Show summary
    print(f"\n=== Overall Summary ===")
    print(f"Total TWAS files: {len(twas_files)}")
    print(f"Total MR files: {len(mr_files)}")
    print(f"Total files: {len(twas_files) + len(mr_files)}")

def main():
    parser = argparse.ArgumentParser(
        description='Combine TWAS (*.twas.tsv.gz) and/or MR (*.mr_result.tsv.gz) result files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine TWAS files
  python combine_results.py /path/to/folder --type twas
  
  # Combine MR files
  python combine_results.py /path/to/folder --type mr
  
  # Combine both TWAS and MR files separately
  python combine_results.py /path/to/folder --type both
  
  # Combine with custom output name
  python combine_results.py /path/to/folder --type twas --output my_results
  
  # List all available files
  python combine_results.py /path/to/folder --list-patterns
        """
    )
    parser.add_argument('directory', help='Directory containing result files')
    parser.add_argument('--type', '-t', choices=['twas', 'mr', 'both'], default='twas',
                       help='Type of files to combine (default: twas)')
    parser.add_argument('--output', '-o', help='Output filename prefix (without extension)')
    parser.add_argument('--list-patterns', action='store_true', 
                       help='List all result file patterns in the directory')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.directory):
        print(f"Error: Directory '{args.directory}' does not exist")
        sys.exit(1)
    
    if not os.path.isdir(args.directory):
        print(f"Error: '{args.directory}' is not a directory")
        sys.exit(1)
    
    if args.list_patterns:
        list_file_patterns(args.directory)
        return
    
    combine_result_files(args.directory, args.type, args.output)

if __name__ == "__main__":
    main()