#!/usr/bin/env python3
"""
Inspect HDF5 file structure and show all available variables/tables.
"""

import argparse
import sys

def inspect_with_h5py(filename):
    """Inspect HDF5 file using h5py"""
    try:
        import h5py
    except ImportError:
        return False
    
    print(f"\n{'='*60}")
    print(f"Inspecting: {filename}")
    print(f"{'='*60}\n")
    
    with h5py.File(filename, 'r') as f:
        def print_structure(name, obj):
            indent = '  ' * (name.count('/'))
            if isinstance(obj, h5py.Dataset):
                shape = obj.shape
                dtype = obj.dtype
                print(f"{indent}{name}")
                print(f"{indent}  Type: Dataset")
                print(f"{indent}  Shape: {shape}")
                print(f"{indent}  Dtype: {dtype}")
                
                # Try to show column names if it's a structured array
                if dtype.names:
                    print(f"{indent}  Columns: {', '.join(dtype.names)}")
                    # Show first few rows as sample
                    if len(shape) > 0 and shape[0] > 0:
                        sample_size = min(3, shape[0])
                        print(f"{indent}  Sample data (first {sample_size} rows):")
                        data = obj[:sample_size]
                        for i, row in enumerate(data):
                            print(f"{indent}    Row {i}: {dict(zip(dtype.names, row))}")
                else:
                    # Show sample values for simple arrays
                    if len(shape) > 0 and shape[0] > 0:
                        sample_size = min(5, shape[0])
                        print(f"{indent}  Sample values: {obj[:sample_size]}")
                print()
            elif isinstance(obj, h5py.Group):
                print(f"{indent}{name}/")
                print(f"{indent}  Type: Group")
                print()
        
        print("File Structure:")
        print("-" * 60)
        f.visititems(print_structure)
        
        # Summary
        print("\n" + "="*60)
        print("Summary:")
        print("="*60)
        datasets = []
        groups = []
        def collect_info(name, obj):
            if isinstance(obj, h5py.Dataset):
                datasets.append((name, obj.shape, obj.dtype))
            elif isinstance(obj, h5py.Group):
                groups.append(name)
        f.visititems(collect_info)
        
        print(f"\nTotal Groups: {len(groups)}")
        print(f"Total Datasets: {len(datasets)}")
        print(f"\nDatasets:")
        for name, shape, dtype in datasets:
            print(f"  {name}: shape={shape}, dtype={dtype}")
    
    return True

def inspect_with_pytables(filename):
    """Inspect HDF5 file using PyTables"""
    try:
        import tables
    except ImportError:
        return False
    
    print(f"\n{'='*60}")
    print(f"Inspecting: {filename}")
    print(f"{'='*60}\n")
    
    hdf = tables.open_file(filename, 'r')
    
    print("File Structure:")
    print("-" * 60)
    
    def print_node(node, indent=0):
        spaces = '  ' * indent
        if isinstance(node, tables.Table):
            print(f"{spaces}{node._v_pathname}")
            print(f"{spaces}  Type: Table")
            print(f"{spaces}  Rows: {node.nrows}")
            print(f"{spaces}  Columns: {', '.join(node.colnames)}")
            
            # Show column info
            print(f"{spaces}  Column details:")
            for colname in node.colnames:
                col = node.cols[colname]
                print(f"{spaces}    {colname}: dtype={col.dtype}, shape={col.shape}")
            
            # Show sample data
            if node.nrows > 0:
                sample_size = min(3, node.nrows)
                print(f"{spaces}  Sample data (first {sample_size} rows):")
                for i, row in enumerate(node[:sample_size]):
                    print(f"{spaces}    Row {i}: {dict(row)}")
            print()
        elif isinstance(node, tables.Group):
            print(f"{spaces}{node._v_pathname}/")
            print(f"{spaces}  Type: Group")
            print()
        elif isinstance(node, tables.Array):
            print(f"{spaces}{node._v_pathname}")
            print(f"{spaces}  Type: Array")
            print(f"{spaces}  Shape: {node.shape}")
            print(f"{spaces}  Dtype: {node.dtype}")
            if node.size > 0 and node.size <= 10:
                print(f"{spaces}  Data: {node[:]}")
            print()
    
    # Walk through all nodes
    for node in hdf.walk_nodes():
        level = node._v_pathname.count('/')
        print_node(node, indent=level)
    
    # Summary
    print("\n" + "="*60)
    print("Summary:")
    print("="*60)
    tables_list = []
    arrays_list = []
    for node in hdf.walk_nodes():
        if isinstance(node, tables.Table):
            tables_list.append((node._v_pathname, node.nrows, node.colnames))
        elif isinstance(node, tables.Array):
            arrays_list.append((node._v_pathname, node.shape))
    
    print(f"\nTotal Tables: {len(tables_list)}")
    for path, nrows, cols in tables_list:
        print(f"  {path}: {nrows} rows, columns={cols}")
    
    print(f"\nTotal Arrays: {len(arrays_list)}")
    for path, shape in arrays_list:
        print(f"  {path}: shape={shape}")
    
    hdf.close()
    return True

def inspect_with_h5ls(filename):
    """Try to use command-line h5ls tool"""
    import subprocess
    try:
        result = subprocess.run(['h5ls', '-r', '-v', filename], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            print("\nUsing h5ls command-line tool:")
            print("="*60)
            print(result.stdout)
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    return False

def main():
    parser = argparse.ArgumentParser(
        description='Inspect HDF5 file structure and show all variables',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python inspect_hdf5.py Gammas.h5
  python inspect_hdf5.py Gammas.h5 --method pytables
  python inspect_hdf5.py Gammas.h5 --method h5py
        """
    )
    parser.add_argument('filename', help='HDF5 file to inspect')
    parser.add_argument('--method', choices=['auto', 'h5py', 'pytables', 'h5ls'],
                       default='auto', help='Method to use for inspection')
    args = parser.parse_args()
    
    success = False
    
    if args.method == 'auto':
        # Try h5py first (usually more detailed)
        if inspect_with_h5py(args.filename):
            success = True
        elif inspect_with_pytables(args.filename):
            success = True
        elif inspect_with_h5ls(args.filename):
            success = True
    elif args.method == 'h5py':
        success = inspect_with_h5py(args.filename)
    elif args.method == 'pytables':
        success = inspect_with_pytables(args.filename)
    elif args.method == 'h5ls':
        success = inspect_with_h5ls(args.filename)
    
    if not success:
        print("Error: Could not inspect file.")
        print("Please install one of: h5py, PyTables (tables), or h5ls command-line tool")
        print("\nInstall with:")
        print("  pip install h5py")
        print("  pip install tables")
        print("  # or install HDF5 tools: conda install hdf5")
        sys.exit(1)

if __name__ == "__main__":
    main()
