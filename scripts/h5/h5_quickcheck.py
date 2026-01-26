#h5ls -r -v gamma_10evt.h5
#h5ls -r gamma_10evt.h5
#python -c "import h5py; f=h5py.File('gamma_100evt.h5','r'); print(list(f.keys())); [print(k, f[k].shape, f[k].dtype) for k in f.keys()]"

#!/usr/bin/env python3
"""Quick check of HDF5 file contents"""
import sys

filename = sys.argv[1] if len(sys.argv) > 1 else "Gamma1.h5"

try:
    import h5py
    print(f"Checking {filename} with h5py...")
    with h5py.File(filename, 'r') as f:
        print(f"\nFile: {filename}")
        print("="*60)
        print("Top-level keys:", list(f.keys()))
        print("\nDetails:")
        for key in f.keys():
            obj = f[key]
            print(f"\n  {key}:")
            print(f"    Type: {type(obj).__name__}")
            if hasattr(obj, 'shape'):
                print(f"    Shape: {obj.shape}")
            if hasattr(obj, 'dtype'):
                print(f"    Dtype: {obj.dtype}")
                if hasattr(obj.dtype, 'names') and obj.dtype.names:
                    print(f"    Columns: {obj.dtype.names}")
                    # Show first row if it's a table
                    if len(obj.shape) > 0 and obj.shape[0] > 0:
                        print(f"    First row: {dict(zip(obj.dtype.names, obj[0]))}")
except ImportError:
    try:
        import tables as tb
        print(f"Checking {filename} with PyTables...")
        hdf = tb.open_file(filename, 'r')
        print(f"\nFile: {filename}")
        print("="*60)
        for node in hdf.walk_nodes():
            if isinstance(node, tb.Table):
                print(f"\n  {node._v_pathname}:")
                print(f"    Rows: {node.nrows}")
                print(f"    Columns: {node.colnames}")
                if node.nrows > 0:
                    print(f"    First row: {dict(node[0])}")
        hdf.close()
    except ImportError:
        print("Error: Need h5py or PyTables installed")
        print("Install with: pip install h5py  or  pip install tables")
        sys.exit(1)
except FileNotFoundError:
    print(f"Error: File '{filename}' not found!")
    print(f"Current directory files:")
    import os
    for f in os.listdir('.'):
        if f.endswith('.h5'):
            print(f"  {f}")
    sys.exit(1)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
