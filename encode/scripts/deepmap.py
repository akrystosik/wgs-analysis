#!/usr/bin/env python3

import pandas as pd
import sys

def verify_files():
    print("=== CCLE Mutations File Analysis ===")
    try:
        # Read mutations file
        mutations = pd.read_csv('data/deepmap/CCLE_mutations.csv')
        print("\nColumns in mutations file:")
        print(mutations.columns.tolist())
        
        print(f"\nTotal number of mutation records: {len(mutations)}")
        
        print("\nFirst mutation record:")
        print(mutations.iloc[0])
        
        # Check for key columns we expect
        expected_cols = ['DepMap_ID', 'Hugo_Symbol', 'Chromosome', 'Start_position']
        missing_cols = [col for col in expected_cols if col not in mutations.columns]
        if missing_cols:
            print("\nWarning: Missing expected columns:", missing_cols)
        else:
            print("\nAll expected mutation columns present")

        print("\n=== Sample Info File Analysis ===")
        # Read sample info file
        samples = pd.read_csv('data/deepmap/sample_info.csv')
        print("\nColumns in sample info file:")
        print(samples.columns.tolist())
        
        # Look for K562
        k562_samples = samples[samples.astype(str).apply(lambda x: x.str.contains('K562', case=False)).any(axis=1)]
        print(f"\nFound {len(k562_samples)} K562 entries")
        if not k562_samples.empty:
            print("\nK562 entries:")
            print(k562_samples)
            
        return True
        
    except FileNotFoundError as e:
        print(f"Error: Could not find file - {str(e)}")
        return False
    except Exception as e:
        print(f"Error analyzing files: {str(e)}")
        return False

if __name__ == "__main__":
    success = verify_files()
    sys.exit(0 if success else 1)