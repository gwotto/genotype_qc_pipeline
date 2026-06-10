#!/usr/bin/env python3
# Fix the WDL file

with open('/Users/gwo/devel/genotype_qc_pipeline/src/genotype_qc_preimputation.wdl', 'r') as f:
    content = f.read()

# Replace the problematic section
old_section = """    String relatedness_log_line = select_first([
        LogStep10.line,
        "Step 10   Relatedness check (KING - samples NOT removed)            [skipped]"
    }
    ])"""

new_section = """    String relatedness_log_line = select_first([
        LogStep10.line,
        "Step 10   Relatedness check (KING - samples NOT removed)            [skipped]"
    ])
    }"""

content = content.replace(old_section, new_section)

with open('/Users/gwo/devel/genotype_qc_pipeline/src/genotype_qc_preimputation.wdl', 'w') as f:
    f.write(content)

print("Fixed WDL file")
