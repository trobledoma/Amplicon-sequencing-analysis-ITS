03_README

Venn diagram of shared taxa

R version 4.4.2 (2024-10-31 ucrt)
packageVersion("dplyr")
‘1.2.1’
packageVersion("tidyr")
‘1.3.2’
packageVersion("stringr")
‘1.6.0’
packageVersion("VennDiagram")
‘1.8.2’
packageVersion("grid")
‘4.4.2’

Each analysis can be performed separately

Be sure the format of the input files follow this format.

-input files required:

feature_table.csv
metadata.csv



Consideration for the feature table (csv format)

Simplified example:
Feature Kingdom Phylum Genus Species sample_001 sample_002 sample_003
ASV_001 k__Fungi p__Ascomycota g__Alternaria 120 85 101
ASV_002 k__Fungi p__Basidiomycota g__Entoloma 30 42 25
ASV_003 k__Fungi p__Ascomycota 10 5

- The taxonomic column name must be exactly "Genus" (case-sensitive).
- Sample names must appear as column headers.
- Abundance values ​​must be non-negative numbers, preferably read counts.
- Unidentified genera may be left blank or labeled with terms such as "g__", "unclassified", or "uncultured". The pipeline will group these under "Unclassified at genus level".
- You may retain other taxonomic columns, such as Kingdom, Phylum, Class, Order, Family, and Species.

Considerations for the metadata (csv format)

This is the metadata table, with one row per sample.
sample_name	location	pH	Temperature
sample_001	Control	7.2	12.5
sample_002	Control	7.4	13.0
sample_003	Treatment	6.8	15.1

These values ​​are for illustrative purposes.

Requirements:
- sample_name: An identifier for each sample. It must match the corresponding column name in `feature_table.csv` exactly.
- location: A label used by the pipeline as a "Condition" for grouping and comparing results.
- Each identifier must appear only once.
- Both `sample_name` and `location` must contain valid values; no empty cells are permitted.
- You may include environmental variables such as pH, temperature, conductivity, nutrient levels, depth, or coordinates. There is no need to remove them from the original file.
Important: If "location" refers solely to a geographic site and you wish to compare a different experimental factor—such as a treatment or habitat type—you will need to adapt the pipeline to use the appropriate column.