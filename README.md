# scramEnricher
Identify regions of enriched siRNAs or miRNAs

## Usage

```
scramEnricher -h
Usage of scramEnricher:
  -input string
        Comma-separated list of input CSV files
  -max-times-aligned int
        Maximum times aligned cutoff (default 5)
  -merge-distance int
        Merge distance for windows (default 100)
  -min-rpmr float
        Minimum average RPMR (default 10)
  -min-unique int
        Minimum number of unique sRNAs (default 5)
  -output string
        Output CSV file (default "enriched_windows.csv")
  -output-fasta string
        Output FASTA file (optional)
  -reference string
        FASTA reference file (optional)
  -window int
        Window size (default 100)
```
