# scramEnricher
Identify regions of enriched siRNAs or miRNAs in scramAligner CSV files

## Releases
You can download the pre-built binaries from the [Releases](https://github.com/username/scramEnricher/releases) page.

## Installation
To build from source (when go is installed), you can follow these steps:

1. Clone the repository:
    ```
    git clone https://github.com/username/scramEnricher.git
    ```

2. Change to the scramEnricher directory:
    ```
    cd scramEnricher
    ```

3. Build the executable:
    ```
    go build
    ```

4. Run the executable:
    ```
    ./scramEnricher -h
    ```

## Usage
Once you have installed scramEnricher, you can use it with the following command line options:

- `-input`: Comma-separated list of input CSV files.
- `-max-times-aligned`: Maximum times aligned cutoff (default 5).
- `-merge-distance`: Merge distance for windows (default 100).
- `-min-rpmr`: Minimum average RPMR (default 10).
- `-min-unique`: Minimum number of unique sRNAs (default 5).
- `-output`: Output CSV file (default "enriched_windows.csv").
- `-output-fasta`: Output FASTA file (optional).
- `-reference`: FASTA reference file (optional).
- `-window`: Window size (default 100).

