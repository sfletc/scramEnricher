package main

import (
	"bufio"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"sync"
)

type sRNAData struct {
	AlignmentHeader string
	Len             int
	sRNA            string
	Position        int
	Strand          string
	TimesAligned    int
	Replicates      []float64
}

type window struct {
	start int
	end   int
}

func main() {
	// Define command-line flags
	inputFiles := flag.String("input", "", "Comma-separated list of input CSV files")
	outputFile := flag.String("output", "enriched_windows.csv", "Output CSV file")
	windowSize := flag.Int("window", 100, "Window size")
	minUniqueSRNAs := flag.Int("min-unique", 5, "Minimum number of unique sRNAs")
	minAvgRPMR := flag.Float64("min-rpmr", 10.0, "Minimum average RPMR")
	maxTimesAligned := flag.Int("max-times-aligned", 5, "Maximum times aligned cutoff")
	mergeDistance := flag.Int("merge-distance", 100, "Merge distance for windows")
	referenceFile := flag.String("reference", "", "FASTA reference file (optional)")
	outputFastaFile := flag.String("output-fasta", "", "Output FASTA file (optional)")
	flag.Parse()

	// Validate command-line flags
	if *inputFiles == "" {
		fmt.Println("Please provide input files using the -input flag")
		os.Exit(1)
	}

	// Parse the list of input files
	fileList := strings.Split(*inputFiles, ",")
	// Parse each SCRAM output CSV file
	sRNADataByFile := make(map[string]map[string][]sRNAData)
	for _, file := range fileList {
		sRNADataByHeader := make(map[string][]sRNAData)
		parseSCRAMFile(file, sRNADataByHeader)
		sRNADataByFile[file] = sRNADataByHeader
	}

	// Collect unique headers from all files
	uniqueHeaders := make(map[string]bool)
	sRNADataByHeader := make(map[string][]sRNAData)
	for _, sRNADataByFile := range sRNADataByFile {
		for header, sRNAData := range sRNADataByFile {
			uniqueHeaders[header] = true
			sRNADataByHeader[header] = append(sRNADataByHeader[header], sRNAData...)
		}
	}

	// Process each unique header concurrently across all files
	var wg sync.WaitGroup
	enrichedRegionsByHeader := make(chan map[string][]window)

	for header := range uniqueHeaders {
		wg.Add(1)
		go func(header string) {
			defer wg.Done()
			regions := processHeader(header, sRNADataByHeader[header], *windowSize, *minUniqueSRNAs, *minAvgRPMR, *maxTimesAligned)
			enrichedRegionsByHeader <- regions
		}(header)
	}

	go func() {
		wg.Wait()
		close(enrichedRegionsByHeader)
	}()

	// Collect enriched regions from all headers
	allEnrichedRegions := make(map[string][]window)
	for regionsByHeader := range enrichedRegionsByHeader {
		for header, regions := range regionsByHeader {
			allEnrichedRegions[header] = append(allEnrichedRegions[header], regions...)
		}
	}
	// Merge adjacent enriched windows for each header across all files
	mergeWindows(allEnrichedRegions, *mergeDistance)

	// Load reference sequences from FASTA file
	referenceSeqs := make(map[string]string)
	if *referenceFile != "" {
		err := loadReferenceSequences(*referenceFile, referenceSeqs)
		if err != nil {
			fmt.Printf("Error loading reference sequences: %v\n", err)
		}
	}

	// Extract sequences in enriched windows and write to FASTA file
	if *outputFastaFile != "" && len(referenceSeqs) > 0 {
		err := writeEnrichedSequencesToFasta(allEnrichedRegions, referenceSeqs, *outputFastaFile)
		if err != nil {
			fmt.Printf("Error writing enriched sequences to FASTA file: %v\n", err)
		}
	}

	// Write windows to CSV
	writeWindowsToCSV(allEnrichedRegions, sRNADataByHeader, *outputFile)

	// Print the total number of windows written to CSV
	totalWindows := 0
	for _, regions := range allEnrichedRegions {
		totalWindows += len(regions)
	}
	fmt.Printf("\nTotal windows written to CSV: %d\n", totalWindows)
}

func writeWindowsToCSV(enrichedRegions map[string][]window, sRNADataByHeader map[string][]sRNAData, outputFile string) {
	file, err := os.Create(outputFile)
	if err != nil {
		fmt.Printf("Error creating CSV file: %s\n", err)
		return
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write CSV header
	header := []string{"Header", "Start", "End", "Length"}
	writer.Write(header)

	for header, regions := range enrichedRegions {
		for _, region := range regions {
			sRNAData := sRNADataByHeader[header]
			seqLength := sRNAData[len(sRNAData)-1].Position
			end := min(region.end, seqLength)
			length := end - region.start

			row := []string{
				header,
				strconv.Itoa(region.start),
				strconv.Itoa(end),
				strconv.Itoa(length),
			}
			writer.Write(row)
		}
	}
}

func parseSCRAMFile(file string, sRNADataByHeader map[string][]sRNAData) {
	log.Printf("Parsing file: %s\n", file)
	f, err := os.Open(file)
	if err != nil {
		fmt.Printf("Error opening file: %s\n", err)
		return
	}
	defer f.Close()

	reader := csv.NewReader(f)

	for {
		row, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			fmt.Printf("Error reading CSV: %s\n", err)
			return
		}

		sRNA := sRNAData{
			AlignmentHeader: row[0],
			Len:             atoi(row[1]),
			sRNA:            row[2],
			Position:        atoi(row[3]),
			Strand:          row[4],
			TimesAligned:    atoi(row[5]),
		}

		for i := 6; i < len(row); i++ {
			rpmr, _ := strconv.ParseFloat(row[i], 64)
			sRNA.Replicates = append(sRNA.Replicates, rpmr)
		}

		sRNADataByHeader[sRNA.AlignmentHeader] = append(sRNADataByHeader[sRNA.AlignmentHeader], sRNA)
	}
}

func processHeader(header string, sRNAData []sRNAData, windowSize, minUniqueSRNAs int, minAvgRPMR float64, maxTimesAligned int) map[string][]window {
	// Count unique sRNAs and calculate average RPMR per window
	windowCounts := make(map[window]struct {
		uniqueCount int
		avgRPMR     []float64
	})
	for _, sRNA := range sRNAData {
		if sRNA.TimesAligned > maxTimesAligned {
			continue // Skip sRNAs that exceed the 'times aligned' cutoff
		}

		windowStart := (sRNA.Position / windowSize) * windowSize
		windowEnd := windowStart + windowSize
		w := window{windowStart, windowEnd}

		counts := windowCounts[w]
		if counts.uniqueCount == 0 {
			counts.avgRPMR = make([]float64, len(sRNA.Replicates))
		}
		counts.uniqueCount++
		for i, rpmr := range sRNA.Replicates {
			counts.avgRPMR[i] += rpmr
		}
		windowCounts[w] = counts
	}

	// Apply thresholds to identify enriched windows
	var enrichedWindows []window
	for w, counts := range windowCounts {
		enriched := counts.uniqueCount >= minUniqueSRNAs
		for _, rpmr := range counts.avgRPMR {
			if rpmr < minAvgRPMR {
				enriched = false
				break
			}
		}
		if enriched {
			enrichedWindows = append(enrichedWindows, w)
		}
	}

	return map[string][]window{header: enrichedWindows}
}

func atoi(s string) int {
	i, _ := strconv.Atoi(s)
	return i
}

func mergeWindows(enrichedRegions map[string][]window, mergeDistance int) {
	for header, windows := range enrichedRegions {
		sort.Slice(windows, func(i, j int) bool {
			return windows[i].start < windows[j].start
		})

		var merged []window
		var lastMerged *window
		for i := 0; i < len(windows); i++ {
			w := windows[i]
			if lastMerged == nil {
				lastMerged = &w
			} else {
				if w.start-lastMerged.end <= mergeDistance {
					lastMerged.end = w.end
					windows = append(windows[:i], windows[i+1:]...)
					i--
				} else {
					merged = append(merged, *lastMerged)
					lastMerged = &w
				}
			}
		}

		if lastMerged != nil {
			merged = append(merged, *lastMerged)
		}

		enrichedRegions[header] = merged
	}
}

func loadReferenceSequences(file string, referenceSeqs map[string]string) error {
	f, err := os.Open(file)
	if err != nil {
		return fmt.Errorf("error opening reference file: %v", err)
	}
	defer f.Close()

	scanner := bufio.NewScanner(f)
	var header string
	var seq strings.Builder

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if header != "" {
				referenceSeqs[header] = seq.String()
				seq.Reset()
			}
			header = strings.TrimPrefix(line, ">")
		} else {
			seq.WriteString(line)
		}
	}

	if header != "" {
		referenceSeqs[header] = seq.String()
	}

	return nil
}

func writeEnrichedSequencesToFasta(enrichedRegions map[string][]window, referenceSeqs map[string]string, outputFile string) error {
	f, err := os.Create(outputFile)
	if err != nil {
		return fmt.Errorf("error creating output FASTA file: %v", err)
	}
	defer f.Close()

	writer := bufio.NewWriter(f)

	for header, regions := range enrichedRegions {
		refSeq, ok := referenceSeqs[header]
		if !ok {
			continue
		}

		for i, region := range regions {
			start := region.start
			end := min(region.end, len(refSeq))
			seq := refSeq[start:end]

			fmt.Fprintf(writer, ">%s_region_%d\n", header, i+1)
			fmt.Fprintf(writer, "%s\n", seq)
		}
	}

	writer.Flush()
	return nil
}
