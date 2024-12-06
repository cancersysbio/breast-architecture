### FIGURE 3D: UCD65 IC2 ECDNA ENRICHMENT IN LINEARLY DIGESTED SAMPLES ############################
# plot avg. log ratio of coverage (digested:parental) for UCD65 IC2 ecDNA against distribution of random ecDNA

### PREAMBLE ######################################################################################

import random
import numpy as np
import matplotlib.pyplot as plt

CHR_REGIONS = {'1': (1, 248956422),
               '2': (1, 242193529),
               '3': (1, 198295559),
               '4': (1, 190214555),
               '5': (1, 181538259),
               '6': (1, 170805979),
               '7': (1, 159345973),
               '8': (1, 145138636),
               '9': (1, 138394717),
               '10': (1, 133797422),
               '11': (1, 135086622),
               '12': (1, 133275309),
               '13': (1, 114364328),
               '14': (1, 107043717),
               '15': (1, 101991189),
               '16': (1, 90338345),
               '17': (1, 83257441),
               '18': (1, 80373285),
               '19': (1, 58617616),
               '20': (1, 64444167),
               '21': (1, 46709983),
               '22': (1, 50818468),
               'X': (1, 156040895)}
BIN_SIZE = 50000

### SUPPORT FUNCTIONS ###

# convert info from raw read count file into average coverage depth dictionary
def convert_raw_read_count_file(raw_read_count_file_path):
    avg_coverage_depth_dict = {}
    with open(raw_read_count_file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            fields = line.strip().split("\t")
            chrom = fields[1]
            start = int(fields[2])
            end = int(fields[3])
            avg_coverage = int(fields[4]) * 150/50000
            if chrom not in avg_coverage_depth_dict:
                avg_coverage_depth_dict[chrom] = []
            avg_coverage_depth_dict[chrom].append((start, end, avg_coverage))
    return avg_coverage_depth_dict

# calculate the median coverage over all regions in an average coverage depth dictionary
def get_median_coverage(avg_coverage_depth_dict):
    all_coverage_values = []
    for chrom in avg_coverage_depth_dict:
        regions = avg_coverage_depth_dict[chrom]
        for region in regions:
            all_coverage_values.append(region[2])
    median_coverage = np.median(all_coverage_values)
    return median_coverage

# condense a list of genomic region tuples into non-overlapping regions
def condense_tuples_list(list):
    merged_list = [list[0]]
    for current_region in list[1:]:
        last_region = merged_list[-1]
        if current_region[0] <= last_region[1]:
            merged_list[-1] = (last_region[0], max(last_region[1], current_region[1]))
        else:
            merged_list.append(current_region)
    return merged_list

# get the copy count from a string of information
def parse_copy_count(string, start_indicator, stop_indicator):
    idx1 = string.index(start_indicator) + 1
    idx2 = idx1 + 1
    while (string[idx2] != stop_indicator):
        idx2 += 1
    return int(string[idx1:idx2])

# make a dictionary to store the segments included on each cycle of ecDNA
def get_cycles(file, amplif_threshold=4):
    f = open(file)
    cycles_dict = {}
    for line in f:
        line = line.strip()
        if line[0] == 'C':
            cycle_info = line.split(";")
            copy_ct = parse_copy_count(cycle_info[1], '=', '.')
            if copy_ct >= amplif_threshold:
                cycle_num_idx = cycle_info[0].index('=') + 1
                cycle_num = cycle_info[0][cycle_num_idx:]
                segments_idx = cycle_info[2].index('=') + 1
                segments = cycle_info[2][segments_idx:]
                cycles_dict[cycle_num] = segments
    return cycles_dict

# make a dictionary to store the regions for each ecDNA segment
def get_segments(file):
    f = open(file)
    all_segments_dict = {}
    for line in f:
        line = line.strip()
        if line[0] == "S":
            segment_info = line.split("\t")
            segment_number = segment_info[1]
            segment_chromosome = segment_info[2]
            segment_start = int(segment_info[-2])
            segment_end = int(segment_info[-1])
            all_segments_dict[segment_number] = [segment_chromosome, segment_start, segment_end]
    return all_segments_dict

# extract all ecDNA genomic regions from the cycle file into a dictionary, stored by chromosome
def get_ecDNA_regions_from_cycle_file(cycle_file):
    amplified_cycles_dict = get_cycles(cycle_file)
    amp_cycle_list = []
    for cycle in amplified_cycles_dict:
        cycle_info = amplified_cycles_dict[cycle]
        cycle_info_list = cycle_info.split(',')
        segment_numbers = [segment[:-1] for segment in cycle_info_list]
        for segment in segment_numbers:
            if segment not in amp_cycle_list and segment != '0':
                amp_cycle_list.append(segment)
    segments_dict = get_segments(cycle_file)
    ecDNA_regions_full = {}
    for segment in amp_cycle_list:
        segment_info = segments_dict[segment]
        segment_chrom = segment_info[0][3:]
        if segment_chrom not in ecDNA_regions_full:
            ecDNA_regions_full[segment_chrom] = []
        ecDNA_regions_full[segment_chrom].append((int(segment_info[1]), int(segment_info[2])))
    for chrom in ecDNA_regions_full:
        ecDNA_regions = ecDNA_regions_full[chrom]
        ecDNA_regions_condensed = condense_tuples_list(ecDNA_regions)
        ecDNA_regions_full[chrom] = ecDNA_regions_condensed
    return ecDNA_regions_full

# calculate the length of an ecDNA from its regions dictionary
def get_ecDNA_length_from_regions(ecDNA_regions_dict):
    segment_lengths = []
    for chrom in ecDNA_regions_dict:
        segments = ecDNA_regions_dict[chrom]
        for segment in segments:
            print(segment)
            segment_lengths.append(int(segment[1]) - int(segment[0]))
            print(str(int(segment[1]) - int(segment[0])))
    return segment_lengths

# get the average coverage across a full raw read count file
def get_average_coverage(raw_read_count_file_path):
    total_read_count = 0
    total_length = 0
    with open(raw_read_count_file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            fields = line.strip().split("\t")
            total_read_count += int(fields[4])
            total_length += 50000
    average_coverage = total_read_count * 150 / total_length
    return average_coverage

# get the average coverage across ecDNA intervals
def get_average_coverage_on_ecDNA_intervals(avg_coverage_depth_dict, ecDNA_intervals):

    total_ecDNA_length = 0
    summed_coverage_depth = 0

    for chrom in ecDNA_intervals:
        intervals = ecDNA_intervals[chrom]
        if chrom in avg_coverage_depth_dict:
            for interval in intervals:
                start = interval[0]
                end = interval[1]
                total_ecDNA_length += end - start + 1
                chrom_coverage = avg_coverage_depth_dict[chrom]
                for sample_start, sample_end, sample_avg_coverage in chrom_coverage:
                    # this coverage interval is before the ecDNA interval
                    if sample_end < start:
                        continue
                    # this coverage interval is after the ecDNA interval
                    elif sample_start > end:
                        continue
                    # this coverage interval starts before, but overlaps with, the ecDNA interval
                    elif sample_start <= start <= sample_end:
                        # the end is in the middle of the ecDNA interval
                        if sample_end <= end:
                            summed_coverage_depth += sample_avg_coverage * (sample_end - start + 1)
                        # the end goes beyond the ecDNA interval
                        elif sample_end > end:
                            summed_coverage_depth += sample_avg_coverage * (end - start + 1)
                    # this coverage interval starts after, but overlaps with, the ecDNA interval
                    elif start < sample_start:
                        # the end goes beyond the ecDNA interval
                        if sample_end >= end:
                            summed_coverage_depth += sample_avg_coverage * (end - sample_start + 1)
                        # the end is inside the ecDNA interval
                        elif sample_end < end:
                            summed_coverage_depth += sample_avg_coverage * (sample_end - sample_start + 1)

    average_coverage_depth_on_ecDNA = summed_coverage_depth / total_ecDNA_length

    return average_coverage_depth_on_ecDNA

# get a list of blacklisted regions to exclude from null ecDNA generation
def get_blacklist_regions(chr):
    blacklist_regions = []
    blacklist_file = "encode_blacklist_regions.bed"
    with open(blacklist_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            chrom = fields[0]
            if chrom == 'chr' + str(chr):
                blacklist_regions.append((int(fields[1]), int(fields[2])))
    return blacklist_regions

# get a list of all genomic regions excluding blacklisted regions and the predicted ecDNA regions
def get_valid_genome_regions(chr, full_ecDNA_regions):
    unmappable_regions = get_blacklist_regions(chr)
    if chr in full_ecDNA_regions:
        ecDNA_region = full_ecDNA_regions[chr]
        unmappable_regions.append(ecDNA_region)
        unmappable_regions = sorted(unmappable_regions, key=lambda x: x[0])
        unmappable_regions = condense_tuples_list(unmappable_regions)
    mappable_regions = []
    if unmappable_regions[0][0] != 0:
        unmappable_start = unmappable_regions[0][0]
        mappable_regions.append((0, unmappable_start - 1))
    start = unmappable_regions[0][1] + 1
    for interval in unmappable_regions[1:]:
        end = interval[0] - 1
        mappable_regions.append((start, end))
        start = interval[1] + 1
    if start != CHR_REGIONS[chr][1]:
        mappable_regions.append((start, CHR_REGIONS[chr][1]))
    return mappable_regions

# get a random region of the genome of a given segment length
def get_random_region(genome_regions, segment_length):
    # Calculate the total length of all ranges combined
    total_range_length = sum(end - start + 1 for start, end in genome_regions)

    # Generate a random starting point within the total length
    random_start = random.randint(0, total_range_length - 1)

    # Find the range that contains the random starting point
    for start, end in genome_regions:
        range_length = end - start + 1
        if random_start > range_length:
            random_start -= range_length
        else:
            starting_bp = start + random_start
            if starting_bp + segment_length < end:
                return (starting_bp, starting_bp + segment_length)

    # if you make it here without returning, this has failed to find a region so we try again
    return get_random_region(genome_regions, segment_length)

# generate a random null ecDNA with equal length and segmenting to an ecDNA prediction
def generate_null_ecDNA(segment_lengths_list, full_ecDNA_regions):
    null_segments = {}
    for segment_length in segment_lengths_list:
        # randomly select a chromosome
        chromosome = str(random.randint(1, 22))
        # get the valid genome regions for the chromosome
        valid_genome_regions = get_valid_genome_regions(chromosome, full_ecDNA_regions)
        # generate ecDNA segments
        null_segment = get_random_region(valid_genome_regions, segment_length)
        if chromosome not in null_segments:
            null_segments[chromosome] = []
        null_segments[chromosome].append(null_segment)
    return null_segments


### MAIN ##########################################################################################

def main():
    UCD12_raw_read_count_file_digested = 'UCD12_EC_raw_readcounts.txt'
    UCD12_raw_read_count_file_parental = 'UCD12_Prt_raw_readcounts.txt'
    UCD12_cycle_file = 'UCD12_amplicon2_cycles.txt'
    UCD12_full_ecDNA_intervals = {'8': (30150001, 42825889), '11': (67048236, 67069082)}

    plot_null_distribution_ratio(UCD12_raw_read_count_file_digested, UCD12_raw_read_count_file_parental,
                                 UCD12_cycle_file, 'UCD12', UCD12_full_ecDNA_intervals)

    UCD65_raw_read_count_file_digested = 'UCD65_EC_raw_readcounts.txt'
    UCD65_raw_read_count_file_parental = 'UCD65_Prt_raw_readcounts.txt'
    UCD65_cycle_file = 'UCD65_amplicon1_cycles.txt'
    UCD65_full_ecDNA_intervals = {'8': (35150001, 42223758), '12': (122300001, 124400000), '17': (27600001, 74490000),
                                  '19': (5684617, 7550000), '20': (57396374, 57417204)}

    plot_null_distribution_ratio(UCD65_raw_read_count_file_digested, UCD65_raw_read_count_file_parental,
                                 UCD65_cycle_file, 'UCD65', UCD65_full_ecDNA_intervals)


# plot figure 3d
def plot_null_distribution_ratio(sample_raw_read_count_file_digested, sample_raw_read_count_file_parental, comparison_sample_cycle_file, sample_name, full_ecDNA_regions):
    # get digested ACD dict, median coverage for normalization
    avg_coverage_depth_dict_digested = convert_raw_read_count_file(sample_raw_read_count_file_digested)
    median_cov_digested = get_median_coverage(avg_coverage_depth_dict_digested)

    # get parental ACD dict, median coverage for normalization
    avg_coverage_depth_dict_parental = convert_raw_read_count_file(sample_raw_read_count_file_parental)
    median_cov_parental = get_median_coverage(avg_coverage_depth_dict_parental)

    # get predicted ecDNA segments dict
    predicted_ecDNA_segments_dict = get_ecDNA_regions_from_cycle_file(comparison_sample_cycle_file)

    # get predicted ACD on ecDNA for digested
    predicted_avg_cov_depth_digested = get_average_coverage_on_ecDNA_intervals(
        avg_coverage_depth_dict_digested,
        predicted_ecDNA_segments_dict) / median_cov_digested

    # get predicted ACD on ecDNA for parental
    predicted_avg_cov_depth_parental = get_average_coverage_on_ecDNA_intervals(
        avg_coverage_depth_dict_parental,
        predicted_ecDNA_segments_dict) / median_cov_parental

    # calculate predicted ACD ratio
    predicted_avg_cov_depth_ratio = predicted_avg_cov_depth_digested / predicted_avg_cov_depth_parental

    ecDNA_segment_lengths = get_ecDNA_length_from_regions(predicted_ecDNA_segments_dict)

    # start collecting null ACDs
    null_average_coverage_ratios = []

    for i in range(1000):
        null_ecDNA_segments = generate_null_ecDNA(ecDNA_segment_lengths, full_ecDNA_regions)
        average_coverage_digested = get_average_coverage_on_ecDNA_intervals(avg_coverage_depth_dict_digested,
                                                                                      null_ecDNA_segments) / median_cov_digested
        average_coverage_parental = get_average_coverage_on_ecDNA_intervals(avg_coverage_depth_dict_parental,
                                                                                      null_ecDNA_segments) / median_cov_parental
        null_average_coverage_ratios.append(average_coverage_digested / average_coverage_parental)

    p_val = 1 - len([num for num in null_average_coverage_ratios if num < predicted_avg_cov_depth_ratio]) / len(null_average_coverage_ratios)
    print(p_val)
    log_null_average_coverage_ratios = [np.log(value) for value in null_average_coverage_ratios]
    plt.figure(figsize=(14, 5))
    plt.hist(log_null_average_coverage_ratios, bins=50, color='black', edgecolor='black')
    plt.axvline(x=np.log(predicted_avg_cov_depth_ratio), color='red', linestyle='--')
    plt.xlabel('Log of Average Coverage Ratio (digested:parental)', fontsize=18)
    plt.ylabel('Frequency', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.xlim(-3.5, 3.75)
    plt.ylim(0, 100)
    plt.legend(fontsize=14)
    plt.show()


if __name__ == '__main__':
    main()