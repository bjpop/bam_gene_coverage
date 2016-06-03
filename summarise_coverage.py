from itertools import groupby
import numpy

# XXX put list of sample IDs in here
samples = []

suffix = ".sorted.bam.genes.coverage"

def proportion(predicate, values):
    num_values = len(values)
    count = 0
    for x in values:
        if predicate(x):
            count += 1
    return float(count) / float(num_values)

def median(values):
    num_values = len(values)
    mid = values[num_values/2]
    return mid

def average(values):
    num_values = len(values)
    return float(sum(values) / float(num_values))

def range(values):
    values = sorted(values)
    low = values[0]
    high = values[-1]
    return low, high

def get_gene(line):
    return line.split()[3]

def get_stats(gene, size, values):
    num_values = len(values)
    assert(num_values == size)
    low, high = range(values)
    return[
        median(values),
        average(values),
        numpy.std(numpy.array(values)),
        low,
        high,
        proportion(lambda x: x == 0, values),
        proportion(lambda x: x >= 5, values),
        proportion(lambda x: x >= 10, values),
        proportion(lambda x: x >= 20, values),
        proportion(lambda x: x >= 40, values)]

NUM_GENES = 34
header = "sample, " + ", ".join(["median", "average", "stddev", "low", "high", "'==0", "'>=5", "'>=10", "'>=20", "'>=40"] * NUM_GENES)

def main():
    print(header)
    for sample in samples:
        filename = sample + suffix
        results_for_sample = []
        with open(filename) as file:
            for gene, group in groupby(file, key=get_gene):
                values = []
                for line in group:
                    fields = line.split()
                    chrom, start, end, gene, index, coverage = fields[:6]
                    start = int(start)
                    end = int(end)
                    coverage = int(coverage)
                    values.append(coverage)
                size = (end - start) 
                results_for_sample.extend(get_stats(gene, size, values))
        stats_str = ", ".join(["{:.2f}".format(x) for x in results_for_sample])
        print("{}, {}".format(sample, stats_str))

if __name__ == '__main__':
    main()
