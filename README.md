[![Continuous Integration](https://github.com/Redmar-van-den-Berg/test-umi-deduplication/actions/workflows/ci.yml/badge.svg)](https://github.com/Redmar-van-den-Berg/test-umi-deduplication/actions/workflows/ci.yml)
[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)
![GitHub release](https://img.shields.io/github/v/release/redmar-van-den-berg/test-umi-deduplication)
![Commits since latest release](https://img.shields.io/github/commits-since/redmar-van-den-berg/test-umi-deduplication/latest)

# test-umi-deduplication
A pipeline to test different methods of UMI deduplication

## Installation
Download the repository from github
```bash
git clone https://github.com/Redmar-van-den-Berg/test-umi-deduplication.git
```

Install and activate the
[conda](https://docs.conda.io/en/latest/miniconda.html)
environment.
```bash
conda env create --file environment.yml
conda activate test-umi-deduplication
```

## Settings
There are three levels where configuration options are set, in decreasing order
of priority.
1. Flags passed to snakemake using `--config`, or in the specified
   `--configfile`.
2. Setting specified in the PEP project configuration, under the key
   `test-umi-deduplication`
3. The default settings for the pipeline, as specified in the `common.smk` file


### Supported settings
| Option               | Type                        | Explanation                                       |
| ---------------------| --------------------------- | ------------------------------------------------- |
| repeat               | Optional integer (default=1)| How many times the benchmark should be repeated   |
| star_index           | Required folder             | STAR-index folder                                 |
| gtf                  | Required file               | GTF file with transcripts, used by STAR           |
| reference            | Required file               | Fasta file for the reference, used by Picard      |
| reference_dict       | Required file               | Dict file for the reference, used by Picard       |
| cluster_method       | "maximum" or "directional" (default) | Clustering method to use                 |
| word_size            | Optional integer (default=24) | Word size for HUMID                             |

## Default and extended mode
There are currently two ways to run the pipeline. By default, it will generate
the most common and useful output files.

By passing the `--until extended` option to Snakemake, the extended version of
the pipeline will run, which adds additional analysis which are usually not
needed:
1. Run UMI-Tools on the output of HUMID, and compare the results
2. Run HUMID on the output of UMI-Tools, and compare the results
3. Analyse and compare the UMI usage for each sample against all possible UMIs
   for a given legth.
4. Generate a separate MultiQC report for 1 and 2

## Tests
You can run the tests that accompany this pipeline with the following commands

```bash
# Check if requirements are installed, and run linting on the Snakefile
pytest --kwd --tag sanity

# Test the pipeline settings in dry-run mode
pytest --kwd --tag dry-run

# Test the performance of the pipeline by running on the test data
pytest --kwd --tag integration
```

## Troubleshooting
If HUMID crashes when using the `maximum` clustering method, try increasing the
stack_size_kb` parameter of the `humid` rule.
