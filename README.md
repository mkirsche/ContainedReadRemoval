# Contained Read Removal

Faster Assemblies with Minimizer-Based Removal of Contained Reads

## Running the filter:

To run the filter, use the executable rect.sh in this directory.

./rect.sh <len_filter> <shared_minimizer_filter> <output_prefix>

## Parameters:

len_filter is the length threshold L such that all reads with length less than L are removed

shared_minimizer_filter is the threshold proportion of shared minimizers such that if a read has that many shared kmers with some other read (both overall and in its ends), the read is removed

## Description of the directories:

assembly_eval is used to assemble filtered datasets with existing assemblers wtdbg2 and Canu as well as evaluate individual assemblies produced by the pipeline


assemblycompare is used to compare assemblies produced by different sets of parameters

figures contains images illustarting the performance of the filter

sim contains the code for running and evaluating the performance of the filter on simulated data

src contains the core source code of the filter


