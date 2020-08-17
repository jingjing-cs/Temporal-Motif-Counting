###
### Efficient Sampling Algorithms for Approximate Temporal Motif Counting
###

### All algorithms were implemented in C++11 compiled by GCC v7.4 with -O3 optimizations, and ran on a single thread.


## For Linux

Here is an example of how to use the code:
```
# Download the code
form https://github.com/jingjing-hnu/Temporal-Motif-Counting 

# Build the cod
cd ES_EWS
make

# Run all algorithms es, es_circle4, ews

# Example
./es ../datasets/example.txt ../motifs/motif11.txt out.txt 300 5

# Parameters
es: Please input five parameters: temporal graph file, motif file, outfile, delta, edge sampling P

es_circle4: please input four parameter: infile, outfile, delta, edge sampling P

ews: please input seven parameters: infile, motiffile, outfile, delta, edge sampling P, wedge sampling Q,seed for P





