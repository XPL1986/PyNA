# PyNA

This pipeline utilizes Genome Analysis Toolkit's best practices workflow for single 
nucleotide polymorphism and insertion/deletions as the basic framework to improve upon. 
The main contribution of this pipeline includes fine search via variant allele frequency,
biotype, annotations, impact, transcript length, and mapping quality associated with each 
variant to increase the overall sensitivity and specificity. Along with variant detection, 
this pipeline also supports the identification of variants demonstrating either allele-specific 
expression or loss of heterozygosity properties. 

Calling Variants in RNAseq \
https://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq

## Features

* Identification of high impact mutations in RNA sequence data.
* Support for allele-specific expression or loss of heterozygosity property detection.
* Pipelines Raw DNA and RNA sequence to Annotated VCF format.
* Detects consequences in resulting peptide in single nucleotide polymorphism.
* Query by name, location, overlapping regions in  GRCh38/hg38.

### Prerequisites

```
Burrows-Wheeler Aligner (BWA) 
http://bio-bwa.sourceforge.net/

Genomic Analysis Toolkit (GATK) version 3.8
https://software.broadinstitute.org/gatk/documentation/version-history.php?id=10063&page=2

Picard
https://broadinstitute.github.io/picard/

STAR: ultrafast universal RNA-seq aligner
https://academic.oup.com/bioinformatics/article/29/1/15/272537

SnpEff
http://snpeff.sourceforge.net/

SnpSift
http://snpeff.sourceforge.net/SnpSift.html
```

### Installation


```
git clone https://github.com/lcwong0928/PyNA.git PyNA
```


## Usage

Query by name about gene of interest based on GRCh38/hg38.

```
pyna name [-h] [-o OUTPUT] [-f FEATURE] gene [genes ...]

# [-o OUTPUT] saves as txt file at given directory
# [-f FEATURE] CDS, exon, five_prime_utr, gene, Selenocysteine,
                  start_codon, stop_codon, three_prime_utr, transcript
```


Query by location based on GRCh38/hg38.

```
pyna location [-h] [-o OUTPUT] coord [coord ...]

# coord formated as chr:coord, i.e. chr2:208236227 (IDH1)
# [-o OUTPUT] saves as txt file at given directory
```


Query by overlapping regions based on GRCh38/hg38.

```
pyna overlap [-h] [-o OUTPUT] coord [coord ...]

# coord formated as chr:coord-coord, i.e. chr2:208236227-208266074 (IDH1)
# [-o OUTPUT] saves as txt file at given directory
```



### Results

Will be documented shortly.

```
Example
```


## Built With

* [PyCharm](https://www.jetbrains.com/pycharm/) - Python IDE used

## Contributing

Please read [CONTRIBUTING.md](https://github.com/lcwong0928/PyNA/blob/master/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.


## Authors

* **Lawrence C. Wong** - *Initial work* - [lcwong0928](https://github.com/lcwong0928)

See also the list of [contributors](https://github.com/lcwong0928/PyNA/graphs/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Jiguang Wang, Principal Investigator, The Hong Kong University of Science and Technology
* Biaobin Jiang, Postdoc Fellow, The Hong Kong University of Science and Technology
* Quanhua Mu, PhD Student, The Hong Kong University of Science and Technology

This computational pipeline was made possible by the International Undergraduate Research Project 
program between Massachusetts Institute of Technology and The Hong Kong University of Science and Technology. 
I would like to express my gratitude to Professor Jiguang Wang at HKUST for his supervision over and support 
on this project. I would also like to thank Postdoc Fellow Biaobin Jiang, PhD Student Quanhua Mu, and members 
of the Wang Lab for their guidance and assistance on this research project.