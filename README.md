# McQueen

This is an MCMC program to estimate selection across the genome in the presence of recombination. 
We consider a large reference dataset, and determine the likelihood of observing a much smaller query dataset.

## Getting Started 

To run this program on OS X, there are a number of requirements:

1.  Install Xcode (from the app store).
2.  Install command line tools (found within the downloads tab under preferences in Xcode).
3.  Ensure [X11](http://xquartz.macosforge.org/landing/) is installed.
4.  Install wget

    * Install wget manually:
    
    ```
    curl -O http://ftp.gnu.org/gnu/wget/wget-1.13.4.tar.gz
    tar -xzf wget-1.13.4.tar.gz
    cd wget-1.13.4
    ./configure --with-ssl=openssl
    make
    sudo make install
    ```

    * Install wget with [MacPorts](http://www.macports.org/install.php):
    
    ```
    sudo port install wget
    ```

    Once this is done, run the Makefile in the [libs](https://github.com/astheeggeggs/mcqueen/tree/master/libs) directory:

    ```
    cd libs
    make
    ```
    
5. Finally, run the Makefile in the [mcqueen](https://github.com/astheeggeggs/mcqueen) directory:

    ```
    cd ..
    make
    ```

There are three main programs, `dmodel`, `dsim` and `dgen`:

* `dmodel` performs MCMC to estimate HLA associated selection along a region.

* `dsim` simulates data under this model.
* `dgen` simulates data under a birth-death process.

Summary of how to use the functions can be obtained by typing `./dmodel`, `./dsim`, and `./dgen` respectively. In more detail:

## dmodel

Usage: `dmodel [options] <ref_seqs.fasta> <query.fasta> <outdir>`

#### Required arguments:

* `ref_seqs.fasta`: Path to reference sequences (saved as a `.fasta` file).
* `query.fasta`: Path to query sequences (saved as a `.fasta` file).
*  If the `-q` flag is used, `ref_seqs.fasta` and `query.fasta` should be the same.

* `outdir`: Path to the directory in which files generated from this run are to be stored.
    * If the directory doesn't yet exist, it is generated.

#### Options:

A long or short flag can be used to pass each of the options below:

* `-H --HLA_inference <HLA.csv>` pass a .csv file with HLA information.

    The collection of HLA types associated to the members of `query.fasta`.
    Samples must be in the same order as the sequences in `query.fasta`.
    The HLA names can be any character string, no specific format is required.
    The format of this `.csv` file should be as follows (a larger example is shown in [test_files](https://raw.githubusercontent.com/astheeggeggs/mcqueen/master/test_files/query.HLA.csv)).

    e.g.

    ```
    "","HLA_1","HLA_2","HLA_3","HLA_4"
    "seq_name_1",0,0,1,0
    "seq_name_2",0,0,1,0
    "seq_name_3",0,0,0,0
    ```
*  `-n --chain_length <N>`: Set chain length _[default: 50,000]_.
    
    An integer: The number of states you wish the MCMC to run for.
* `-t --runtime_hours <T>`: run for `T` hours.
    
    An integer: An alternative to `-n`, instead setting the time limit on the MCMC run.
* `-s --sample_every <S>`: Sample every `S` steps _[default: 100]_.
    
    An integer: How often do you wish to write a state to file.
* `-c --consensus_fasta <C>`: Consensus set as consensus of this `.fasta` file _[default: <ref_seqs.fasta>]_.
    
    A filepath to a collection of sequences, or single sequence. The 'wild type' strain is set to the consensus sequence of this collection of sequences. If the `-c` flag is not used, the 'wild type' is set as the consensus of the collection of reference sequences (sequences in `ref_seqs.fasta`).
* `-r --resume_json <state.json>`: Resume state.
    
    If you wish to continue an MCMC run from the end of a previous run, pass the `.json` file that was generated at the end of that run, the program will set this as the initial state.
* `-l --closest_n <L>`: Closest L sequences _[default: 100]_.
    
    It is impossible to consider all members of the reference set (`ref_seqs.fasta`) if it is very large. Using Hamming distance as a distance metric, we consider the closest L sequences to each member of the query set (`query.fasta`). 
    
    If the `-q` flag is passed, then the closest L sequences excluding the query sequence currently under consideration are used.
* `-p --parameters_json <parameters.json>`: `.json` file containing parameter information.

    A `.json` file containing the codon usage, transition transversion ratio and proportion of codon changes which involve at least two nucleotide substitutions, if the user requires different parameters to the hardcoded defaults.
* `-m --true_mosaic_fasta <M>`: True mosaic of sequences from simulation:
    Used for part of a simulation study. If the collection of sequences are known to be recombinants from members of the reference set, they can be passed. This was used to test the validity of the Hamming distance metric in the presence of high recombination rates.
* `-q --separate_reference_fasta <Q>`: Separate reference set _[default: True]_.
    
    Flag for whether we consider a separate reference set or not. If `-q` is passed, them the reference sequence file (`ref_seqs.fasta`) should be the same as the query sequence file (`query.fasta`). When determining the closest L sequences, we don't consider the query sequence itself.
* `-i --no_HLA_inference`: No HLA associated selection inference is to be performed _[default: False]_.
    
    Use this flag to turn off inference of HLA associated selection parameters in the model. This will turn off all window moves in the MCMC scheme, and set the HLA associated selection to 1 for all HLA types.
* `-j --no_rev_inference`: No reversion inference _[default: False]_.

    Use this flag to turn off inference of the reversion scaling. This will turn off all reversion window moves in the MCMC scheme, and set the reversion scaling to 1 at each site. 

* `-o --no_recombination_inference`: No recombination _[default: False]_.

* `-b --sample_prior`: Sample all parameters from the specified priors _[default: False].

* `-z --no_simulated_annealing`: Don't perform simulated annealing during the initial portion of the MCMC _[default: False]_.
    
    Use this flag to turn off inference of the recombination probability between each pair of sites. When flagged, no recombination moves are proposed in the MCMC scheme and the recombination probability is set to 0 between all neighbouring sites.

The following options allow the user to easily switch between a small set of parameter priors in the inference:

* `-v --omega_gamma`: Gamma prior on omega coefficients _[default: exponential prior]_.

* `-e --omega_log_normal`: Log Normal prior on omega coefficients _[default: exponential prior]_.

* `-w --omega_uniform`: Improper uniform prior on omega coefficients _[default: exponential prior]_.

* `-g --HLA_coeff_esc_gamma`: Gamma prior on HLA coefficients _[default: log normal prior]_.

* `-k --HLA_coeff_esc_normal`: Normal prior on HLA coefficients _[default: log normal prior]_.

* `-f --HLA_coeff_esc_uniform`: Flat prior on HLA coefficients _[default: log normal prior]_.

* `-a --coeff_rev_gamma`: Gamma prior on HLA coefficients _[default: log normal prior]_.

* `-q --coeff_rev_normal`: Normal prior on HLA coefficients _[default: log normal prior]_.

* `-u --coeff_rev_uniform`: Flat prior on HLA coefficients _[default: log normal prior]_.

#### Estimated parameters:
* Selection coefficients _dN/dS_ at each site in the genome (these are not HLA associated).

* Probability of recombination between every pair of sites (unless using the `-o` flag).

* Synonymous transition rate, mu.

* Selection parameters within windows for every HLA type associated to the collection of query sequences, which we denote 'escape' coefficients (non-synonymous base changes away from some defined consensus strain), unless using the `-i` flag.

* Selection parameters independent of HLA type which we denote 'reversion' coefficients (non-synonymous base changes towards some defined consensus strain), unless using the `-i` or `-j` flags.

#### Move proposals:
At each iteration of the MCMC one of the following collection of moves is proposed:

1.  Mutation move: Change the overall mutation rate (through a scaling move).

2.  Selection move: Change the selection coefficient (_dN/dS_) at a site chosen uniformly at random through a scaling move.

3.  Recombination move (ignored when using the `-o` flag): Change a recombination rate between a pair of sites (chosen uniformly at random) through a scaling move.

4.  Window moves (ignored when using the `-h` flag):

    *  Change the selection parameter associated to escape (reversion) within a window chosen uniformly at random from the collection of windows.
    *  Merge two windows.
    *  Split a window.
    *  Extend a window.

## dsim

Used for the first simulation study in our paper, this function simulates query sequence data under our modelling regime. Query sequences are generated through recombination and mutation, where mutation is dependent on host HLA profiles.

Usage: `dsim [options] <parameters.json> <ref_seqs.fasta> <outdir>`

#### Required arguments:
* `parameters.json`: path to `.json` file containing parameters of the simulation run. The following parameters are required in the `.json` file, and the relevant error is thrown if they are not present. An example is provided in the `test_files` folder. _Note: the functions in `creating_simulation_json.r` allow the user to automatically generate a correctly formatted `.json` file for use in dgen or dsim_.

    * `kappa`: The transition/transversion ratio.
    * `mu`: The mutation rate of synonymous transversions.
    * `codon_sequence_length`: The length of the region to be simulated in codons.
    * `n_HLA`: A vector of the number of HLA types in each HLA gene (must be of length n_genes).
    * `total_n_HLA`: The total number of HLA types.
    * `ploidy`: The number of each HLA gene that an individual has.
    * `n_genes`: The number of HLA genes that an individual has.
    * `HLA_prevalences`: A vector of HLA prevalences for the `total_n_HLA` HLA types (must be of length `total_n_HLA`).
    * `sites`: A struture of length `codon_sequence_length` which contains the following enties for each codon:
        * `omega`: The _dN/dS_ ratio at the site.
        * `R`: The probability of recombination between adjacent sites.
        * `reversion`: Selection back to wildtype from amino-acids a non-synonymous single base change from the wildtype strain.
        * `HLA`: A vector of length `total_n_HLA` that defines selection away from wildtype to amino-acids that are a single non-synonymous base change away from the wildtype at this codon.
        
* `ref_seqs.fasta`: Path to reference sequences (saved as a `.fasta` file).
  
* `outdir`: Path to the directory in which files generated from this run are to be stored.
    If the directory doesn't yet exist, it is generated.
#### Options:
* `-s --n_sims <S>`: Set number of simulated sequences _[default: 500]_.
* `-c --consensus_fasta <C>`: Set a consensus sequence _[default: consensus of ref_seqs.fasta]_.
    A filepath to a collection of sequences, or single sequence. The 'wild type' strain is set to the consensus sequence of this collection of sequences. If the `-c` flag is not used, the 'wild type' is set as the consensus of the collection of reference sequences (sequences in `ref_seqs.fasta`).
* `-n --prop_no_HLA <N>`: Set a proportion of the simulated query sequences to not have associated HLA information _[default: 0.2]_.

## dgen

Used for the second simulation study in our paper, we consider a generative birth-death process to create reference and query sequence data to test our method. This function creates an instance of a birth-death process according to the specified parameters, and simulates sequence data down the resultant ancestry.

Usage: `dgen [options] <parameters.json> <N> <M> <lambda> <mu_tree> <outdir>`

### Required arguments:
* `parameters.json`: path to `.json` file containing parameters of the simulation run. The following parameters are required in the `.json` file, and the relevant error is thrown if they are not present. An example is provided in the `test_files` folder. _Note: the functions in `creating_simulation_json.r` allow the user to automatically generate a correctly formatted `.json` file for use in dgen or dsim_.

    * `kappa`: The transition/transversion ratio.
    * `mu`: The mutation rate of synonymous transversions.
    * `codon_sequence_length`: The length of the region to be simulated in codons.
    * `n_HLA`: A vector of the number of HLA types in each HLA gene (must be of length n_genes).
    * `total_n_HLA`: The total number of HLA types.
    * `ploidy`: The number of each HLA gene that an individual has.
    * `n_genes`: The number of HLA genes that an individual has.
    * `HLA_prevalences`: A vector of HLA prevalences for the `total_n_HLA` HLA types (must be of length `total_n_HLA`).
    * `sites`: A struture of length `codon_sequence_length` which contains the following enties for each codon:
        * `omega`: The _dN/dS_ ratio at the site.
        * `R`: The probability of recombination between adjacent sites.
        * `reversion`: Selection back to wildtype from amino-acids a non-synonymous single base change from the wildtype strain.
        * `HLA`: A vector of length `total_n_HLA` that defines selection away from wildtype to amino-acids that are a single non-synonymous base change away from the wildtype at this codon.

* `N`: The sampled number of individuals at the present. 

* `M`: The total number of individuals at the present.

* `lambda`: The birth rate of the birth-death process generating the tree (death rate backwards in time).

* `mu_tree`: The death rate of the birth-death process generating the tree (birth rate backwards in time).

* `outdir`: : Path to the directory in which files generated from this run are to be stored.
    If the directory doesn't yet exist, it is generated.

### Optional arguments:
* `-s --prop_past_sampling <S>`: Set sampling proportion for historic sampling _[default: M/N]_.

* `-c --consensus_fasta <C>`: Pass `.fasta` file from which to generate the consensus sequence.

* `-r --root_fasta <R>`: Pass `.fasta` file giving the sequence at the root node.

* `-n --write_newick_tree`: Write the resultant tree to a Newick tree file.

* `-g --write_tree_matrix`: Write the resultant tree to a tree file.

* `-q --num_queries`: Provide the number of queries sequences you wish to generate.

* `-f --prop_query`: Provide the query fraction _[default: 0.01]_.
