# Bacmeta 
Version 1.1.0


Bacmeta is free for use, under the terms of BSD 3-clause License. See LICENSE.txt. 


It has been described in article: 

Aleksi Sipola, Pekka Marttinen, Jukka Corander; Bacmeta: simulator for genomic evolution in bacterial metapopulations, Bioinformatics, , bty093, https://doi.org/10.1093/bioinformatics/bty093

Please cite the article if you are publishing any work done with help of Bacmeta!

--- 

Developer welcomes all forms of feedback, help inquiries, developing suggestions and general correspondence about the program via email!



aleksi.sipola 

helsinki 

fi

---
##  WHAT?

Simulator for forward-simulating evolution of bacterial genomes in multiple, diversely connected populations. 

Simulation starts with a single genotype, which can consist of one or multiple loci of desired length, and proceeds in discrete time, "generation-by-generation". Single nucleotide mutations, indels of Zipf distributed length and recombinations of whole loci or segments of geometrically distributed lengths are modelled as Poisson processes. Migration and microepidemics are also modelled as Poisson processes, with options for migration probabilities and microepidemic sizes. Metadata about each simulated level and object, i.e. locus, individual, population or metapopulation, can be considered and data of various aspects are provided, together with the main output files of sequences and loci-wise pairwise distances.

See `ModelManual.pdf` and `InputOutputManual.pdf` for further information.

---
##  HOW?


### Preparing: 

Only requirements are C++11 and g++ compiler. Windows users need cygwin.

Download bacmeta from ` bitbucket.com/aleksisipola/bacmeta `. Extract and place the acquired folder `/bacmeta/` into a preferred location. Note that all the output files will be created into the subfolder `***/bacmeta/outputs/`. All source code files and make file for compilation are included in `/bacmeta/src/`. Subfolder `/bacmeta/outputs/` will be created automatically when the simulator is run for the first time.

 
### Compiling: 
	
Simulator must be compiled with the included makefile to create an executable simulator program `simu`. This needs to be done only once, or when updating. Open a terminal and cd into `***/bacmeta/src/`, and enter `make`. Successful make should result in text: `g++ Bac.o Gene.o Para.o Pop.o Simu.o Superpop.o -o ../simu`. After compiling, all files needed for running the simulator, i.e. input files and simulator executable itself, are located in main folder `/bacmeta/`.

 
### Setting up simulation: 

Parameters for the simulation run are set with the plain text file:  

`simu.input`

or with help of command line options, variation of the filename can be used:

`simu123.input` , 

where in place of `123` any short, alphanumeric string can be used. Use of the  command line options is described further in this README. Normal text editor or command line programs such as `sed` can be used to edit the input files. All parameter options are shown in `default.params` and those lines can be copied into and modified in `simu.input`. Note that `default.params` should not be edited as it is always used for initializing all necessary parameters. Furthermore, only contents of `simu.input` are documented for simulation runs.

When using non-identical migration rates between populations, the rates are set with file:

`migration.input`

or, like with general parameters file, custom name can be used with command line options.

Note that, if parameters are set without consideration, output files might get excessively numerous or large.

See `InputOutputManual.pdf` for more information.


### Simulating:

Simulation is run on terminal in the main folder `/bacmeta/` with command 

`./simu`  

 or with command line options e.g.  
 
`./simu -p 123 -m 321 -o 213`


If you want to **STOP** the simulation, press `Ctrl+c`.

Notice that any previous output files with same filename modifier parameter value are overwritten! 

In the beginning Bacmeta goes through initializing, after which simulation progress is shown in intervals of 1000 generations. Successful run will end with line: `Simulation finished and output files created.`



Previous simulations can be replicated using the same seed value that is documented into `rundetails***.txt` using the same platform. However, the results can differ if runs are made in different platforms. This is due to the freedom of random number generation implementations in the C++ standard. 
    
    

####Command line options:  

  `-p`   Specify 'simu---.input' file for parameter input.  The argument replaces '---'. 
        For example './simu -p 123' would use simu123.input.  
        
  `-m`   Specify 'migration---.input' file for  migration rate matrix input.  
        For example './simu -p 123 -m 333' would use migration333.input (and simu123.input).  
        
  `-b`   Specify both of the files above with same suffix.  
        For example './simu -b 123' would use simu123.input AND migration123.input.  
        
  `-o`   Specify modifier for the output filenames.  
        For example './simu -p 123 -o 456' would use simu123.input BUT output files would have filename modifier 456, ignoring the value in simu123.input.  


  Arguments work at least with numbers and alphabetic letters. Symbols have NOT been tested and are not recommended.  


  (These command line option descriptions are also shown with the command `./simu -h`. No simulation is run when using `-h` command line option.) 



### Outputs:

Variety of output files can be created by a simulation run. Output files are created into folder `/bacmeta/outputs/`, which is created automatically if one does not already exist. All files have the filename modifier suffix given in `simu.input` or with the flag `-o`. Some of the output files can have multiple occurrences at intervals. To exclude certain documentation, give negative number for documentation interval, or toggle off with `0` according to the comments in input file. 

#### Note:

These files can reach relatively large sizes if parameters for documentation are not set accordingly for big simulations. 


More thorough descriptions of outputs are given in `InputOutputManual.pdf`, and examples are provided in `FeatureExample.zip`. 



---

## WHERE?

This program is written for and tested in UNIX-environment. Main developing platform currently is Linux Mint 18 (Sarah). 

Current version has also been tested successfully on OS X Mavericks and Windows 10 with cygwin. Though current version has not been tested on these, previous versions have run well on Debian (Jessie), Ubuntu, OS X Yosemite and Windows 7 and 10 both with cygwin and VMware player running Ubuntu virtual machine.



---


