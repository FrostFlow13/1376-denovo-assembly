# I highly suggest viewing this document in a program like Notepad++ - it will make navigating/looking at the document a lot nicer.

# Additionally, as a general rule of thumb, if a section of code DOESN'T start with
source /[HOMEPATH]/miniconda3/bin/activate
# that means it's very likely something that needs to be done via the terminal! This is because the jobs need conda itself activated first.

-------------------------------------------------------------------

{Conda Information
# GORT CONDA (installed it into my /[HOMEPATH]/):
	platform : linux-64
	conda version : 23.5.0
	python version : 3.10.12.final.0
	user-agent : conda/23.5.0 requests/2.31.0 CPython/3.10.12 Linux/4.4.0-210-generic ubuntu/16.04.7 glibc/2.23
	channels:
	  - conda-forge
	  - bioconda
	  - defaults
	  channel_priority: strict
	mamba version : 1.4.9


# OHIO SUPERCOMPUTER (OSC (Owens cluster or Pitzer cluster depending on availability)) CONDA (installed it into my /[HOMEPATH]/):
	platform : linux-64
	conda version : 23.5.2
	python version : 3.11.4.final.0
	user-agent : conda/23.5.2 requests/2.29.0 CPython/3.11.4 Linux/3.10.0-1160.90.1.el7.x86_64 rhel/7.9 glibc/2.17
	channels:
	  - conda-forge
	  - bioconda
	  - defaults
	  channel_priority: strict
	mamba version : 1.4.9
}

{Job Setup
#	For all tools, I created separate environments to prevent any weird incompatabilities between packages, following instructions
#	as provided by each group. Pilon and DeconSeq and Clair3 had special installation instructions - see their sections for setup.
	mamba create -n [NAME]_env -c [CHANNELS] [package]
#		NOTE: Medaka more or less REQUIRES mamba for installation! It's incredibly slow installing it via conda.


#	For running conda environments inside OSC (Owens cluster) via [NAMEHERE].sh job scripts:
		#!/bin/bash
		#SBATCH --time=[HH:mm:ss]
		#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=28
		#SBATCH --account=PAS[####]
		#SBATCH --job-name=[NAMEHERE]
		#SBATCH --mail-user=[name.#]@buckeyemail.osu.edu
		#SBATCH --mail-type=BEGIN,END,FAIL
		#SBATCH --export=ALL
		#SBATCH --output=./joblogs/[NAMEHERE].out.%j
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate [NAMEHERE]_env
		cd /[STRAINPATH]/
		SCRIPT HERE

#			For "--time", give it your best estimate. For everything I've done so far, I've set it at "--time=02:00:00", with everything
#			ACTUALLY finishing in 2 hours or less (if I am remembering correctly). Giving it a shorter time means you're more
#			likely to get slotted in, but it also runs the risk of a process getting interrupted.
#				Highly suggest setting it to "--time=02:00:00" - definitely get stuff going sooner, and everything finishes in time.
#				Except for Flye (~3 hours) and Hairsplitter (~7 hours). They need more time.
#			For the line with "nodes" and whatnot, DEFINITELY double check the right arguments/inputs before running. Just to be sure.
#				Plus, if you can use fewer cores/don't need as many, you'll get stuff run sooner!

#	After this, just run the script in the supercomputer via terminal "sbatch [NAME].sh" and it'll set the job up automatically.
#	If necessary, you can use "scancel [JOB#]" and that will halt the job running.
}

-------------------------------------------------------------------

# Some general terms explained:
#	For paths, you must update these manually to your own locations.
	"/[HOMEPATH]/" # default user path with miniconda3  (OSC e.g., "/users/PAS[####]/woodruff207/" ; GORT e.g., "/home/woodruff.207/")
	"/[STRAINPATH]/" # server path with your data in it (OSC e.g., "/fs/ess/PAS[####]/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1376/" ; 
#                                                             GORT e.g., "/srv/work/Andrew_Woodruff/2023_06_15-MAY1376_TLOKOs_LongRead/1376/")
	"[NAME]" # the name of your strain or some identifier you want the files started with

#	Example:
		sbatch script.sh "/users/PAS[####]/woodruff207" "/fs/ess/PAS[####]/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1376" "1376"
#			For directories, EXCLUDE the last '/', since the script should have the next / in it.

	[TOOL] # RAN ON GORT = it was ran on the lab's server as a simple script.
	[TOOL] # RAN ON OSC = it was ran on OSC as a JOB (more complex .sh file, unless stated otherwise).


{Sequencing & Basecalling & Demultiplexing
#	Already performed for us by University of Wisconsin-Madison.

#	RUN SETUP
	Flow cell type			= FLO-PRO114M (R10.4.1)
	Flow cell type alias	= FLO-PRO114M (R10.4.1)
	Flow cell ID			= PAO85317
	Kit type				= SQK-NBD114-24 (V14)
	Q20+ chemistry, essentially?

#	RUN SETTINGS
	Specified run length		= 12 hrs
	Active channel selection	= On
	Pore scan freq.				= 1.5 hrs
	Reserved pores				= On
	Minimum read length			= 200 bp
	Read splitting				= On
	Basecalling					= High-accuracy model, 400 bps
	Modified basecalling		= Off
	Trim barcodes				= Off
	Mid-read barcode filtering	= Off

#	SOFTWARE VERSIONS
	MinKNOW			= 22.12.5
	Bream			= 7.4.8
	Configuration	= 5.4.7
	Guppy			= 6.4.6
	MinKNOW Core	= 5.4.3
	
#	BARCODES
	MAY1376 = BC15
	MAY1739 = BC16 = 1699-1-13-2
	MAY1740 = BC17 = 1700-1-5-1
	MAY1741 = BC18 = 1700-1-9-1
	
#	MISC (from report_PAO85317_20230615_1429_373a141b.md)
	{
	"asic_id": "0004A30B001D61E8", 
	"asic_id_eeprom": "0004A30B001D61E8", 
	"asic_temp": "50.171822", 
	"asic_version": "Unknown", 
	"auto_update": "0", 
	"auto_update_source": "https://cdn.oxfordnanoportal.com/software/MinKNOW/", 
	"bream_is_standard": "0", 
	"configuration_version": "5.4.7", 
	"device_id": "3A", 
	"device_type": "promethion", 
	"distribution_status": "stable", 
	"distribution_version": "22.12.5", 
	"exp_script_name": "sequencing/sequencing_PRO114_DNA_e8_2_400T:FLO-PRO114M:SQK-NBD114-24:400", 
	"exp_script_purpose": "sequencing_run", 
	"exp_start_time": "2023-06-15T14:29:16.558430-05:00", 
	"flow_cell_id": "PAO85317", 
	"flow_cell_product_code": "FLO-PRO114M", 
	"guppy_version": "6.4.6+ae70e8f", 
	"heatsink_temp": "33.683838", 
	"host_product_code": "PRO-PRC024", 
	"host_product_serial_number": "PROMETHION24.SECURE.BIOTECH.WISC.EDU", 
	"hostname": "promethion24.secure.biotech.wisc.edu", 
	"hublett_board_id": "01354ef4eb2e6ac1", 
	"hublett_firmware_version": "2.1.9", 
	"installation_type": "nc", 
	"ip_address": "", 
	"local_firmware_file": "1", 
	"mac_address": "", 
	"operating_system": "ubuntu 20.04", 
	"protocol_group_id": "4101P", 
	"protocol_run_id": "373a141b-9c93-4839-88b1-1929688a1e41", 
	"protocol_start_time": "2023-06-15T14:23:31.685977-05:00", 
	"protocols_version": "7.4.8", 
	"run_id": "79755bc0304042014d0862a0877c7e9b7d4d40da", 
	"sample_id": "4101P", 
	"satellite_board_id": "01345842774e0617", 
	"satellite_firmware_version": "2.1.9", 
	"sequencer_hardware_revision": "", 
	"sequencer_product_code": "PRO-SEQ024", 
	"sequencer_serial_number": "", 
	"usb_config": "fx3_0.0.0#fpga_0.0.0#unknown#unknown", 
	"version": "5.4.3"
	}

#	ADDITIONAL MISC (from report_PAO85317_20230615_1429_373a141b.json)
	basecalling_config_filename: "dna_r10.4.1_e8.2_400bps_hac_prom.cfg"
}


{Adapter Removal
	Porechop_ABI v0.5.0 # RAN ON GORT
		https://github.com/bonsai-team/Porechop_ABI | https://anaconda.org/bioconda/porechop_abi

#		The following will only use the known adapter sequences (not predict any), process a specified .fastq file (the concatenated
#		long-read sequencing file), '-b' demultiplex it again while trimming the adapters AND trim off an extra 8 bases from the ends it
#		does trim (this helps to remove the extra 5'-CAGCACCT-3'(top)/5'-AGGTGCTG-3'(bottom) sequence that are internal from the
#		barcodes on either end), then place the output files into a directory called "1_demul_adtrim".
#			i.e. TOP = 5'-[adapter]CAGCACCT[read]-3' ; BOTTOM = 5'-[read]AGGTGCTG[adapter]-3'
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate porechop_env
		cd /[STRAINPATH]/
		porechop_abi -i [NAMEHERE].fastq -b ./1_demul_adtrim --extra_end_trim 8 -t 28

#		This also uses adapters we identified via consesus (ClustalW) for the SQK-NBD114.24 kit and added into the adapters.py file
#		(we added it directly after the first adapter listed):
            Adapter('SQK-NBD114_24',
                    start_sequence=('SQK-NBD114_24_Y_Top', 'TCGTTCAGTTACGTATTGCTAAGGTTAA'),
                    end_sequence=('SQK-NBD114_24_Y_Bottom', 'TTAACCTTAGCAATACGTA')),
}


{Removal of Smaller Reads
	SeqKit v2.5.0 # RAN ON OSC
		https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit

#		The following will remove all reads smaller than 10/25 kb from the [PORECHOPPED].fastq file in preparation for the
#		Flye assembly (expects 1 kb minimum, prefers 5+ kb)
#		Ideally, you don't want to lose TOO much coverage, so check the remaining bases from the "stats" output and run the math.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate seqkit_env
		cd /[STRAINPATH]/1_demul_adtrim/
		seqkit stats [PORECHOPPED].fastq
		seqkit seq [PORECHOPPED].fastq -m 10000 > [PORECHOPPED]-10kbmin.fastq
		seqkit stats [PORECHOPPED]-10kbmin.fastq
		seqkit seq [PORECHOPPED].fastq -m 25000 > [PORECHOPPED]-25kbmin.fastq
		seqkit stats [PORECHOPPED]-25kbmin.fastq
}


{Genome Assembly
	Flye v2.9.2 & Bandage # RAN ON OSC
		https://github.com/fenderglass/Flye | https://anaconda.org/bioconda/flye & https://rrwick.github.io/Bandage/

#		The following will run Flye using assumptions for ONT data basecalled with Q20 data, grab the porechopped file,
#		then run (Flye collapses the assembly down to haploid if it can, so it's not split between haplotypes). 
#		It'll also use the number of threads listed (be sure your job request has the same number). After running, it places 
#		the resulting directory  up into the folder of the strain being processed for ease of access in the next step. 
#		Probably allocate at least 3 hours in a job for this! Then, use Bandage to visualize .gfa afterwards.
#			***NOTE: Flye apparently expects reads above 1 kb, with 5kb+ being recommended.***
#			***NOTE2: 25 kb makes EXTREMELY consistent graphs with Flye! 
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/1_demul_adtrim/
		flye --nano-hq [PORECHOPPED]-25kbmin.fastq --out-dir ../2_flye_assembly --threads 28
#			The "PORECHOPPED" refers to whatever name Porechop_ABI gives the file post-demultiplexing/adapter trimming (e.g. BC##).

#		After, I suggest identifying which contigs in assembly.fasta correspond to which chromosomes, then rename the fasta
#		headers for ease of identification in the following steps. Use Bandage for visualizing the .gfa assembly to get
#		the best visual reference of each assembly run.

#		Upload the file with chromosome names re-added as
		[NAME]-iden_assembly.fasta
		
#		***Check for telomeric sequences - I've noticed a few contigs that hit telomeric sequence, then fuse something else 
#		   on AFTER. Trim those (and the improper extra sequence) out of the resulting .fasta files in the chromosome-level 
#		   assembly step coming up!***
#			   TRUE 23-bp telomerase repeat (correct order for looping, and with 5'-3' order relative to reference genome)
				5'-ACACCAAGAAGTTAGACATCCGT-3' (left arm)
				5'-ACGGATGTCTAACTTCTTGGTGT-3' (right arm)

#			SPECIAL NOTES:
#				*For visualization of the assembly, load the "assembly_graph.gfa" into "Bandage" (https://rrwick.github.io/Bandage/),
#				 then check to see how the assembly looks and make sure it's GENERALLY close to Assembly 21/22 sizes.*
#				**Sometimes, you might need to run it a few times to get a nice looking map, or to have everything present.**
#				***With ~400-500X (haploid), my runtime for this was 1.5 hours. If you have deeper coverage, you might need to use the
#				   additional arguments,
				   --genome-size 14.9m --asm-coverage X 
#				   where X = reduced coverage for initial disjointig assembly - prioritizes longest reads.***
#				****With Q20 data, we theoretically could be adding "--read-error 0.03" to the command in order to further
#				    improve generation of the assemblies, but I never tested this.****
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 & SeqKit v2.5.0 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast &
		https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit

#		The following will grab the Flye assembly file, then run BUSCO and QUAST for calculating the quality of the assembly (QUAST
#		can use a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz"). You
#		must provide QUAST the reference if you want it to use a reference. You can also provide it genomic feature files (like .txt  
#		or .gff) for detecting those! Here, I've split the Ca22 assembly so that there's only the Ca22 A homologs being used for QUAST.
#		It also drops ChrM.
#			For the wget lines, if you can't find them at the provided links, delete the file name from the link and check the webpage
#			to see if there's an updated version of the assembly. If you definitively want to use the same version as we did, go up one
#			directory, go into the "archive" directory, then find the corresponding files, and update the links. This should also work
#			for ChrM further below during the DeconSeq step.
#				Additionally, if you aren't using Candida albicans (or even SC5314), you might need to upload your own files
#				at this point, then continue.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/
		mkdir ./busco_quast/
		cd /[STRAINPATH]/busco_quast/
		busco -i ../2_flye_assembly/[NAME]-iden_assembly.fasta -l saccharomycetes_odb10 -o ./2_busco_flye_results -m genome -c 28
		conda deactivate
		conda activate seqkit_env
		cd /[STRAINPATH]/busco_quast/
		wget 'http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/archive/C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz'
		gzip -d C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz
		seqkit grep -nrp 'Ca22chr[1234567R]A_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -w 70 -n > Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta
		seqkit grep -nrp 'Ca22chr[1234567R]B_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -w 70 -n > Ca_SC5314_A22-s07-m01-r187_chrs-B.fasta		
		seqkit grep -nrp 'Ca22chrM_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -w 70 -n  > Ca_SC5314_A22-s07-m01-r187_chrs-ChrM.fasta
		wget 'http://www.candidagenome.org/download/gff/C_albicans_SC5314/archive/C_albicans_SC5314_version_A22-s07-m01-r187_features.gff'
		grep -E '^[#*]|^Ca22.*A_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-A.gff
		grep -E '^[#*]|^Ca22.*B_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-B.gff
		grep -E '^[#*]|^Ca22.*M_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-ChrM.gff
		conda deactivate
		conda activate quast_env
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../2_flye_assembly/[NAME]-iden_assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta -g Ca_SC5314_A22-s07-m01-r187_feats-A.gff -o ./2_quast_flye_results -t 28
}


{Mitochondrial DNA Removal
	DeconSeq v0.4.3 # RAN ON OSC
		https://github.com/kiwiroy/deconseq | https://deconseq.sourceforge.net/

#		This requires slightly more complicated setup! I followed these instructions (with some modifications): https://web.archive.org/web/20220715225827/https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/deconseq.htm
#		Run the script as follows in a terminal to properly install and set up DeconSeq with Candida albicans ChrM 
#		(that you must provide - I fed it Ca22 ChrM from the A22-s07-m01-r187 chromosomes .fasta). You may need
#		to upload your own files if it's not Candida albicans (or even SC5314).
#			***EMPHASIS: this is OUTSIDE of a job!! Do this in a terminal.***
#			   Additionally, this starting part only needs to be run ONCE for all strains, due to /[HOMEPATH]/.
		mamba create -n deconseq_env -c bioconda -c conda-forge perl-cpan-shell perl
		conda activate deconseq_env
		cd /[HOMEPATH]/
		cpan Data::Dumper Getopt::Long Pod::Usage File::Path Cwd FindBin # If it asks you to configure automatically, say yes.
		wget 'https://sourceforge.net/projects/deconseq/files/standalone/deconseq-standalone-0.4.3.tar.gz'
		tar -zxvf deconseq-standalone-0.4.3.tar.gz
		ln -s deconseq-standalone-0.4.3 deconseq
		cd deconseq
		chmod +x *.pl
		chmod +x bwa64
		rm bwaMAC
		rm bwasw_modified_source/ -rf
		mkdir ./databases/
		nano DeconSeqConfig.pm # From here, edit the sections as follows:
			use constant DB_DIR => '/[HOMEPATH]/deconseq/databases/';
			use constant PROG_DIR => '/[HOMEPATH]/deconseq/';
			use constant DBS => {camito => {name => 'Candida albicans - Ca22 ChrM (mitochondrial DNA)', #database name used for display and used as input for -dbs and -dbs_retain
											db =>   'ca22_ChrM'}};                                      #database name as defined with -p for "bwa index -p ..." (specify multiple database chunks separated with commas without space; e.g. hs_ref_s1,hs_ref_s2,hs_ref_s3)
			use constant DB_DEFAULT => 'camito';
# 			After making the edits, hit "CTRL+O", then hit "ENTER" to confirm the changes.
		cd /[HOMEPATH]/
		cp .bash_profile bash_profile-BACKUP # This is just in case something goes very, very wrong.
		nano .bash_profile # From here, add to the very end of the file (right below "fi"):
			export PATH=$PATH:/[HOMEPATH]/deconseq
# 			After making the edits, hit "CTRL+O", then hit "ENTER" to confirm the changes.
		cd /[HOMEPATH]/deconseq/databases/
		bwa64 index -p ca22_ChrM -a is /[STRAINPATH]/busco_quast/Ca_SC5314_A22-s07-m01-r187_chrs-ChrM.fasta >out.txt 2>&1 &
#		There should now be a group of 8 files (technically 10 if you count the .fasta and out.txt)

#		***SETUP COMPLETE! For any future runs, just run the script below as a job (as per usual).***

#		The following will run deconseq and remove mitochondrial DNA sequences from the Flye assembly (back to jobs! It will
#		generate two files - one with your non-mitochondrial DNA ([NAME]_clean.fa) and one with your mito DNA ([NAME]_cont.fa).
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate deconseq_env
		cd /[STRAINPATH]/
		mkdir ./3_deconseq/
		cd /[STRAINPATH]/3_deconseq/
		perl /[HOMEPATH]/deconseq/deconseq.pl -f ../2_flye_assembly/[NAME]-iden_assembly.fasta -dbs camito -out_dir ./ -id [NAME]
}


{Chromosome-Level Assembling & Checking
	CGD BLASTN & SnapGene Viewer & IGV

#		Use a combination of them for identifying overlaps/places where you can combine chromosome contigs so that you can
#		make a full-chromosome contig (as best as you can, working around any inversions/issues observed).
		Download the [STRAINPATH]/3_deconseq/[NAME]_clean.fa assembly, then finish assembly by hand.
		Try to also remove any contigs that you are CERTAIN belong to unincorporated haplotypes.

#			You can use (in the conda seqkit_env)
			seqkit seq -w 0 [NAME]_clean.fa > [NAME]-[something].fasta
#			to linearize your sequences (since SnapGene won't let you copy more than 1,000,000 bp at a time). 

#			Then, after stitching them together, use
			seqkit seq -w 70 [NAME]-[something].fasta > [NAME]-assembled.fasta
#			to reintroduce line breaks, with 70 characters per line.

#		This step is very user-dependent. A lot of it is assemble by hand, then re-align long-reads to it and see whether things
#		generally agree or not. You just want ONE chromosome homolog, and so long as it is GENERALLY rightish, the later steps
#		SHOULD work work for phasing it into the two homologs.

#		After making your haploid contigs and removing any contigs you KNOW are likely haplotypes, upload your new 
		[NAME]-assembled.fasta
#		assembly into a new "/[STRAINPATH]/4_manualassembly/" directory.
		cd /[STRAINPATH]/
		mkdir ./4_manualassembly/

#		After uploading, run the next section's BUSCO/QUAST quality checks, as well as a minimap2 long-read alignment to the assembly to 
#		make sure there's nothing that completely disagrees. There'll be SOME disagreement, but when visualizing the long-read alignment 
#		in IGV you should clearly see at mostly normal coverage for joint regions between contig-joined parts. If you see a complete
#		lack of coverage, or odd coverage, this means it probably needs a different join, and you'll need to try a new arrangement
#		of contigs (or that the region is heterozygous or even heterozygous and inverted). 
#			For example, adjacent to the rDNA locus or the MRS loci there can sometimes be oddities with how Flye put things
#			together, because it'll try to put parts of the internal repeats at the border of the locus, causing some misalignments
#			for long reads. Fixing those is fairly easy, thankfully, as the misalignment should be obvious in IGV. After, run another
#			minimap2 alignment to see if that fits a bit better.

#			QUAST results from the manual assembly's contigs are very helpful for seeing how the contigs stack together!

#		Something else that can be helpful is running the long reads against a partially assembled chromosome, then opening the
#		alignment in IGV and setting it so you can see soft-clipped bases, then searching for those unaligned nucleotides to
#		figure out which contig should be next. You can even do this at the END of contigs in IGV to see if there are any little
#		gaps, and base-crawl your way to proper joinings (or hope that Medaka or Pilon can fix it up after).
#			You can also copy paste full read sequences into SnapGene. It'll be a touch messy, but it'll be enough to get it done.

#		For long-read alignment using minimap2 against your manual assembly, run the following:
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/4_manualassembly/
		flye-minimap2 -ax map-ont -t 28 [NAME]-assembled.fasta ../1_demul_adtrim/[PORECHOPPED]-10kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-assembled-10kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-assembled-10kbmin-mapping.bam
			
#		Then visualize in IGV.
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the haploid chromosome-level assembly file, then run BUSCO and QUAST for calculating the quality of the assembly
#		(QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz"). You must
#		provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for detecting those!
#			***NOTE: If you followed the installation steps above for the references, you shouldn't need to repeat the installation
#			   steps here or for any of the future BUSCO/QUAST steps.***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../4_manaulassembly/[NAME]-assembled.fasta -l saccharomycetes_odb10 -o ./4_busco_manual_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../4_manualassembly/[NAME]-assembled.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta -g Ca_SC5314_A22-s07-m01-r187_feats-A.gff -o ./4_quast_manual_results -t 28
}


{Polishing/ONT Long-Read Error Correction
	Medaka v1.8.0 # RAN ON OSC
		https://github.com/nanoporetech/medaka | https://anaconda.org/bioconda/medaka

#		The following will run medaka & subtools on the assembly, using the pruned basecalls ([PORECHOPPED]-25kb.fastq) for polishing.
#		The settings listed after "-m" are based on their suggestions for models.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate medaka_env
		cd /[STRAINPATH]/4_manualassembly/
		medaka_consensus -i ../1_demul_adtrim/[PORECHOPPED]-25kbmin.fastq -d ./[NAME]-assembled.fasta -o ../5_medaka_consensus -t 28 -m r1041_e82_400bps_hac_g632
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the medaka consensus assembly file, then run BUSCO and QUAST for calculating the quality of the
#		assembly (QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz").
#		You must provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for 
#		detecting those! 
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../5_medaka_consensus/consensus.fasta -l saccharomycetes_odb10 -o ./5_busco_medaka_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../5_medaka_consensus/consensus.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta -g Ca_SC5314_A22-s07-m01-r187_feats-A.gff -o ./5_quast_medaka_results -t 28
}


{Polishing/Illumina Short-Read Error Correction
	Pilon v1.24 & Bowtie2 v2.5.1 & Samtools v1.6 # RAN ON OSC
		https://github.com/broadinstitute/pilon/wiki | https://anaconda.org/bioconda/pilon

#		First, create an environment for Pilon (this requires multiple tools up front) via the terminal:
		mamba create -n pilon_env -c bioconda pilon bowtie2 samtools

#		***SETUP COMPLETE! For any future runs, just run the script below as a job (as per usual).***

#		The following will use Bowtie2 to index the medaka consensus assembly, then use Bowtie2 to align our trimmed Illumina
#		2x150 paired-end NextSeq2000 reads to the index, use samtools to convert the .sam into .bam, sort the .bam, then generate 
#		a .bai index of the sorted .bam file. Then, Pilon will call the consensus sequence, then use the pre-aligned short reads 
#		to further polish the sequence.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate pilon_env
		cd /[STRAINPATH]/
		mkdir ./6_pilon/
		cd /[STRAINPATH]/6_pilon/
		bowtie2-build ../5_medaka_consensus/consensus.fasta consensus_ref -p 28
		bowtie2 -3 1 -x consensus_ref -1 ../short_read/[NAME]_R1_paired.fq.gz -2 ../short_read/[NAME]_R2_paired.fq.gz -S [NAME]-short_vs_consensus.sam -p 28
		samtools view -bS -@ 28 [NAME]-short_vs_consensus.sam > [NAME]-short_vs_consensus.bam
		samtools sort -@ 28 [NAME]-short_vs_consensus.bam -o [NAME]-short_vs_consensus.sorted.bam
		samtools index -@ 28 [NAME]-short_vs_consensus.sorted.bam [NAME]-short_vs_consensus.sorted.bam.bai
		cd /[STRAINPATH]/6_pilon/
		java -Xmx120G -jar /[HOMEPATH]/miniconda3/envs/pilon_env/share/pilon-1.24-0/pilon.jar --genome ../5_medaka_consensus/consensus.fasta --changes --frags [NAME]-short_vs_consensus.sorted.bam --output [NAME]-r[#] --outdir ./pilon_r[#]

#		NOTE: check the .changes file after to see what was corrected. The ideal method (apparently) is to run Pilon until it
#		stops making changes. To do this, map Illumina reads against the Pilon-polished consensus, run Pilon again, check changes,
#		and repeat until it looks like Pilon is making minimal changes. r[#] stands for the polishing round - increase for each polishing 
#		run. r[#+1] means increase that number by 1 for tracking.
#			*Follow up: 7 rounds of polishing seemed to work well. Checking the r7 *.changes file, it looks like most of those remaining 
#			 changes are flipping back and forth between possible heterozygous positions, maybe.*
#			Example:
			source /[HOMEPATH]/miniconda3/bin/activate
			conda activate pilon_env
			cd /[STRAINPATH]/6_pilon/
			bowtie2-build ./pilon_r[#]/[NAME]-r[#].fasta [NAME]-r[#]_ref -p 28
			bowtie2 -3 1 -x [NAME]-r[#]_ref -1 ../short_read/[NAME]_R1_paired.fq.gz -2 ../short_read/[NAME]_R2_paired.fq.gz -S [NAME]-short_vs_r[#].sam -p 28
			samtools view -bS -@ 28 [NAME]-short_vs_r[#].sam > [NAME]-short_vs_r[#].bam
			samtools sort -@ 28 [NAME]-short_vs_r[#].bam -o [NAME]-short_vs_r[#].sorted.bam
			samtools index -@ 28 [NAME]-short_vs_r[#].sorted.bam [NAME]-short_vs_r[#].sorted.bam.bai
			cd /[STRAINPATH]/6_pilon/
			java -Xmx120G -jar /[HOMEPATH]/miniconda3/envs/pilon_env/share/pilon-1.24-0/pilon.jar --genome ./pilon_r[#]/[NAME]-r[#].fasta --changes --frags [NAME]-short_vs_r[#].sorted.bam --output [NAME]-r[#+1] --outdir ./pilon_r[#+1]


#		***WARNING - Pilon might not properly fix homopolymers above 4 bp! Some conflicting information online, but it's possible.
#		   It seems to be related to a --nanopore option that isn't explicitly noted anywhere... If there's a better short-read 
#		   corrector than Pilon, let me know!***
#		****It might actually be related to long reads, so short-read correction should theoretically work fine...?****
#			https://github.com/broadinstitute/pilon/releases/tag/v1.23 (this is the source of my information)
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the Pilon polished assembly file, then run BUSCO and QUAST for calculating the quality of the assembly
#		(QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz"). You must
#		provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for detecting those! 
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../6_pilon/pilon_r[#]/[NAME]-r[#].fasta -l saccharomycetes_odb10 -o ./6_busco_pilon_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../6_pilon/pilon_r[#]/[NAME]-r[#].fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta -g Ca_SC5314_A22-s07-m01-r187_feats-A.gff -o ./6_quast_pilon_results -t 28
}


{Long-Read Alignment - Preparing for Variant Identification/Phasing
	SeqKit v2.5.0 & Minimap2 v2.24 + SAMtools v1.9 (both contained within Flye v2.9.2) # RAN ON OSC
		https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit & https://github.com/fenderglass/Flye | https://anaconda.org/bioconda/flye

#		The following will run seqkit to alphanumerically organize the chromosome-level contigs (assuming you named them already),
#		then run minimap2, mapping the [PORECHOPPED]-10kbmin long reads against the polished assembly.
#			***NOTE: for the 'sed' step, replace the number of pilons with the number you have in your last polished file.***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate seqkit_env
		cd /[STRAINPATH]/6_pilon/pilon_r[#]/
		sed 's/pilon_pilon_pilon_pilon_pilon_pilon_pilon/pilon-r7/' [NAME]-r[#].fasta > [NAME]-r[#]-rename.fasta
		seqkit sort -w 70 -n [NAME]-r[#]-rename.fasta > [NAME]-collapsed-r[#].fasta
		conda deactivate
		conda activate flye_env
		cd /[STRAINPATH]/6_pilon/pilon_r[#]/
		flye-minimap2 -ax map-ont -t 28 [NAME]-collapsed-r[#].fasta ../../1_demul_adtrim/[PORECHOPPED]-10kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-collapsed-r[#]-10kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-collapsed-r[#]-10kbmin-mapping.bam

#		Afterwards, check the regions you fixed previously (i.e. flanking rDNA sequences or MRS sequences) to make sure that they
#		didn't get """fixed""" during the polishing (i.e. reverted back to the less-harmonious sequence). Correct them in the
#		[NAME]-collapsed-r[#].fasta file if necessary - just remember that if you DO correct anything, be sure to run the alignment
#		again on the new .fasta assembly so that the .bam files are aligned to the correct assembly.
}


{SNP Calling - Obtaining Variant Positions
	Clair3 v1.0.4 & LongPhase v1.5.1 # RAN ON OSC
		https://github.com/HKU-BAL/Clair3 | https://anaconda.org/bioconda/clair3 & https://github.com/twolinin/longphase

#		Prior to running, follow the installation guidelines as described under "Option 4. Build an anaconda virtual environment".
#		I especially suggest using mamba for this - speeds things up a lot. After, install "Rerio" (https://github.com/nanoporetech/rerio),
#		and use it to get the specific basecaller model used for clair3. In my case, it was "r1041_e82_400bps_hac_g632". Then, 
#		move that into the Clair3 models folder.

# 		After that, install longphase v1.5.1 (Clair3 only has v1.3) into your "/HOMEPATH/" directory as directed in the github.

#		The following will run Claire3 & subtools on the assembly, then call where it finds variants. Then, it'll use whatshap and
#		longphase to phase the called variants (indels will not be phased, as ONT has issues with homopolymers, which often end up
#		getting called as indels). It will still HAVE the called indels and phase them, it just won't ENABLE the phases.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate clair3_env
		cd /[STRAINPATH]/
		mkdir ./7_clair3-phasing/
		cd /[STRAINPATH]/6_pilon/pilon_r[#]/
		rm ./[NAME]-collapsed-r[#].fasta.fai
		samtools faidx [NAME]-collapsed-r[#].fasta
		/[HOMEPATH]/Clair3/run_clair3.sh --bam_fn=[NAME]-collapsed-r[#]-10kbmin-mapping.bam --ref_fn=[NAME]-collapsed-r[#].fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="ont" --model_path=/[HOMEPATH]/Clair3/models/r1041_e82_400bps_hac_g632/ --output=/[STRAINPATH]/7_clair3-phasing/
		cd /[STRAINPATH]/7_clair3-phasing/
		/[HOMEPATH]/longphase/longphase phase -s merge_output.vcf.gz -b ../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#]-10kbmin-mapping.bam -r ../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#].fasta -t 8 -o longphase-snp-phased_merge_output --ont
		bgzip -k longphase-snp-phased_merge_output.vcf && tabix longphase-snp-phased_merge_output.vcf.gz
		whatshap phase -o whatshap-snp-phased_merge_output.vcf --reference=../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#].fasta merge_output.vcf.gz ../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#]-10kbmin-mapping.bam --tag=PS --ignore-read-groups
		bgzip -k whatshap-snp-phased_merge_output.vcf && tabix whatshap-snp-phased_merge_output.vcf.gz
		
#		Afterwards, load the "[NAME]-collapsed-r[#].fasta" file, the "[NAME]-collapsed-r[#]-10kbmin-mapping.bam" file (with the ".bai"
#		in the same folder), and the various vcf files into IGV to verify that it has generally called most variants accurately.
#		The settings here SHOULD get most everything (except for larger indels/complex regions), but there's always a chance
#		it could miss something. 

#		Then, go to any large homozygous region (or region you particularly care about) in IGV, and try to identify variants from the
#		.vcf files that look real (i.e. 20+ quality score). If you can find one, try to figure out which allele it belongs to, then
#		manually annotate that in a copy of the .vcf file. This will help things down the line properly pull down reads! Additionally,
#		try to use some LOH'd strains to figure out alleles that fall FAR outside the region where phasing can be successfully done (i.e.
#		for SC5314, the far right arm of Chr3 or the far arms of Chr7 and ChrR).

#		I did this by loading up the whatshap .vcf in Excel, pasting the longphase .vcf "SAMPLE" column next to the whatshap one,
#		then going through the homozygous or homozygous-adjacent regions and making sure the two agree, and changing anything that
#		looked wrong or that looked real in the file/IGV. By changing, I mean turning "/" (off) to "|" (on) if I thought it was real (and 
#		assigning it to the same phase block as the rest of the chromosome if it was different), or turning it from "|" to "/" if it seemed 
#		like it wasn't real in IGV. If you can figure out WHICH homolog the variants out past homozygous regions belong to (i.e. near the ends),
#		make sure that the "#|#" is assigned correctly, too! The order tells programs which haplotype the REF|ALT belong to, and if it had
#		put them in their own phaseblock, it might have the wrong order for the OTHER phaseblock. Afterwards, I copied the information into
#		a new .vcf file (a copy of the original), excluding the extra "SAMPLE" column from longphase. 
#			*Don't forget to either add (turning on phasing) or remove (turning off phasing) the ":PS" in the "FORMAT" column, too! Plus the
#			 associated phase SET (the set of 5 numbers after the final colon, like "0|1:12:257:131,122:0.4747:39523"), where :39523 = PS.
#			**For my stuff, it looked like REAL INDELS usually had a quality greater than 20, and real SNPs had a quality greater
#			  than 5 or 6. HOWEVER, some INDELS that are real do fall below the 20 quality - I've just found that if the indel is
#			  above 20, it is almost guaranteed to be real and not to be a homopolyer/repetitive issue.
#			***Lastly, if you find any regions where there is clear heterozygosity (and probably some sort of inversion) and you
#			   still see a large chunk of reads not assigned to the right haplotype, it might be calling the reads incorrectly 
#			   because the reads inside the inversion begin getting mixed with the reads outside the inversion (you'll usually
#			   see a "wall" if this happens). You might be able to turn on/turn off called SNPs in that region to get better 
#			   assignments. For example, TCA4-4 on Chr4.

#		Upload the new file into the "/7_clair3-phasing/" directory with the name, 
		dual-man-phased_merge_output.vcf
}



At this point, you should have a PRETTY good collapsed assembly that can be used for most things, including knowing where MOST of the
smaller variants are in the strain relative to your collapsed reference.

Beyond this point is the process for generating a phased assembly.



{Haplotagging Reads & Read Splitting - Assigning Haplotypes to Each Long Read & Splitting Original Reads Based on Haplotype Assignment
	Clair3 v1.0.4 # RAN ON OSC
		https://github.com/HKU-BAL/Clair3 | https://anaconda.org/bioconda/clair3

#		The following will run whatshap on the collapsed assembly to tag each long read inside "[NAME]-collapsed-r[#]-10kbmin-mapping.bam" 
#		with a haplotype, as well as make a LIST of which reads are which haplotype. This list will then be used to split the original unmapped
#		reads into two entirely separate files called "hap1" and "hap2" (thus avoiding any issues with clipping), allowing FULL reassembly 
#		of the homologous chromosomes in the next steps! It will also throw all unassigned reads into BOTH phased hap files, ensuring 
#		that homozygous regions still get properly represented during assembly. There MAY be some reads that are from the other
#		haplotype, but that shouldn't be TOO big of an issue (and it should be very few).
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate clair3_env
		cd /[STRAINPATH]/7_clair3-phasing/
		bgzip -k dual-man-phased_merge_output.vcf && tabix dual-man-phased_merge_output.vcf.gz
		cd /[STRAINPATH]/
		mkdir ./8_clair3-read_split/
		cd /[STRAINPATH]/8_clair3-read_split/
		whatshap haplotag -o [NAME]-haplotagged.bam --reference ../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#].fasta ../7_clair3-phasing/dual-man-phased_merge_output.vcf.gz ../6_pilon/pilon_r[#]/[NAME]-collapsed-r[#]-10kbmin-mapping.bam --output-threads 28 --ignore-read-groups --output-haplotag-list [NAME]-haplotagged_list.tsv
		samtools index -@ 28 [NAME]-haplotagged.bam
		whatshap split --output-h1 [NAME]-hap1_reads.fastq --output-h2 [NAME]-hap2_reads.fastq --add-untagged ../1_demul_adtrim/[PORECHOPPED]-10kbmin.fastq [NAME]-haplotagged_list.tsv

#		Afterwards, load the "[NAME]-collapsed-r[#].fasta" file, the "[NAME]-haplotagged.bam" file (with the ".bai" in same folder), and the 
#		"dual-man-phased_merge_output.vcf.gz" file (with the ".tbi" in the same folder) into IGV to verify that it has generally called 
#		most variants accurately. The settings here SHOULD get most everything (except for unphased indels/complex regions), but there's
#		always a chance it could miss something.

#		YOU CAN NOW VISUALIZE WHICH READS HAVE BEEN GIVEN WHICH HAPLOTYPE TAG! Right click the track in IGV, do "Color reads by" and select 
#		"tag" and input "HP", then do "Group reads by" and select "phase". You should see semi-even coverage between the haplotypes 
#		for each group of reads, but sometimes one or the other might get more simply because this assembly isn't perfect.
#			***MOST IMPORTANTLY, MAKE SURE IT DIDN'T GET GREEDY WITH ASSIGNMENTS!! Look in the heterozygous regions flanking large homozygous
#			   regions - you should see generally an equal balance of each tag! Additionally, make sure you have the "add-untagged" setting
#			   turned on, otherwise the right arm of Chr3 and similar large homozygous regions will NOT have any sequence (it by default drops
#			   anything that hasn't been tagged, which is, again, a problem for homozygous regions).***
#		Additionally, CHECK CAREFULLY to make sure that you don't have any oddities in haplotype switching! Make sure reads agree ALL THE way
#		across the genome, and that the phased SNPs don't appear to swap at some point!
}


{Genome Reassembly - Phased
	Flye v2.9.2 & SeqKit v2.5.0 # RAN ON OSC
		https://github.com/fenderglass/Flye | https://anaconda.org/bioconda/flye & https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit

#		The following will run Flye using assumptions for ONT data basecalled with Guppy5+ SUP, grab the [NAME]-hap[1/2]_reads.fastq
#		files made in the last step, then run Flye. After running, it places the resulting directory  up into the folder of the strain
#		being processed for  ease of access in the next step. Probably allocate at least 4-5 hours in a job for this! Additionally, 
#		it ends up using minimum-length 25kb reads to minimize risk of messing things up, but it still will keep the 10kbmin split .fastq 
#		files for use in future steps/alignments.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate seqkit_env
		cd /[STRAINPATH]/
		mkdir ./9_flye_assembly-phased/
		cd /[STRAINPATH]/8_clair3-read_split/
		seqkit seq [NAME]-hap1_reads.fastq -m 25000 > [NAME]-hap1_reads-25kbmin.fastq
		seqkit seq [NAME]-hap2_reads.fastq -m 25000 > [NAME]-hap2_reads-25kbmin.fastq
		conda deactivate
		conda activate flye_env
		flye --nano-hq [NAME]-hap1_reads-25kbmin.fastq --out-dir ../9_flye_assembly-phased/flye_hap1 --threads 28
		flye --nano-hq [NAME]-hap2_reads-25kbmin.fastq --out-dir ../9_flye_assembly-phased/flye_hap2 --threads 28

#		After, I suggest identifying which contigs in assembly.fasta correspond to which chromosomes, then rename the fasta
#		headers for ease of identification in the following steps. Additionally, at this point you should be able to use an
#		older semi-phased assembly (e.g., Ca22 for us) to get a general idea of which contig from hap1 and hap2 are which 
#		homolog! i.e., sort out which contigs from hap1 and hap2 are Chr1, Chr2, etc, AND which one is the A vs B homolog.
#		Reminder - A and B will probably be mixed between the two hap files. Don't unmix them - the split reads need to be
#		associated with the same pool of contigs!
#			NOTE: I think it actually should consistently split the reads. From what I can see right now, for my .vcf...
#				hap1 = Chr1B, Chr2A, Chr3A, Chr4B, Chr5B, Chr6B, Chr7B, ChrRB
#				hap2 = Chr1A, Chr2B, Chr3B, Chr4A, Chr5A, Chr6A, Chr7A, ChrRA
#					Note0: These are specifically for MY .vcf file! Will probably differ if you make your own.
#					Note1: You could theoretically make it so that all of the assignments are corrected to always phase the
#					variants into the same split... Hm.
#					Note2: The telomere-proximal portion of Chr3R and Chr7L are currently randomly assigned, as the variants
#					(SNPs for Chr3R, LEU2+/- for Chr7L) are too far away from any other variants to assign a haplotype, and
#					we have no strains for Chr7L that could be used to determine which homolog has which variant. BWP17 (MAY307)
#					could be used for Chr3R, but I have not had time to do this yet.

#		Upload the file with chromosome names and homologs re-added as
		[NAME]-hap[1/2]-iden_assembly.fasta

#		***NOTE2: Check for telomeric sequences - I've noticed a few contigs that hit telomeric sequence, then fuse something else 
#		   on AFTER. Trim those (and the improper extra sequence) out of the .fasta in the chromosome-level assembly step coming up!***
#			   TRUE 23-bp telomerase repeat (correct order for looping, and with 5'-3' order relative to reference genome)
				5'-ACACCAAGAAGTTAGACATCCGT-3' (left arm)
				5'-ACGGATGTCTAACTTCTTGGTGT-3' (right arm)

#			SPECIAL NOTES:
#				*For visualization of the assembly, load the "assembly_graph.gfa" into "Bandage" (https://rrwick.github.io/Bandage/),
#				 then check to see how the assembly looks and make sure it's GENERALLY close to Assembly 21/22 sizes.*
#				**Sometimes, you might need to run it a few times to get a nice looking map, or to have everything present.**
}


{Quality Assessment - Phased
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the chromosome-level assembly files, then run BUSCO and QUAST for calculating the quality of the assembly
#		(QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz"). You must
#		provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for detecting those!
#			I'd suggest trying to pair each hap file with the assembly that the majority of the contigs agree with.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../9_flye_assembly-phased/flye_hap1/[NAME]_hap1-iden_assembly.fasta -l saccharomycetes_odb10 -o ./9_busco_flye-hap1_results -m genome -c 28
		busco -i ../9_flye_assembly-phased/flye_hap2/[NAME]_hap2-iden_assembly.fasta -l saccharomycetes_odb10 -o ./9_busco_flye-hap2_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../9_flye_assembly-phased/flye_hap1/[NAME]_hap1-iden_assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./9_quast_flye-hap1_results -t 28
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../9_flye_assembly-phased/flye_hap2/[NAME]_hap2-iden_assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./9_quast_flye-hap2_results -t 28}
}


{Mitochondrial DNA Removal - Phased
	DeconSeq v0.4.3 # RAN ON OSC
		https://github.com/kiwiroy/deconseq | https://deconseq.sourceforge.net/

#		The following will run deconseq and remove mitochondrial DNA sequences from the Flye assembly. It will generate two files - 
#		one with your non-mitochondrial DNA ([NAME]_hap[1/2]_clean.fa) and one with your mito DNA ([NAME]_hap[1/2]_cont.fa).
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate deconseq_env
		cd /[STRAINPATH]/
		mkdir ./10_deconseq-phased/
		cd /[STRAINPATH]/10_deconseq-phased/
		perl /[HOMEPATH]/deconseq/deconseq.pl -f ../9_flye_assembly-phased/flye_hap1/[NAME]_hap1-iden_assembly.fasta -dbs camito -out_dir ./ -id [NAME]_hap1
		perl /[HOMEPATH]/deconseq/deconseq.pl -f ../9_flye_assembly-phased/flye_hap2/[NAME]_hap2-iden_assembly.fasta -dbs camito -out_dir ./ -id [NAME]_hap2
}


{Chromosome-Level Assembling & Checking
	CGD BLASTN & SnapGene Viewer & IGV

#		Use a combination of them for identifying overlaps/places where you can combine chromosome contigs so that you can
#		make a full-chromosome contig (as best as you can, working around any inversions/issues observed).
		Download the /[STRAINPATH]/10_deconseq-phased/[NAME]_clean.fa assembly, then finish assembly by hand.
		Try to also remove any contigs that you are CERTAIN belong to unincorporated haplotypes.

#			You can use (in the conda seqkit_env)
			seqkit seq -w 0 [NAME]_hap[1/2]-clean.fa > [NAME]_hap[1/2]-[something].fasta
#			to linearize your sequences (since SnapGene won't let you copy more than 1,000,000 bp at a time). 

#			Then, after stitching them together, use
			seqkit seq -w 70 [NAME]_hap[1/2]-[something].fasta > [NAME]_hap[1/2]-assembled.fasta
#			to reintroduce line breaks, with 70 characters per line.

#		This step is very user-dependent. A lot of it is assemble by hand, then re-align long-reads to it and see whether things
#		generally agree or not. This will also probably require more work than the original haploid assembly, too, as it's very likely
#		some regions have been slightly missassembled due to carryover of reads that were either unmapped or mapped incorrectly, especially
#		from repetitive regions like the MRS loci. Be ESPECIALLY sure to look at the MRS or similar repetitive regions at this stage using
#		IGV and minimap2 mapped long reads. 
#			I found that while most everything was assembled correctly, Chr1A's MRS and ChrRA's MRS were both misassembled.

		USE THIS STEP TO ASSEMBLE THE TELOMERES! Use the original assembly files (especially those assembled with telomeres) and the 
		IGV alignments to try to determine the best version of the telomeres for each haplotype. You can even grab reads from IGV AND
		generate a consensus sequence. I personally included 3 repeats of the 23 bp sequence at each end.

#		After making your haploid contigs and removing any contigs you KNOW are orphans (i.e. small fragments), upload your new 
		[NAME]_hap[1/2]-assembled.fasta
#		assembly into a new "/[STRAINPATH]/11_manualassembly-phased/" directory.
		cd /[STRAINPATH]/
		mkdir ./11_manualassembly-phased/

#		After uploading, run the next section's BUSCO/QUAST quality checks, as well as a minimap2 long-read alignment to the assembly to 
#		make sure there's nothing that completely disagrees. There'll be SOME disagreement, but when visualizing the long-read alignment 
#		in IGV you should clearly see at mostly normal coverage for joint regions between contig-joined parts. If you see a complete
#		lack of coverage, or odd coverage, this means it probably needs a different join, and you'll need to try a new arrangement
#		of contigs (or that the region is heterozygous or even heterozygous and inverted). If you fix things, be sure to rerun your
#		quality checks (ESPECIALLY IGV with reads aligned to your new sequence) to make sure it actually fixed things.
#			For example, adjacent to the rDNA locus or the MRS loci there can sometimes be oddities with how Flye put things
#			together, because it'll try to put parts of the internal repeats at the border of the locus, causing some misalignments
#			for long reads. Fixing those is fairly easy, thankfully, as the misalignment should be obvious in IGV. After, run another
#			minimap2 alignment to see if that fits a bit better.

#			QUAST results from the manual assembly's contigs are very helpful for seeing how the contigs stack together!

#		Something else that can be VERY helpful is running the long reads against a partially assembled chromosome, then opening the
#		alignment in IGV and setting it so you can see soft-clipped bases, then searching for those unaligned nucleotides to
#		figure out which contig should be next. You can even do this at the END of contigs in IGV to see if there are any little
#		gaps, and base-crawl your way to proper joinings (or hope that Medaka or Pilon can fix it up after).
#			You can also copy paste full read sequences into SnapGene. It'll be a touch messy, but it'll be enough to get it done.

#		For long-read alignment using minimap2 against your manual assembly, run the following:
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/11_manualassembly-phased/
		flye-minimap2 -ax map-ont -t 28 [NAME]_hap1-assembled.fasta ../8_clair3-read_split/[NAME]-hap1_reads.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]_hap1-assembled-mapping.bam
		flye-samtools index -@ 28 [NAME]_hap1-assembled-mapping.bam
		flye-minimap2 -ax map-ont -t 28 [NAME]_hap2-assembled.fasta ../8_clair3-read_split/[NAME]-hap2_reads.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]_hap2-assembled-mapping.bam
		flye-samtools index -@ 28 [NAME]_hap2-assembled-mapping.bam

#		Then visualize in IGV.
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the chromosome-level assembly files, then run BUSCO and QUAST for calculating the quality of the assemblies
#		(QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz"). You must
#		provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for detecting those!
#			***NOTE: If you followed the installation steps above for the references, you shouldn't need to repeat the installation
#			   steps here or for any of the future BUSCO/QUAST steps.***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../11_manaulassembly-phased/[NAME]_hap1-assembled.fasta -l saccharomycetes_odb10 -o ./11_busco_manual_hap1_results -m genome -c 28
		busco -i ../11_manaulassembly-phased/[NAME]_hap2-assembled.fasta -l saccharomycetes_odb10 -o ./11_busco_manual_hap2_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../11_manualassembly-phased/[NAME]_hap1-assembled.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./11_quast_manual_hap1_results -t 28
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../11_manualassembly-phased/[NAME]_hap2-assembled.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./11_quast_manual_hap2_results -t 28
}


{Polishing/ONT Long-Read Error Correction & Manual Corrections
	Medaka v1.8.0 # RAN ON OSC
		https://github.com/nanoporetech/medaka | https://anaconda.org/bioconda/medaka

#		The following will run medaka & subtools on the assembly, using the 25+ kb phased read pools for polishing. The settings listed 
#		after "-m" are based on their suggestions for models. It will then resort the assemblies into alphanumerical order.
#			***BE SURE TO CHECK THE ASSEMBLY AFTERWARDS WITH MINIMAP2!!! It's very possible medaka may have "goofed" some of the MRS loci, 
#			   especially any that you had to fix during the manual assembly step. Perform long-read alignment with flye-minimap2 as 
#			   described above for mapping the phased long reads against your newly polished assemblies.***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate medaka_env
		cd /[STRAINPATH]
		mkdir ./12_medaka_consensus-phased/
		cd /[STRAINPATH]/11_manualassembly-phased/
		medaka_consensus -i ../8_clair3-read_split/[NAME]-hap1_reads-25kbmin.fastq -d ./[NAME]_hap1-assembled.fasta -o ../12_medaka_consensus-phased/hap1_consensus -t 28 -m r1041_e82_400bps_hac_g632
		medaka_consensus -i ../8_clair3-read_split/[NAME]-hap2_reads-25kbmin.fastq -d ./[NAME]_hap2-assembled.fasta -o ../12_medaka_consensus-phased/hap2_consensus -t 28 -m r1041_e82_400bps_hac_g632

#		In addition to realignment of long reads to your assemblies, a good FIRST pass is aligning your polished assemblies against your unpolished 
#		assemblies! This can let you REALLY quicky scan through them in IGV and see what all it corrected, THEN hop to the same place in your 
#		long-read IGV sessions (for both your unpolished AND polished versions) and see whether the reads seem to generally support the changes. 
#		Repeat as necessary for comparing the assemblies to each other AND the phased long reads to the new versions.
#			*NOTE1: I suggest only really doing a deep-dive if it's either a base change OR in a region you have previousy corrected. DON'T 
#			 manually fixing anything in HOMOPOLYMERS unless it is abdundantly clear it needs fixing. Medaka knows better than you for this one.*
#			**NOTE2: Don't be afraid to manually revert changes if the a region of the original assembly seems more correct against the long 
#			  reads! I've seen it "miscorrect" the MRS loci that I had manually entered before, due to aberrant reads aligning from OTHER MRS 
#			  loci. Don't dismiss ALL of the changes entirely though: for example, it DID correctly polish a small indel in Chr5B MRS, but 
#			  miscorrected a bunch of SNPs in other parts of the MRS.**
#			***NOTE3: Honestly, trust the long reads. If the haplophased long reads tell you the assembly should be something else, fix it to
#			   that. If you get weird bases that are heterozygous for BOTH "SNPs", sort by strand in IGV and see if it's a strand-specific SNP!***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/12_medaka_consensus-phased/
		flye-minimap2 -ax asm5 -t 28 --eqx	../11_manualassembly-phased/[NAME]_hap1-assembled.fasta ./hap1_consensus/consensus.fasta | flye-samtools sort -@ 28 -m 4G > [NAME]_hap1_pol-vs-1376_hap1_unpol.bam
		flye-samtools index -@ 28 [NAME]_hap1_pol-vs-1376_hap1_unpol.bam
		flye-minimap2 -ax asm5 -t 28 --eqx	../11_manualassembly-phased/[NAME]_hap2-assembled.fasta ./hap2_consensus/consensus.fasta | flye-samtools sort -@ 28 -m 4G > [NAME]_hap2_pol-vs-1376_hap2_unpol.bam
		flye-samtools index -@ 28 [NAME]_hap2_pol-vs-1376_hap2_unpol.bam

#		After you finish any remaining edits, sort your fasta file and make the final named file like this (rename the "[CONSENSUS].fasta" below to
#		whatever your actual final version of the unsorted consensus file is named):
		conda activate seqkit_env
		cd /[STRAINPATH]/12_medaka_consensus-phased/hap1_consensus/
		seqkit sort -w 70 -n [HAP1CONSENSUS].fasta > [NAME]_hap1-consensus.fasta
		cd /[STRAINPATH]/12_medaka_consensus-phased/hap2_consensus/
		seqkit sort -w 70 -n [HAP2CONSENSUS].fasta > [NAME]_hap2-consensus.fasta

#		At this point, you should sort the chromosomes out into two new .fasta files for their haplotype phasing! Do this by cutting and pasting within
#		the files themselves, then repackaging the chromosome sequences into the haplotype-organized versions. Keep the  original [...]-consensus.fasta 
#		files, just in case! You don't want to find your cutting tool/program you're using can't "remember" that much text at once.
#			Also, make sure each line only has 60-70 characters! The above seqkit sorting should do that automatically, but some tools like YMAP strongly 
#			prefer 60-70 characters per line, and may fail with 80+.
		C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta
		C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta

#		You should now realign your [10/25] kb+ [PORECHOPPED] and short read files against the properly haplotype-split assemblies, for validating that 
#		everything is nice and neatly sorted (and for any last checks). You can also use these to generate a new variant calling file, if you so desire.
#		Aligning the short reads can also help to finalize any questionable homopolymers, though be VERY careful in repetitive regions/homologous genes,
#		as there very likely will be short reads mismapping to the wrong homologs. As above, if you make any final edits, realign to be sure it improved.
#		I also suggest if you do this that you align it against BOTH homologs, just to be sure any heterozygous indels (especially large ones) are 
#		accurate for BOTH homologs. 
#			I personally would suggest the 25kbmin files, assuming you have high depth. It minimizes crosstalk between regions like high-similarity
#			subtelomeres.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/12_medaka_consensus-phased/
		mkdir ./alignments/
		cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/
		flye-minimap2 -ax map-ont -t 28 ../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta ../../1_demul_adtrim/[PORECHOPPED]-[10/25]kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-A_Chrs-[10/25]25kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-A_Chrs-[10/25]kbmin-mapping.bam
		flye-minimap2 -ax map-ont -t 28 ../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta ../../1_demul_adtrim/[PORECHOPPED]-[10/25]kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-B_Chrs-[10/25]25kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-B_Chrs-[10/25]kbmin-mapping.bam
		conda deactivate
		conda activate pilon_env
		mkdir ./refs/
		mkdir ./shortread/
		cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/refs/
		bowtie2-build ../../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta [NAME]_A-ref -p 28
		bowtie2-build ../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta [NAME]_B-ref -p 28
		cd ../shortread/
		bowtie2 -3 1 -x ../refs/[NAME]_A-ref -1 ../../../short_read/[NAME]_R1_paired.fq.gz -2 ../../../short_read/[NAME]_R2_paired.fq.gz -S [NAME]_short-vs-[NAME]_A.sam -p 28
		samtools view -bS -@ 28 [NAME]_short-vs-[NAME]_A.sam > [NAME]_short-vs-[NAME]_A.bam
		samtools sort -@ 28 [NAME]_short-vs-[NAME]_A.bam -o [NAME]_short-vs-[NAME]_A.sorted.bam
		samtools index -@ 28 [NAME]_short-vs-[NAME]_A.sorted.bam [NAME]_short-vs-[NAME]_A.sorted.bam.bai
		bowtie2 -3 1 -x ../refs/[NAME]_B-ref -1 ../../../short_read/[NAME]_R1_paired.fq.gz -2 ../../../short_read/[NAME]_R2_paired.fq.gz -S [NAME]_short-vs-[NAME]_B.sam -p 28
		samtools view -bS -@ 28 [NAME]_short-vs-[NAME]_B.sam > [NAME]_short-vs-[NAME]_B.bam
		samtools sort -@ 28 [NAME]_short-vs-[NAME]_B.bam -o [NAME]_short-vs-[NAME]_B.sorted.bam
		samtools index -@ 28 [NAME]_short-vs-[NAME]_B.sorted.bam [NAME]_short-vs-[NAME]_B.sorted.bam.bai
}


{Quality Assessment
	BUSCO v5.4.7 & QUAST v5.0.2 # RAN ON OSC
		https://busco.ezlab.org/ | https://anaconda.org/bioconda/busco & https://github.com/ablab/quast | https://anaconda.org/bioconda/quast

#		The following will grab the medaka consensus assembly file, then run BUSCO and QUAST for calculating the quality of the
#		assembly (QUAST can a reference genome - in this case, I used "C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz").
#		You must provide QUAST the reference if you want it. You can also provide it genomic feature files (like .txt or .gff) for 
#		detecting those! 
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate busco_env
		cd /[STRAINPATH]/busco_quast/
		busco -i ../12_medaka_consensus-phased/hap1_consensus/consensus.fasta -l saccharomycetes_odb10 -o ./12_busco_medaka_hap1_results -m genome -c 28
		busco -i ../12_medaka_consensus-phased/hap2_consensus/consensus.fasta -l saccharomycetes_odb10 -o ./12_busco_medaka_hap2_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /[STRAINPATH]/busco_quast/
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../12_medaka_consensus-phased/hap1_consensus/consensus.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./5_quast_medaka_results -t 28
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../12_medaka_consensus-phased/hap2_consensus/consensus.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./5_quast_medaka_results -t 28
}



{
-|-|-|-|-|-|-|-
'vIN PROGRESSv'
-|-|-|-|-|-|-|-

{Assembling Mitochondria!
source /[HOMEPATH]/miniconda3/bin/activate
conda activate flye_env
cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/
flye-minimap2 -ax map-ont -t 28 ../../busco_quast/Ca_SC5314_A22-s07-m01-r187_chrs-ChrM.fasta ../../1_demul_adtrim/[PORECHOPPED]-25kbmin.fastq | flye-samtools sort -@ 28 -m 4G > Ca22_ChrM-[NAME]_25kbmin-mapping.bam
flye-samtools index -@ 28 Ca22_ChrM-[NAME]_25kbmin-mapping.bam

# grab several reads that are 40+ kb and at least map over a good chunk of the Ca22_ChrM sequence, then compare the ChrM contigs
# from Flye against the reads to figure out the organization. NOTE: the mtDNA in C. albicans might actually be in multiple states,
# including NOT circular, but what it assembled is at least A unit of mtDNA.

# assemble by hand using the reads from the FLye assembly, using reads grabbed from the alignment to help guide it. then, realign
# reads to the new assembly and see how it looks! Name the file,
C_albicans-[NAME]-ChrM-denovo.fasta
# and place it in the same folder as the medaka-polished sequences. at least for C. albicans, no need to ACTUALLY polish with
# medaka, as there's nothing complex enough to truly "confuse" flye as far as i can tell. otherwise, follow the above instructions
# for polishing with medaka, adapted to your work!

source /[HOMEPATH]/miniconda3/bin/activate
conda activate flye_env
cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/
flye-minimap2 -ax map-ont -t 28 ../C_albicans-[NAME]-ChrM-denovo.fasta ../../1_demul_adtrim/BC15-25kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]_ChrM-25kbmin-mapping.bam
flye-samtools index -@ 28 [NAME]_ChrM-25kbmin-mapping.bam
cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/refs/
conda deactivate
conda activate pilon_env
bowtie2-build ../../C_albicans-[NAME]-ChrM-denovo.fasta [NAME]_ChrM-ref -p 28
cd ../shortread/
bowtie2 -3 1 -x ../refs/[NAME]_ChrM-ref -1 ../../../short_read/[NAME]_R1_paired.fq.gz -2 ../../../short_read/[NAME]_R2_paired.fq.gz -S [NAME]_short-vs-[NAME]_ChrM.sam -p 28
samtools view -bS -@ 28 [NAME]_short-vs-[NAME]_ChrM.sam > [NAME]_short-vs-[NAME]_ChrM.bam
samtools sort -@ 28 [NAME]_short-vs-[NAME]_ChrM.bam -o [NAME]_short-vs-[NAME]_ChrM.sorted.bam
samtools index -@ 28 [NAME]_short-vs-[NAME]_ChrM.sorted.bam [NAME]_short-vs-[NAME]_ChrM.sorted.bam.bai

# then check in IGV. if things disagree, edit and run the above again (delete the previously generated files first!)

# Finally, format it to be 70 characters per line like the other assembly files:
conda activate seqkit_env
cd /[STRAINPATH]/12_medaka_consensus-phased
seqkit sort -w 70 -n C_albicans-[NAME]-ChrM-denovo.fasta > C_albicans-[NAME]-ChrM-denovo_manual.fasta
}


{Polishing Homopolymers!
# used a combination of our reads, a recent PacBio dataset of SC5314, and short-reads to identify homopolymers that were likely
# incorrect, then use 1376-specific short reads to verify the regions for MAY1376

# start with calling variants using your long reads like last time, will use this to haplotype-split
source /[HOMEPATH]/miniconda3/bin/activate
conda activate clair3_env
cd /[STRAINPATH]/12_medaka_consensus-phased/
mkdir ./homo_resolve/
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/
mkdir ./presplit/
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/presplit/
/[HOMEPATH]/Clair3/run_clair3.sh --bam_fn=../../alignments/[NAME]-B_Chrs-25kbmin-mapping.bam --ref_fn=../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="ont" --model_path=/[HOMEPATH]/Clair3/models/r1041_e82_400bps_hac_g632/ --output=/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/presplit/
/[HOMEPATH]/longphase/longphase phase -s merge_output.vcf.gz -b ../../alignments/[NAME]-B_Chrs-25kbmin-mapping.bam -r ../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta -t 8 -o longphase-snp-phased_merge_output --ont
bgzip -k longphase-snp-phased_merge_output.vcf && tabix longphase-snp-phased_merge_output.vcf.gz
whatshap phase -o whatshap-snp-phased_merge_output.vcf --reference=../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta merge_output.vcf.gz ../../alignments/[NAME]-B_Chrs-25kbmin-mapping.bam --tag=PS --ignore-read-groups
bgzip -k whatshap-snp-phased_merge_output.vcf && tabix whatshap-snp-phased_merge_output.vcf.gz

# used B because the subtelomeres fortuitously always had a YRF1 gene if a YRF1 gene was at the end of that chromosome.
# load in whatshap in excel, paste longphase 'sample' next to whatshap 'sample', then check to see that they agree. you also NEED
# to swap any 1/0 and 1|0 to 0/1 and 0|1, respectively, because REFERENCE (B) is "0" in the GT, because it's the reference. therefore,
# if it's SWAPPED with 1 in the first position, it's saying the REFERENCE allele is "1", which will split the reads up without ALL of
# them being from the same haplotype for every chromosome!! tl;dr - just find a way to do a swap. should be able to do it pretty easily
# in excel. also, i found that longphase mucked up Chr5 by calling half of it as 0|1 (correct) and half of it as 1|0 (incorrect) - it
# got confused due to the MRS locus. this emphasizes the importance of making sure you correct the calls in excel!

# also, if you had previously enabled variant calls such as the SINGULAR SNP on the right arm of ChrR or some of the SNPs near the
# LEU2+/- locus in Chr7L, you should go find them again and ensure they're enabled! i.e. ":PS" tag and "|" to denote it's ON (/ = off)
# additionally, find your MRS loci and MAKE SURE IT ISN'T PULLING IN OTHER MRS LOCI AND GOOFING VARIANT CALLS UP! MRS-5 for me has
# repeatedly been problematic, and I'm probably going to need to disable a lot of the variant calls in it and manually add indels.
# surprisingly, the rest have been fine.

# after, paste it all into a new tab-delimited variant calling file named
dual-man-denophased_merge_output.vcf
# and use that for haplotagging and read splitting. upload it into the directory,
/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/presplit/


# below splits hoyer SC5314 pacbio reads: no reason to split OUR reads here because this is now working for HOMOPOLYMER polishing,
# and we already know that our ONT reads get iffy in long stretches, and that short reads are hard to be phased. so why not
# use longish, high quality SC5314 reads and then as mentioned above just use our short-read to VERIFY what the pacbio finds?

# pacbio reads came from data associated with https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/38362496/, which were then placed into
# the following folder 
/[STRAINPATH]/Hoyer_SC5314ATCC_PacBio/

# after, they were unzipped and then subset to have minimum 15 kb reads (made sure it still had at least 100-200X coverage).
# ran this in terminal command by ommand to actively check things
cd /[STRAINPATH]/Hoyer_SC5314ATCC_PacBio/
gunzip -c SRR23724250.fastq.gz > SRR23724250.fastq
conda activate seqkit_env
seqkit stats SRR23724250.fastq
seqkit seq SRR23724250.fastq -m 15000 > SC5314_ATCC-15kbmin.fastq
seqkit stats SC5314_ATCC-15kbmin.fastq

# then run the following for the mapping, then spot-checked in IGV to make sure it looked good alongside OUR long reads/agreed
# 	be sure to TURN OFF "hide small indels" - pacbio doesn't need it enabled! very accurate.
source /[HOMEPATH]/miniconda3/bin/activate
conda activate flye_env
cd /[STRAINPATH]/12_medaka_consensus-phased/alignments/
flye-minimap2 -ax map-hifi -t 28 ../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta ../../Hoyer_SC5314ATCC_PacBio/SC5314_ATCC-15kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-A_Chrs-Hoyer15kbmin_mapping.bam
flye-samtools index -@ 28 [NAME]-A_Chrs-Hoyer15kbmin_mapping.bam
flye-minimap2 -ax map-hifi -t 28 ../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta ../../Hoyer_SC5314ATCC_PacBio/SC5314_ATCC-15kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-B_Chrs-Hoyer15kbmin_mapping.bam
flye-samtools index -@ 28 [NAME]-B_Chrs-Hoyer15kbmin_mapping.bam

# after checking in IGV, actually run the splitting of B-aligned pacbio reads using the previously generated variant calling file
# again, specifically against B because it had all YRF1 copies
source /[HOMEPATH]/miniconda3/bin/activate
conda activate clair3_env
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/presplit/
bgzip -k dual-man-denophased_merge_output.vcf && tabix dual-man-denophased_merge_output.vcf.gz
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/
mkdir ./tagged/
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/tagged/
whatshap haplotag -o [NAME]-B_Chrs-Hoyer15kbmin_haplotagged.bam --reference ../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta ../presplit/dual-man-denophased_merge_output.vcf.gz ../../alignments/[NAME]-B_Chrs-Hoyer15kbmin_mapping.bam --output-threads 28 --ignore-read-groups --output-haplotag-list [NAME]-B_Chrs-Hoyer15kbmin_haplotagged_list.tsv
samtools index -@ 28 [NAME]-B_Chrs-Hoyer15kbmin_haplotagged.bam
whatshap split --output-h1 Hoyer15kbmin-B_reads.fastq --output-h2 Hoyer15kbmin-A_reads.fastq --add-untagged ../../../Hoyer_SC5314ATCC_PacBio/SC5314_ATCC-15kbmin.fastq [NAME]-B_Chrs-Hoyer15kbmin_haplotagged_list.tsv
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/
mkdir ./A_homs/
mkdir ./B_homs/
mv ./tagged/Hoyer15kbmin-A_reads.fastq ./A_homs/
mv ./tagged/Hoyer15kbmin-B_reads.fastq ./B_homs/

# after that, then it's time to realign the split reads and then call variants again! i decided to map both homopolymers AND general
# variants, just in case it can help me find any differences between SC5314 and MAY1376.
#		use short reads for SPECIFICALLY homopolymer-position calling, because otherwise it'll just call ALL variants...


#mamba create -n freebayes_env -c bioconda freebayes

# runs the stuff for realigning and calling new variants
source /[HOMEPATH]/miniconda3/bin/activate
conda activate flye_env
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/
flye-minimap2 -ax map-hifi -t 28 ../../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta ./Hoyer15kbmin-A_reads.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-A_Chrs-Hoyer15kbmin-A_reads-mapping.bam
flye-samtools index -@ 28 [NAME]-A_Chrs-Hoyer15kbmin-A_reads-mapping.bam
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/
flye-minimap2 -ax map-hifi -t 28 ../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta ./Hoyer15kbmin-B_reads.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-B_Chrs-Hoyer15kbmin-B_reads-mapping.bam
flye-samtools index -@ 28 [NAME]-B_Chrs-Hoyer15kbmin-B_reads-mapping.bam
conda deactivate
conda activate clair3_env
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/
mkdir ./varcalls/
#mkdir ./homocalls/
/[HOMEPATH]/Clair3/run_clair3.sh --bam_fn=./[NAME]-A_Chrs-Hoyer15kbmin-A_reads-mapping.bam --ref_fn=../../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="hifi" --model_path=/users/PAS1802/woodruff207/Clair3/models/hifi_sequel2/ --output=/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/varcalls/
cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/
mkdir ./varcalls/
#mkdir ./homocalls/
/[HOMEPATH]/Clair3/run_clair3.sh --bam_fn=./[NAME]-B_Chrs-Hoyer15kbmin-B_reads-mapping.bam --ref_fn=../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="hifi" --model_path=/users/PAS1802/woodruff207/Clair3/models/hifi_sequel2/ --output=/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/varcalls/
#conda deactivate
#conda activate freebayes_env
#cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/homocalls/
#freebayes 
#--bam_fn=./[NAME]-A_Chrs-Hoyer15kbmin-A_reads-mapping.bam --ref_fn=../../../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="ilmn" --model_path=/users/PAS1802/woodruff207/Clair3/models/ilmn/ --output=/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/homocalls/ --bed_fn=/[STRAINPATH]/12_medaka_consensus-phased/annotations/A_8Homo.bed
#cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/homocalls/
#freebayes 
#--bam_fn=./[NAME]-B_Chrs-Hoyer15kbmin-B_reads-mapping.bam --ref_fn=../../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta --threads=28 --include_all_ctgs --var_pct_full=1 --ref_pct_full=1 --var_pct_phasing=1 --platform="ilmn" --model_path=/users/PAS1802/woodruff207/Clair3/models/ilmn/ --output=/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/homocalls/ --bed_fn=/[STRAINPATH]/12_medaka_consensus-phased/annotations/B_8Homo.bed

# after all of the above, FINALLY figure out what all the variants are, if they're real based on short reads, and then 
# correct them and make it all a new assembly! delete the WRONG lines out of a NEW .vcf file before the following step, as well
# as any lines you didnt check as it could be an INCORRECT call!

# Name the new variant calling files 
"A_varcalls-man-merge_output.vcf" and "B_varcalls-man-merge_output.vcf" 
# and place them into the respective directories:
/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/varcalls/A_varcalls-man-merge_output.vcf
/[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/varcalls/B_varcalls-man-merge_output.vcf


USE BCFTOOLS CONSENSUS - YOU CAN USE IT TO MAKE AN ASSEMBLY WITH THE """HAPLOTYPE 2""" ALLELES (i.e. the variants you confirm) AND
THEN APPLY THEM TO THE SEQUENCE TO AUTOMATICALLY INCORPORATE THEM!!! THAT WOULD WORK PERFECTLY FOR WHAT IM DOING AND MEANS THERES
NO NEED TO APPLY THEM ALL BY HAND!!! Just use short-read to verify the SNPs, and just means i need to figure out how to check our
short read for any indels that werent present in the Hoyer reads!

#		The following will run BCFtools on the assembly, using the manually verified variant file (with unchecked or incorrect
#		variants either modified or deleted as described above) to generate a version of the assembly with the updated Variant
#		positions! Both will use "haplotype 2", as that is the haplotype possessing the corrections! It'll then save it as a 
#		new file.
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate clair3_env
		cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/varcalls/
		bgzip -k A_varcalls-man-merge_output.vcf && tabix A_varcalls-man-merge_output.vcf.gz
		cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/varcalls/
		bgzip -k B_varcalls-man-merge_output.vcf && tabix B_varcalls-man-merge_output.vcf.gz
		conda deactivate
		conda activate bcftools_env
		cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/varcalls/
		bcftools consensus -f ../../../C_albicans-[NAME]-A_Chrs-denovo_medaka.fasta -H 2 ./A_varcalls-man-merge_output.vcf.gz > C_albicans-[NAME]-A_Chrs.fasta
		cd /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/varcalls/
		bcftools consensus -f ../../../C_albicans-[NAME]-B_Chrs-denovo_medaka.fasta -H 2 ./B_varcalls-man-merge_output.vcf.gz > C_albicans-[NAME]-B_Chrs.fasta

#		Afterwards, do another alignment and check to make sure read-assembly agreement has improved!

# After one last verification step (i.e. scanning through with short and long reads aligned), move the FINALIZED ASSEMBLIES into a new folder!
# Follow the naming scheme below:

cd /[STRAINPATH]/
mkdir ./13_finalassembly/
cp /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/A_homs/varcalls/C_albicans-[NAME]-A_Chrs.fasta ./13_finalassembly/
cp /[STRAINPATH]/12_medaka_consensus-phased/homo_resolve/B_homs/varcalls/C_albicans-[NAME]-B_Chrs.fasta ./13_finalassembly/
cp /[STRAINPATH]/12_medaka_consensus-phased/C_albicans-[NAME]-ChrM-denovo_manual.fasta ./13_finalassembly/C_albicans-[NAME]-ChrM.fasta
}

-|-|-|-|-|-|-|-
'^IN PROGRESS^'
-|-|-|-|-|-|-|-
}


'MAKE THE ABOVE CLEANED ASSEMBLY A NEW STEP, i.e. 13, THEN HAVE THE BELOW PUT INTO IT'
{Annotating Assemblies via Lift Over and More!
	Liftoff v1.6.3 & SeqKit v2.5.0 # RAN ON OSC
		https://github.com/agshumate/Liftoff | https://anaconda.org/bioconda/liftoff & https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit

#		Prior to running, follow the installation guidelines as follows (some troubleshooting was required to get this to work).
		cd /[HOMEPATH]/
		git clone https://github.com/agshumate/Liftoff liftoff
		mamba create -n liftoff_env -c bioconda minimap2=2.17 python=3.12.0
		conda activate liftoff_env
		cd liftoff
		python setup.py install

#		Next, generate a file called "otherfeats.txt" with the following lines in it. This will provide Liftoff with a list of features to lift over, 
#		just to make sure it is pulling nearly everything over. Otherwise, it might miss some.
#			NOTE: If there are any additional features you know exist in your files and want to include, be sure to add them! Just follow the section
#			      called "Feature Types" from the GitHub page to make sure formatting is good.
		long_terminal_repeat
		retrotransposon
		tRNA
		ncRNA
		repeat_region
		blocked_reading_frame
		snoRNA
		snRNA
		centromere
		rRNA
		pseudogene

#		Next, make text files called "chroms_[A/B/ChrM].txt" and place text as follows in them to map Ca22 assemblies to yours 
#		(the first names come from the .gff files, specifically from how each chromosome is labelled in them):
		Ca22chr1A_C_albicans_SC5314,Chr1A
		Ca22chr2A_C_albicans_SC5314,Chr2A
		Ca22chr3A_C_albicans_SC5314,Chr3A
		Ca22chr4A_C_albicans_SC5314,Chr4A
		Ca22chr5A_C_albicans_SC5314,Chr5A
		Ca22chr6A_C_albicans_SC5314,Chr6A
		Ca22chr7A_C_albicans_SC5314,Chr7A
		Ca22chrRA_C_albicans_SC5314,ChrRA
#		Do the same for B and ChrM, too.

#		Then, place the "otherfeats.txt" and "chroms_[A/B/ChrM].txt" files into the following directory.
		cd /[STRAINPATH]/13_finalassembly/
		mkdir ./annotations/
		/[STRAINPATH]/13_finalassembly/annotations/otherfeats.txt
		/[STRAINPATH]/13_finalassembly/annotations/chroms_A.txt
		/[STRAINPATH]/13_finalassembly/annotations/chroms_B.txt
		/[STRAINPATH]/13_finalassembly/annotations/chroms_ChrM.txt
		
#		Then, install the LATEST annotations and sequence for your the reference assembly! This just makes sure whatever you're pulling in is
#		the most up-to-date.
		conda activate seqkit_env
		cd /[STRAINPATH]/13_finalassembly/annotations/
		wget 'http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/current/C_albicans_SC5314_A22_current_chromosomes.fasta.gz'
		gzip -d C_albicans_SC5314_A22_current_chromosomes.fasta.gz
		seqkit grep -nrp 'Ca22chr[1234567R]A_C_albicans_SC5314' ./C_albicans_SC5314_A22_current_chromosomes.fasta.gz | seqkit sort -w 70 -n > Ca_SC5314_A22-current_chrs-A.fasta
		seqkit grep -nrp 'Ca22chr[1234567R]B_C_albicans_SC5314' ./C_albicans_SC5314_A22_current_chromosomes.fasta.gz | seqkit sort -w 70 -n > Ca_SC5314_A22-current_chrs-B.fasta
		seqkit grep -nrp 'Ca22chrM_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -w 70 -n > Ca_SC5314_A22-current_chrs-ChrM.fasta
		wget 'http://www.candidagenome.org/download/gff/C_albicans_SC5314/archive/C_albicans_SC5314_A22_current_features.gff'
		grep -E '^[#*]|^Ca22.*A_C_alb*' C_albicans_SC5314_A22_current_features.gff > Ca_SC5314_A22-current_feats-A.gff
		grep -E '^[#*]|^Ca22.*B_C_alb*' C_albicans_SC5314_A22_current_features.gff > Ca_SC5314_A22-current_feats-B.gff
		grep -E '^[#*]|^Ca22.*M_C_alb*' C_albicans_SC5314_A22_current_features.gff > Ca_SC5314_A22-current_feats-ChrM.gff
		
#		This step is OPTIONAL, but makes things nicer - transfer over the gene names from A to B for ease of identification/avoiding having to
#		search a gene to find out what it's called! In RStudio (https://posit.co/download/rstudio-desktop/), perform the following to transfer
#		the [...]_feats-A.gff gene names to the [...]_feats-B.gff file.
			{	# Loads libraries for use.
				library(dplyr)
				library(data.table)
				library(emmeans)
				library(DescTools)
				library(car)
				
				# Pulls in the file for usage. Needs to be formatted the exact same every time. For the first file picker, use the 
				# "Ca_SC5314_A22-current_feats-A.gff" file, and for the second, "Ca_SC5314_A22-current_feats-B.gff".
				chrAdf = read.table(file.choose(), 
									header = FALSE,
									sep="\t")
				chrBdf = read.table(file.choose(), 
									header = FALSE,
									sep="\t")
				chrAdf
				chrBdf
				
				# The following will extract IDs from both files.
				get_id <- function(annotation) {
					id_match <- gsub(".*ID=(C[1234567R]_\\d+[A-Z]?).*", "\\1", annotation)
					return(id_match)
				}
				chrAdf$ID <- sapply(chrAdf$V9, get_id)
				chrBdf$ID <- sapply(chrBdf$V9, get_id)
				chrAdf
				chrBdf
				# The following will identify matching IDs between the two files and associate the lines with each other.
				matching_lines <- merge(chrAdf, chrBdf, by.x=c("V3", "ID"), by.y=c("V3", "ID"))
				matching_lines
				
				# The following transfer the "Gene=" section into the B versions.
				chrBupdf <- chrBdf
				for (i in 1:nrow(matching_lines)) {
					gene_annotation <- matching_lines[i, "V9.x"]
					existing_gene <- gsub(".*Gene=([^;]*);.*", "\\1", gene_annotation)
					start_section <- gsub("(ID=.*?);Note=.*", "\\1", matching_lines[i, "V9.y"], perl = TRUE)
					note_section <- gsub(".*Note=([^*]+).*", "\\1", matching_lines[i, "V9.y"])
					# Check if "Gene=" annotation exists in ChrA.
					if (any(grepl("Gene=", trimws(gene_annotation)))) {
						# Add "Gene=" annotation to ChrB between "Name=" and "Note=" sections.
						updated_annotation <- paste(start_section, ";Gene=", existing_gene, ";Note=", note_section, sep="")
						# Update ChrB with the updated annotation.
						chrBupdf[chrBupdf$V3 == matching_lines[i, "V3"] & chrBupdf$ID == matching_lines[i, "ID"], "V9"] <- updated_annotation
					}
				}
				chrBupdf$ID <- NULL
				chrBupdf
				
				# Saves down a new GFF3 file and removes the final "ID" column.
				write.table(chrBupdf, file="Ca_SC5314_A22-current_feats-B_gene.gff", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
			
#		Then, manually add the header back to the start of the file:
			##gff-version 3
			# File name: 
			# Organism: Candida albicans SC5314
			# Genome version: A22-s07-m01-r187
			# Date created: Sun Jul 16 07:00:57 2023
			# Created by: The Candida Genome Database (http://www.candidagenome.org/)
			# Contact Email: candida-curator AT lists DOT stanford DOT edu
			# Funding: NIDCR at US NIH, grant number 1-R01-DE015873-01
			#
			[rest of file]
				
#		Lastly, place the "Ca_SC5314_A22-current_feats-B_gene.gff" file into the 
		/[STRAINPATH]/13_finalassembly/annotations/
#		directory, where the other .gff files are located.

#		The following will run Liftoff for lifting over Assembly 22A annotations onto your "A" assembly, and the Assembly 22B annotations onto
#		your "B" assembly. Additionally, it will save unmapped features (ones which didn't weren't found) to their own file, plus it will perform
#		some polishing on start/stop positions in an attempt to restore in-frame coding sequencing if the liftover causes loss of a start/stop codon
#		or introduces an in-frame stop (it will make an unpolished file, too). Lastly, it will also first map features to their absolute best matches,
#		THEN search for any features that are 95+% identical and annotate them as copies (just in case some things were missed in the original
#		assembly). It is currently set to output ".gff" (GFF3) files. 
#			NOTE: I've checked the results, and it does a fairly good job for Candida albicans! It even correctly placed one copy of LEU2 into the 
#			      unmapped file, which was correct given that was the homolog with the deleted LEU2 locus! HOWEVER, it CAN sometimes get confused by
#			      expanded gene families such as the TLOs, where it DID map a few of them to the wrong locus (but it was always a TLO locus).
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate liftoff_env
		cd /[STRAINPATH]/13_finalassembly/annotations/
		liftoff -g ./Ca_SC5314_A22-current_feats-A.gff -o [NAME]_A_features.gff -u [NAME]_A_features-unmapped.txt -dir intermediate_A_files -f otherfeats.txt -polish -chroms chroms_A.txt -copies -sc 0.95 ../C_albicans-[NAME]-A_Chrs.fasta ./Ca_SC5314_A22-current_chrs-A.fasta
		liftoff -g ./Ca_SC5314_A22-current_feats-B_gene.gff -o [NAME]_B_features.gff -u [NAME]_B_features-unmapped.txt -dir intermediate_B_files -f otherfeats.txt -polish -chroms chroms_B.txt -copies -sc 0.95 ../C_albicans-[NAME]-B_Chrs.fasta ./Ca_SC5314_A22-current_chrs-B.fasta
		liftoff -g ./Ca_SC5314_A22-current_feats-ChrM.gff -o [NAME]_ChrM_features.gff -u [NAME]_ChrM_features-unmapped.txt -dir intermediate_ChrM_files -f otherfeats.txt -polish -chroms chroms_ChrM.txt -copies -sc 0.95 ../C_albicans-[NAME]-ChrM.fasta ./Ca_SC5314_A22-current_chrs-ChrM.fasta

#		Afterwards, I'd highly suggest going in and checking a few annotations that are likely to be problematic for assignments, such as TLOs and 
#		YRF1 genes (especially YRF1, as it isn't properly annotated in Assembly 22 and will require more in depth manual reannotation). If you find
#		that most everything looks good (i.e. annotations agree with start and stop codons of your new assembly), don't begin "fixing" things apart
#		from very clearly misannotated genes (again, like the TLOs). Our run got the TLOs mixed up a LOT. 
#			NOTE: I did find that using the "-chroms" option helped a LOT in resolving the mismapping of TLOs.
#			NOTE: I emphasize not "fixing" things because sometimes, a gene genuinely might not properly code for a protein/might have issues because
#			      the original assembly either got it wrong or missed some form of heterozygosity. Just be careful to not fix real events, basically.
#		For example...
		Copy "[NAME]_A_features.gff_polished" to a renamed "[NAME]_A_features_polished-man.gff", then manually change THAT file. This preserves the 
		original file, just in case you need to refer back or muck something up. Additionally, I suggest manually adding the chromosomes themselves 
		into the file as a new feature, just so that if youre using it in IGV it will show you where the ends of the chromosomes are. 
#			Example:
			Chr1A	Manual	chromosome	1	3240790	.	.	.	ID=Chr1A;Name=Chr1A

		***OPTIONAL PARTS PAST HERE***
#		If you REALLY, REALLY want to go the distance for checking the assembly/annotating, I highly suggest mapping both the homopolymer positions
#		AND finding ORFs de novo. Homopolymers of 8 or more CAN be iffy, so if a gene has a nonsense mutation due to a "frameshift" but there's a 
#		large homopolymer in it, it might actually be intact and the polishing got it wrong! GENERALLY trust it, but if you hit like 10-12 in the 
#		homopolymer, regard the "ORF/not ORF" ambiguity with a bit of caution. Finding ORFs de novo can be helpful to find what the ACTUAL coding
#		sequences are if your liftover annotations say something is wrong.
#			Homopolymers:
#				This will search for any stretches of 8+ of the same nucleotide. By separating each search then concatenating, it resolved a prior
#				niche issue due to the non-greedy system. Additionally, by renaming "location" to "homoloca", it allows the ORF prediction and 
#				homopolymer identification files to be loaded into different tracks. Also makes a .bed file for any variant-calling needs.
				source /[HOMEPATH]/miniconda3/bin/activate
				conda activate seqkit_env
				cd /[STRAINPATH]/13_finalassembly/annotations/
				seqkit locate --gtf -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-A_Chrs.fasta > A_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.gtf
				sed 's/location/homoloca/' A_8Homo-pre.gtf > A_8Homo.gtf
				rm ./A_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-B_Chrs.fasta > B_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.gtf
				sed 's/location/homoloca/' B_8Homo-pre.gtf > B_8Homo.gtf
				rm ./B_8Homo-pre.gtf
				cd /[STRAINPATH]/13_finalassembly/annotations/
				seqkit locate --bed -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-A_Chrs.fasta > A_8Homo-pre.bed
				seqkit locate --bed -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.bed
				seqkit locate --bed -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.bed
				seqkit locate --bed -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-A_Chrs.fasta >> A_8Homo-pre.bed
				sed 's/location/homoloca/' A_8Homo-pre.bed > A_8Homo.bed
				rm ./A_8Homo-pre.bed
				seqkit locate --bed -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-B_Chrs.fasta > B_8Homo-pre.bed
				seqkit locate --bed -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.bed
				seqkit locate --bed -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.bed
				seqkit locate --bed -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-B_Chrs.fasta >> B_8Homo-pre.bed
				sed 's/location/homoloca/' B_8Homo-pre.bed > B_8Homo.bed
				rm ./B_8Homo-pre.bed
				cd /[STRAINPATH]/13_finalassembly/annotations/
				seqkit locate --gtf -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-ChrM.fasta > ChrM_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.gtf
				seqkit locate --gtf -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.gtf
				sed 's/location/homoloca/' ChrM_8Homo-pre.gtf > ChrM_8Homo.gtf
				rm ./ChrM_8Homo-pre.gtf
				seqkit locate --bed -GrP -p "AAAAAAAAA*" ../C_albicans-[NAME]-ChrM.fasta > ChrM_8Homo-pre.bed
				seqkit locate --bed -GrP -p "TTTTTTTTT*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.bed
				seqkit locate --bed -GrP -p "GGGGGGGGG*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.bed
				seqkit locate --bed -GrP -p "CCCCCCCCC*" ../C_albicans-[NAME]-ChrM.fasta >> ChrM_8Homo-pre.bed
				sed 's/location/homoloca/' ChrM_8Homo-pre.bed > ChrM_8Homo.bed
				rm ./ChrM_8Homo-pre.bed

#			ORFs:
#				This will search for ANYTHING starting with an ATG and ending with any of the three canonical stop codons, as long as it is longer
#				than 222 nucleotides/74 amino acids (which could miss some smaller ORFs, but anything smaller is getting close to dubious). 
#				The initial "[...].txt" files contain positions for ALL possible ORFs, no filtering. The code for RStudio that follows is for
#				downstream filtering and producing a FILTERED .gtf file that follows the size limits I mentioned above (and fixes a niche issue
#				revolving around stop codons immediately after start codons being ignored).
#					There is a very, very, VERY important note to consider! This method will NOT predict most INTRON-CONTAINING ORFs, as it 
#					cannot account for the second exon laying anywhere out-of-frame - any intron-containing genes will need more manual verification 
#					by visuall checking the predictions or looking at RNA-Seq/related data. Thankfully, most lifted over annotations for 
#					intron-containing genes have seemed fairly accurate. There might be other tools out there to look for unnanotated intron genes!
				seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" -r ../C_albicans-[NAME]-A_Chrs.fasta > A_ORF_Pred-seqs.txt
				seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" -r ../C_albicans-[NAME]-B_Chrs.fasta > B_ORF_Pred-seqs.txt
				seqkit locate -i -p "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" -r ../C_albicans-[NAME]-ChrM.fasta > ChrM_ORF_Pred-seqs.txt

#				In RStudio (https://posit.co/download/rstudio-desktop/), perform the following to filter out ORFs smaller than 74 amino acids AND 
#				simultaneously remove overlapping ORFs that share the exact same stop position (i.e. assuming only the largest ORF with that stop
#				site is the real one). This might not always be true, but it's easy enough to find the next ATG if necessary.
				{	# Loads libraries for use.
					library(dplyr)
					library(data.table)
					library(emmeans)
					library(DescTools)
					library(car)
					
					# Pulls in "[A/B]_ORF_Pred-seqs.txt" using a file browser.
					orfdf = read.table(file.choose(), 
									   header = FALSE,
									   sep="\t",
									   stringsAsFactors = FALSE)
					orfdf
					
					# The following will remove any calls that have a stop codon IMMEDIATELY after the start codon (resolving a niche issue with the pattern matching method).
					trim_orfdf <- filter(orfdf, !grepl("^(ATGTAG.*|ATGTAA.*|ATGTGA.*)", orfdf$matched))
					trim_orfdf
					
					# The following will remove anything less than 222 nucleotides long (i.e. anything less than 74 amino acids gets removed).
					trim_orfdf$abs <- abs(trim_orfdf$start - trim_orfdf$end)
					subtrim_orfdf <- subset(trim_orfdf, abs > 221)
					subtrim_orfdf
					
					# The following will sort by seqID, strand, and start/end (based on direction), then only keep the largest ORF with the same stop (start/end) position. Also combines the two dataframes.
					subtrimpos_orfdf <- subset(subtrim_orfdf, strand == "+")
					subtrimpos_orfdf
					subtrimposcull_orfdf <- subtrimpos_orfdf %>% group_by(seqID, end) %>% top_n(1, abs(abs))
					subtrimposcull_orfdf
					
					subtrimneg_orfdf <- subset(subtrim_orfdf, strand == "-")
					subtrimneg_orfdf
					subtrimnegcull_orfdf <- subtrimneg_orfdf %>% group_by(seqID, start) %>% top_n(1, abs(abs))
					subtrimnegcull_orfdf
					
					finaltrim_orfdf <- rbind(subtrimposcull_orfdf, subtrimnegcull_orfdf)
					finaltrim_orfdf
					
					# Creates a new dataframe that mimics a .gtf file.
					chr <- c(finaltrim_orfdf$seqID)
					start <- c(finaltrim_orfdf$start)
					end <- c(finaltrim_orfdf$end)
					strand <- c(finaltrim_orfdf$strand)
					
					printme_orfdf <- data.frame(chr)
					printme_orfdf['meth'] <- "SeqKit"
					printme_orfdf['loca'] <- "location"
					printme_orfdf$start <- start
					printme_orfdf$end <- end
					printme_orfdf['zero'] <- "0"
					printme_orfdf$strand <- strand
					printme_orfdf['dot'] <- "."
					printme_orfdf['pattern'] <- "gene_id A[TU]G(?:.{3})+?[TU](?:AG|AA|GA);"
					printme_orfdf
					
					# Saves down .csv file - rename as necessary (make it .gtf).
					write.table(printme_orfdf, "[A/B/ChrM]_ORF_PredFilt.gtf", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)}
}