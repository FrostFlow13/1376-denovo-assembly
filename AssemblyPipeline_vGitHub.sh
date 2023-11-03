# I highly suggest viewing this document in a program like Notepad++ - it will make navigating/looking at the document a lot nicer.

# Additionally, as a general rule of thumb, if a section of code DOESN'T start with "source /[HOMEPATH]/miniconda3/bin/activate", that
# means it's very likely something that needs to be done via the terminal! This is because the jobs need conda itself activated first.

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
#	as provided by each group. Pilon and DeconSeq and HapDup had special installation instructions - see their sections for setup.
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
#				Except for Flye (~3 hours), as it sometimes needs more time.
#			For the line with "nodes" and whatnot, DEFINITELY double check the right arguments/inputs before running. Just to be sure.
#				Plus, if you can use fewer cores/don't need as many, you'll get stuff running sooner!

#	After this, just run the script in the supercomputer via terminal "sbatch [NAME].sh" and it'll set the job up automatically.
#	If necessary, you can use "scancel [JOB#]" and that will halt the job running.
}

-------------------------------------------------------------------

{Moving Files
# For moving stuff back and forth between OSC and GORT:
#	GORT -> OSC
	scp src /srv/work/[NAME]/[PATH&FILE] [name#]@sftp.osc.edu:/fs/ess/PAS[####]/[NAME]/[PATH]/[DESTINATION]/
#		e.g., scp src /srv/work/Andrew_Woodruff/2023_06_15-MAY1376_TLOKOs_LongRead/1741/demul_adtrim/BC18.fastq woodruff207@sftp.osc.edu:/fs/ess/PAS[####]/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1741/demul_adtrim/
#		Note: doesn't necessarily have to be server from GORT. Can also be frome /home/[PATH]/ stuff.

#	OSC -> GORT
	scp src [name#]@sftp.osc.edu:/fs/ess/PAS[####]/[NAME]/[PATH&FILE] /srv/work/[NAME]/[PATH]/[DESTINATION]/
#		e.g., scp src woodruff207@sftp.osc.edu:/fs/ess/PAS[####]/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1741/demul_adtrim/BC18.fastq /srv/work/Andrew_Woodruff/2023_06_15-MAY1376_TLOKOs_LongRead/1741/demul_adtrim/	
}

-------------------------------------------------------------------

# Some general terms explained:
#	For paths, you must update these manually to your own locations.
	"/[HOMEPATH]/" # default user path with miniconda3  (OSC e.g., "/users/PAS[####]/woodruff207/" ; GORT e.g., "/home/woodruff.207/")
	"/[STRAINPATH]/" # server path with your data in it (OSC e.g., "/fs/ess/PAS[####]/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1376/" ; 
#                                                             GORT e.g., "/srv/work/Andrew_Woodruff/2023_06_15-MAY1376_TLOKOs_LongRead/1376/")
	"[NAME]" # the name of your strain or some identifier you want the files started with
	Anything Else = # whatever else needs to be input = look at the script to tell what you should put in!

	[TOOL] # RAN ON GORT = it was ran on GORT as a simple script.
	[TOOL] # RAN ON OSC = it was ran on OSC as a JOB (more complex .sh file, unless stated otherwise).


{Sequencing & Basecalling & Demultiplexing
#	Already performed for us by UW-M.

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
#		Flye assembly (expects 1 kb minimum, prefers 5+ kb).
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
#			***NOTE2: 25 kb makes EXTREMELY consistent graphs with Flye! However, check for telomeric sequences - I've noticed
#			   a few contigs that hit telomeric sequence, then fuse something else on AFTER. Trim those (and the improper extra
#			   sequence) out of the resulting .fasta files in the chromosome-level assembly step coming up!***
#				   TRUE 23-bp telomerase repeat (correct order for looping, and with 5'-3' order relative to reference genome)
					5'-ACACCAAGAAGTTAGACATCCGT-3' (left arm)
					5'-ACGGATGTCTAACTTCTTGGTGT-3' (right arm)
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/1_demul_adtrim/
		flye --nano-hq [PORECHOPPED]-25kbmin.fastq --out-dir ../2_flye_assembly --threads 28
#			The "PORECHOPPED" refers to whatever name Porechop_ABI gives the file post-demultiplexing/adapter trimming (e.g. BC##).

#		After, I suggest identifying which contigs in assembly.fasta correspond to which chromosomes, then rename the fasta
#		headers for ease of identification in the following steps. Use Bandage for visualizing the .gfa assembly to get
#		the best visual reference of each assembly run.

#		Rename the original file to something else like
		assembly-ORIGINAL.fasta
#		, then upload the file with chromosome names re-added as
		assembly.fasta

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
		busco -i ../2_flye_assembly/assembly.fasta -l saccharomycetes_odb10 -o ./2_busco_flye_results -m genome -c 28
		conda deactivate
		conda activate seqkit_env
		cd /[STRAINPATH]/busco_quast/
		wget 'http://www.candidagenome.org/download/sequence/C_albicans_SC5314/Assembly22/archive/C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz'
		gzip -d C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta.gz
		seqkit grep -nrp 'Ca22chr[1234567R]A_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -n > Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta
		seqkit grep -nrp 'Ca22chr[1234567R]B_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -n > Ca_SC5314_A22-s07-m01-r187_chrs-B.fasta		
		seqkit grep -nrp 'Ca22chrM_C_albicans_SC5314' ./C_albicans_SC5314_version_A22-s07-m01-r187_chromosomes.fasta | seqkit sort -n  > Ca_SC5314_A22-s07-m01-r187_chrs-ChrM.fasta
		wget 'http://www.candidagenome.org/download/gff/C_albicans_SC5314/archive/C_albicans_SC5314_version_A22-s07-m01-r187_features.gff'
		grep -E '^[#*]|^Ca22.*A_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-A.gff
		grep -E '^[#*]|^Ca22.*B_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-B.gff
		grep -E '^[#*]|^Ca22.*M_C_alb*' C_albicans_SC5314_version_A22-s07-m01-r187_features.gff > Ca_SC5314_A22-s07-m01-r187_feats-ChrM.gff
		conda deactivate
		conda activate quast_env
		python /[HOMEPATH]/miniconda3/envs/quast_env/bin/quast.py ../2_flye_assembly/assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-A.fasta -g Ca_SC5314_A22-s07-m01-r187_feats-A.gff -o ./2_quast_flye_results -t 28
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
		perl /[HOMEPATH]/deconseq/deconseq.pl -f ../2_flye_assembly/assembly.fasta -dbs camito -out_dir ./ -id [NAME]
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
#		then run minimap2, mapping the [PORECHOPPED]-[10/25]kbmin long reads against the polished assembly.
#			***NOTE: for the 'sed' step, replace the number of pilons with the number you have in your last polished file.***
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate seqkit_env
		cd /[STRAINPATH]/6_pilon/pilon_r[#]/
		sed 's/pilon_pilon_pilon_pilon_pilon_pilon_pilon/pilon-r7/' [NAME]-r[#].fasta > [NAME]-r[#]-rename.fasta
		seqkit sort -n [NAME]-r[#]-rename.fasta > [NAME]-collapsed-r[#].fasta
		conda deactivate
		conda activate flye_env
		cd /[STRAINPATH]/6_pilon/pilon_r[#]/
		flye-minimap2 -ax map-ont -t 28 [NAME]-collapsed-r[#].fasta ../../1_demul_adtrim/[PORECHOPPED]-10kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-collapsed-r[#]-10kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-collapsed-r[#]-10kbmin-mapping.bam
		flye-minimap2 -ax map-ont -t 28 [NAME]-collapsed-r[#].fasta ../../1_demul_adtrim/[PORECHOPPED]-25kbmin.fastq | flye-samtools sort -@ 28 -m 4G > [NAME]-collapsed-r[#]-25kbmin-mapping.bam
		flye-samtools index -@ 28 [NAME]-collapsed-r[#]-25kbmin-mapping.bam

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

#		Then, go to any large homozygous region (or region you particularly care about) in IGV, and try to identify indels from the
#		.vcf files that look real (i.e. 20+ quality score). If you can find one, try to figure out which allele it belongs to, then
#		manually annotate that in a copy of the .vcf file. This will help things down the line properly pull down reads! Additionally,
#		try to use some LOH'd strains to figure out alleles that fall FAR outside the region where phasing can be successfully done (i.e.
#		for SC5314, the far right arm of Chr3 or the far arms of Chr7 and ChrR).
#			The way I did this was by loading up the whatshap .vcf, pasting the longphase .vcf "SAMPLE" column next to the whatshap one,
#			then going through the homozygous or homozygous-adjacent regions and making sure the two agree, and changing anything that
#			looked wrong or that looked real in the file/IGV. Afterwards, I copied the information into a new .vcf file (a copy of the
#			original), excluding the extra "SAMPLE" column from longphase. By changing, I mean turning "/" to "|" if I thought it was real
#			(and assigning it to the same phase block as the rest of the chromosome if it was different), or turning it from "|" to "/" if 
#			it seemed like it wasn't real in IGV. If you can figure out WHICH homolog the variants out past homozygous regions belong to
#			(i.e. near the ends), make sure that the "#|#" is assigned correctly, too! The order tells programs which haplotype the 
#			REF|ALT belong to, and if it had put them in their own phaseblock, it might have the wrong order for the OTHER phaseblock.
#				*Don't forget to either add or remove or add the ":PS" in the "FORMAT" column, too! 
#				**For my stuff, it looked like REAL INDELS usually had a quality greater than 20, and real SNPs had a quality greater
#				  than 5 or 6. HOWEVER, some INDELS that are real do fall below the 20 quality - I've just found that if the indel is
#				  above 20, it is almost guaranteed to be real and not to be a homopolyer/repetitive issue.
#				***Lastly, if you find any regions where there is clear heterozygosity (and probably some sort of inversion) and you
#				   still see a large chunk of reads not assigned to the right haplotype, it might be calling the reads incorrectly 
#				   because the reads inside the inversion begin getting mixed with the reads outside the inversion (you'll usually
#				   see a "wall" if this happens). You might be able to turn on/turn off called SNPs in that region to get better 
#				   assignments. For example, TCA4-4 on Chr4.

#		Upload the new file into the "/7_clair3-phasing/" directory with the name, 
		dual-man-phased_merge_output.vcf
}



At this point, you should have a PRETTY good collapsed assembly that can be used for most things, including knowing where MOST of the
smaller variants are in the strain relative to your collapsed reference.

Beyond this point is the process for generating a phased assembly.



{Haplotagging Reads & Read Splitting - Assigning Haplotypes to Each Long Read & Splitting Reads Based on Haplotype Assignment
	Clair3 v1.0.4 # RAN ON OSC
		https://github.com/HKU-BAL/Clair3 | https://anaconda.org/bioconda/clair3

#		The following will run whatshap on the collapsed assembly to tag each long read inside "[NAME]-collapsed-r[#]-10kbmin-mapping.bam" 
#		with a haplotype. This will then be used to split the "hap1" and "hap2" reads into two entirely separate files, thus allowing FULLY
#		reassembly of the homologous chromosomes in the next steps! It will also throw all unassigned reads into BOTH phased hap files,
#		ensuring that homozygous regions still get properly represented during assembly. There MAY be some reads that are from the other
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
		whatshap split --output-h1 [NAME]-hap1_reads.bam --output-h2 [NAME]-hap2_reads.bam --add-untagged [NAME]-haplotagged.bam [NAME]-haplotagged_list.tsv

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
}


{Genome Reassembly - Phased
	Flye v2.9.2 & SeqKit v2.5.0 # RAN ON OSC
		https://github.com/fenderglass/Flye | https://anaconda.org/bioconda/flye & https://github.com/shenwei356/seqkit | https://anaconda.org/bioconda/seqkit

#		The following will run Flye using assumptions for ONT data basecalled with Guppy5+ SUP, grab the [NAME]-hap[1/2]_reads.bam
#		files, convert them to .fastq, then run Flye. After running, it places the resulting directory  up into the folder of the 
#		strain being processed for  ease of access in the next step. Probably allocate at least 4-5 hours in a job for this! Additionally, 
#		it ends up using minimum-length 25kb reads to minimize risk of messing things up, but it still will keep the 10kbmin split .fastq 
#		files for use in future steps/alignments.
#			***NOTE2: Check for telomeric sequences - I've noticed a few contigs that hit telomeric sequence, then fuse something else 
#			   on AFTER. Trim those (and the improper extra sequence) out of the .fasta in the chromosome-level assembly step coming up!***
#				   TRUE 23-bp telomerase repeat (correct order for looping, and with 5'-3' order relative to reference genome)
					5'-ACACCAAGAAGTTAGACATCCGT-3' (left arm)
					5'-ACGGATGTCTAACTTCTTGGTGT-3' (right arm)
		source /[HOMEPATH]/miniconda3/bin/activate
		conda activate flye_env
		cd /[STRAINPATH]/
		mkdir ./9_flye_assembly-phased/
		cd /[STRAINPATH]/8_clair3-read_split/
		flye-samtools fastq [NAME]-hap1_reads.bam > [NAME]-hap1_reads.fastq
		flye-samtools fastq [NAME]-hap2_reads.bam > [NAME]-hap2_reads.fastq
		conda deactivate
		conda activate seqkit_env
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
#			NOTE2: The telomere-proximal portion of Chr3R and Chr7L are currently randomly assigned, as the variants
#			(SNPs for Chr3R, LEU2+/- for Chr7L) are too far away from any other variants to assign a haplotype, and
#			we have no strains for Chr7L that could be used to determine which homolog has which variant. BWP17 (MAY307)
#			could be used for Chr3R, but I have not had time to do this yet.

#		Upload the file with chromosome names and homologs re-added as
		[NAME]-hap[1/2]-iden_assembly.fasta

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
		source /users/PAS1802/woodruff207/miniconda3/bin/activate
		conda activate busco_env
		cd /fs/ess/PAS1802/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1376/busco_quast/
		busco -i ../9_flye_assembly-phased/flye_hap1/1376_hap1-iden_assembly.fasta -l saccharomycetes_odb10 -o ./9_busco_flye-hap1_results -m genome -c 28
		busco -i ../9_flye_assembly-phased/flye_hap2/1376_hap2-iden_assembly.fasta -l saccharomycetes_odb10 -o ./9_busco_flye-hap2_results -m genome -c 28
		conda deactivate
		conda activate quast_env
		cd /fs/ess/PAS1802/ALW/2023_06_15-MAY1376_TLOKOs_LongRead/1376/busco_quast/
		python /users/PAS1802/woodruff207/miniconda3/envs/quast_env/bin/quast.py ../9_flye_assembly-phased/flye_hap1/1376_hap1-iden_assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./9_quast_flye-hap1_results -t 28
		python /users/PAS1802/woodruff207/miniconda3/envs/quast_env/bin/quast.py ../9_flye_assembly-phased/flye_hap2/1376_hap2-iden_assembly.fasta -r Ca_SC5314_A22-s07-m01-r187_chrs-[A/B].fasta -g Ca_SC5314_A22-s07-m01-r187_feats-[A/B].gff -o ./9_quast_flye-hap2_results -t 28}
}


Beyond this, you can follow the "Chromosome-Level Assembling & Checking" section near the beginning to reassemble your genomes, and you can 
also use your phased .fastq files for long-read polishing using Medaka. However, I have not yet solved phasing short-reads for final polishing
of phased assemblies, as I have concerns regarding expanded gene families and the potential for accidentally introducing variants that are
actually from extremely similar paralogs.