## Installing Seqr on a fresh VM (Ubuntu 16.04)
Brian Uapinyoying
07/26/2018

These notes are based off of the [official instructions for installing on a MacOS machine](https://github.com/macarthur-lab/seqr/tree/master/deploy/mac_osx#loading-your-own-datasets)

1. Update and upgrade apt-get: 
```bash
sudo apt-get update && sudo apt-get upgrade
```

2. Install python 2.7: `sudo apt-get python`

3. Install python pip and upgrade it
```bash
sudo -H apt-get install python-pip
sudo -H pip2 install --upgrade pip
```
4. Make an install directory for seqr.  I am going to try it in the home directory of my VM.
	- `/home/<user>/seqr`

5. Create an bash variable or alias:
```bash
export SEQR_INSTALL_DIR="/home/<user>/seqr"
```
	- or more perminantly add it as a line to the `.bash_profile` and run `source .bash_profile`

6. Create a few directories inside the seqr folder
```bash
mkdir code data data/reference_data data/projects
```

7. Download the seqr reference data, and uncompress
```bash
cd ${SEQR_INSTALL_DIR}/data/reference_data
wget https://storage.googleapis.com/seqr-reference-data/seqr-resource-bundle.tar.gz
tar -xzf seqr-resource-bundle.tar.gz
```

8. Install MongoDB
```bash
sudo apt-get install mongodb
mongo # to test if it works, Ctr+D to quit
```

9. Install Postgres
[Instructions on how to install postgres on Ubuntu 16.04](https://www.linode.com/docs/databases/postgresql/how-to-install-postgresql-on-ubuntu-16-04/)
```bash
sudo apt-get install postgresql postgresql-contrib
# Don't change the password to the postgres user, causes problems
```

10. Edit`/etc/postgresql/9.5/main/pg_hba.conf` to change permission settings to make postgres work
```bash
# Database administrative login by Unix domain socket
local   all             postgres                                peer

# TYPE  DATABASE        USER            ADDRESS                 METHOD

# "local" is for Unix domain socket connections only
local   all             all                                     trust
# IPv4 local connections:
host    all             all             127.0.0.1/32            trust
# IPv6 local connections:
host    all             all             ::1/128                 trust
# Allow replication connections from localhost, by a user with the
# replication privilege.
#local   replication     postgres                                peer
#host    replication     postgres        127.0.0.1/32            md5
#host    replication     postgres        ::1/128                 md5
```

11. Restart the service, create a seqr database and test out a connection to the postgresql database
```bash
sudo service postgresql restart

createdb seqrdb

psql -U postgres # Check connection with db, quit with ctr+D

```

12. Installing the official Oracle Java Development kit 8 (JDK v1.8) for PhenoTips to work properly.  Don't get the java runtime enviornment (JRE), its not enough. Although you can also use the default-jdk.
```bash
sudo apt-get update && apt-get upgrade
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java8-installer

# OR
sudo apt-get install default-jdk

# This will **NOT** work:
#sudo apt-get install default-jre
```

11. Install PhenoTips for storing structured phenotype information. I installed it to the seqr home directory
```bash
wget https://nexus.phenotips.org/nexus/content/repositories/releases/org/phenotips/phenotips-standalone/1.2.6/phenotips-standalone-1.2.6.zip
# if its not installed
# sudo apt-get install unzip
unzip phenotips-standalone-1.2.6.zip
rm phenotips-standalone-1.2.6.zip
cd phenotips-standalone-1.2.6
./start.sh # run it
```

12. Clone the seqr repo from github into the code directory and add the PYTHONPATH variable to .bash_profile
```bash
cd ${SEQR_INSTALL_DIR}/code/seqr
git clone https://github.com/macarthur-lab/seqr.git

# add these lines to .bash_profile
export PYTHONPATH=${SEQR_INSTALL_DIR}/code/seqr:$PYTHONPATH
export PYTHONPATH=${SEQR_INSTALL_DIR}/code/seqr/deploy/mac_osx/xbrowse_settings:$PYTHONPATH

# source it
source .bash_profile
```

12. Install seqr's python dependencies
```bash
cd ${SEQR_INSTALL_DIR}/code/seqr
sudo -H pip install -r requirements.txt
```

### Part III. Initialize Seqr

1. Run these commands to initialize postgres and load reference data:
```bash
cd ${SEQR_INSTALL_DIR}/code/seqr
# django command for initializing seqr database tables
./manage.py migrate # Broke here

# create a new django admin user
./manage.py createsuperuser
```

2. In order to do the next step, `./manage.py load_resources` properly, you will need to do several steps to fix the seqr code and provide missing files so Clinvar can load properly. Start off by copying over everything in the reference_data folder from spinup1 to the new VM: `/seqr/data/reference_data/`

	+ This `refeence_data` folder should include:

		- ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.decomposed.with_popmax.vcf.gz
		- cleaned_exac_with_pHI_march16_pLI.csv.  **# in resource bundle**
		- clinvar_alleles.single.b37.tsv
		- clinvar_alleles.single.b37.vcf.gz
		- clinvar_alleles.single.b37.vcf.gz.tbi
		- **clinvar.tsv** -> clinvar_alleles.single.b37.tsv
		- **clinvar.vcf.gz** -> clinvar_alleles.single.b37.vcf.gz
		- **clinvar.vcf.gz.tbi** -> clinvar_alleles.single.b37.vcf.gz.tbi
			+ Above bolded files are softlinks, make sure they are still linked properly
		- ExAC.r0.3.final.vep.popmax.clinvar.vcf.gz
		- ExAC.r0.3.final.vep.popmax.clinvar.vcf.gz.tbi
		- fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
		- forweb_cleaned_exac_r03_2015_03_16_z_data_missense.csv. # in resource bundle
		- gencode.v19.annotation.gtf.gz
		- gencode.v27lift37.annotation.gtf.gz  # in resource bundle
		- gene_constraint_scores.csv  **# in resource bundle**
		- genemap2.txt
		- GTEx_Data_V6_Annotations_SampleAttributesDS.txt  **# in resource bundle**
		- GTEx_samples.txt
		- high_variability.genes.txt  **# in resource bundle**
		- omim/genemap2.txt  # in resource bundle
		- GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz  **# in resource bundle**
			+ Can be downloaded from https://www.gtexportal.org/home/datasets and a login is required. 
			+ Navigate to `gtex_analysis_v6` → `RNA-Seq Data section` → `Gene RPKM (1.9 gb)`
			+ Leave gzipped, the name is different if you uncompress

	+ Some files can be found online:
```bash
wget https://github.com/macarthur-lab/clinvar/blob/master/output/b37/single/clinvar_alleles.single.b37.vcf.gz
wget https://github.com/macarthur-lab/clinvar/blob/master/output/b37/single/clinvar_alleles.single.b37.vcf.gz.tbi
```

3. Add a new line to: `/seqr/code/seqr/deploy/utils/constants.py`
```python
REFERENCE_DATA_FILES = {  # Line 79 of the script
    'gencode': 'gencode.v27lift37.annotation.gtf.gz',
    'high_variability_genes': 'high_variability.genes.txt',
    'constraint_scores': 'constraint_scores.csv',
    'gtex_expression': 'GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz',
    'gtex_samples': 'GTEx_Data_V6_Annotations_SampleAttributesDS.txt',
    'omim_genmap': 'omim/genemap2.txt',
    'clinvar': 'clinvar.tsv',
    'clinvar_vcf': 'clinvar.vcf.gz',   # <-- ADD line
    'dbnsfp': 'dbNSFP3.5_gene'
}
```

4. Add new lines to: `/seqr/code/seqr/deploy/mac_osx/xbrowse_settings/reference_settings.py`
```python
gene_test_statistic_tags = [    # line 29 of script
    {
        'slug': 'lof_constraint',
        'data_field': 'pLI',
        'file_path': xbrowse_reference_data_dir + '/cleaned_exac_with_pHI_march16_pLI.csv' # <-- add line
    },
    {
        'slug': 'missense_constraint',
        'data_field': 'mis_z',  # DELETE this line
        'data_field': 'zscore', # ADD this line
        'file_path': xbrowse_reference_data_dir + '/forweb_cleaned_exac_r03_2015_03_16_z_data_missense.csv' # <--add line
    }
]

# add the below line near the bottom of the script, ~ line 45
clinvar_vcf_file = os.path.join(xbrowse_reference_data_dir, REFERENCE_DATA_FILES['clinvar_vcf'])
```

5. Add and remove lines from : `/seqr/code/seqr/xbrowse/reference/reference.py`
```python
def _load_gene_test_statistic_tags(self):   # Function starts on line 243 in the script
    
    score_data = pandas.DataFrame.from_csv(self.settings_module.constraint_scores_file) # <-- DELETE this line

    for gene_tag in self.settings_module.gene_test_statistic_tags:
        score_data = pandas.DataFrame.from_csv(gene_tag['file_path']) # <-- add this line
        tag_id = gene_tag['slug']
        tag_field = gene_tag['data_field']

        # first set all to None
        self._db.genes.update({}, {'$set': {'tags.'+tag_id: None}}, multi=True)

        scores = getattr(score_data, tag_field)
        ranks = scores.rank(ascending=False)
        for gene_id, score in scores.iteritems():
            self._db.genes.update({'gene_id': gene_id}, {'$set': {'tags.'+tag_id: score}})
        for gene_id, rank in ranks.iteritems():
            self._db.genes.update({'gene_id': gene_id}, {'$set': {'tags.'+tag_id+'_rank': [int(rank), len(ranks)]}})
```

6. Edit: `/seqr/code/seqr/deploy/mac_osx/xbrowse_settings/annotator_settings.py`

	+ Remove 1000 Genome reference
		- Comment out everything from line 13-26 in the script

	+ Change the word 'sites' in these two lines
		- 'file_path': '%(xbrowse_install_dir)s/data/reference_data/ExAC.r0.3.**sites**.vep.popmax.clinvar.vcf.gz' % locals(),
		- 'file_path': '%(xbrowse_install_dir)s/data/reference_data/ExAC.r0.3.**sites**.vep.popmax.clinvar.vcf.gz' % locals(),

	+ Into 'final' like so:
		- 'file_path': '%(xbrowse_install_dir)s/data/reference_data/ExAC.r0.3.**final**.vep.popmax.clinvar.vcf.gz' % locals(),
		- 'file_path': '%(xbrowse_install_dir)s/data/reference_data/ExAC.r0.3.**final**.vep.popmax.clinvar.vcf.gz' % locals(),


7. Now load the resources
```bash
# load reference data - this may take over 1 hour
./manage.py load_resources

# **Optional**: If the script gets interrupted at some point and you need to fix a bug, 
# you can comment out the functions for resources that have successfully been loaded to save time.
# File: ./seqr/code/seqr/xbrowse/reference/reference.py

def load(self):  # function starts line #75
#    self._load_clinvar()
#    self._load_genes()
#    self._load_additional_gene_info()
#    self._load_tags()
    self._load_gtex_data()
    self._reset_reference_cache()
```

8. Create a test project and give seqr a test drive
```bash
# create an empty new project
./manage.py add_project test_project

# start the development server
./manage.py runserver <your.ip>:8000
# example the command for my VM was:
# ./manage.py runserver 10.5.34.3:8000

# Open the page in the browser window. You can log in with the username and password 
# provided to the createsuperuser command in step 1.
open http://<your.ip>:8000
# e.g. http://10.5.34.3:8000
```

### Part IV. Loading your own data set.

NOTE: seqr expects VCFs to have the following genotype format: GT:AD:DP:GQ:PL

Before a VCF can be loaded into mongodb, it must be annotated with Variant Effect Predictor (VEP) to add a specific set of annotations, including those provided by the dbNSFP and LoFTEE plugins. Run the command below to install VEP. Also, it's useful to install tabix in order to optimize the VEP cache.

1. Install tabix and VEP, this is probably something you want to do on a cluster. It needs a lot of memory:
```bash
sudo apt-get install tabix

# Consider installing VEP on a cluster to process variants
cd ${SEQR_INSTALL_DIR}
wget https://github.com/Ensembl/ensembl-tools/archive/release/85.zip
unzip 85.zip 
mv ensembl-tools-release-85/scripts/variant_effect_predictor .
rm -rf 85.zip ensembl-tools-release-85
cd variant_effect_predictor
perl INSTALL.pl --AUTO acf --CACHEDIR ../vep_cache --SPECIES homo_sapiens --ASSEMBLY  GRCh37 --CONVERT

# Once vep and the plugins have been installed, you can use the following command to annotate your VCF file (here called my_data.vcf.gz):
perl ./vep/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl --everything --vcf --allele_number --no_stats --cache --offline --dir ./vep_cache/ --force_overwrite --cache_version 81 --fasta ./vep_cache/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa --assembly GRCh37 --tabix --plugin LoF,human_ancestor_fa:./loftee_data/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15 --plugin dbNSFP,./reference_data/dbNSFP/dbNSFPv2.9.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred -i my_data.vcf.gz -o my_data.vep.vcf.gz
```
	- NOTE: VEP tends to run out of memory on large VCFs, so it's best to split the vcf into chuncks with 5000 or fewer variants in each, run VEP on each chunk in parallel, and then recombine. The grabix indexing tool is very helpful for the splitting step as it lets you extract an arbitrary range of lines from the vcf, and these can be piped into VEP.

2. Once you have an annotated VCF and a .fam/.ped file describing the pedigree of your samples, load it
```bash
# Add families and individuals to the project, could also be a '.ped' file instead of '.fam'
./manage.py add_individuals_to_project test_project --ped 'test_project_samples.fam'   

# I am going to test this out using some data Monkol got from 1k exomes that he 
# already processed through VEP, I placed it in a new folder i created 
# mkdir seqr/data/projects/1kg_project <-- 1kg.ped + 1kg.vep_minimal.vcf.gz
./manage.py add_project 1kg_project
./manage.py add_individuals_to_project 1kg_project --ped '/home/<user>/seqr/data/projects/1kg_project/1kg.ped'
```

3. Now set the VCF locations, doesn't load the data yet. Make sure to use VEP generated VCF
```bash
./manage.py add_vcf_to_project test_project my_data.vep.vcf.gz

# Again, my personal command: 
./manage.py add_vcf_to_project test_project /home/<user>/seqr/data/projects/1kg_project/1kg.vep_minimal.vcf.gz
```

4. This step actually loads the VCF data and takes a while
```bash
./manage.py load_project test_project

# Mine:
./manage.py load_project 1kg_project
```

5. Load an additional mongodb collection that is used for gene search - this also takes a while
```bash
./manage.py load_project_datastore test_project

# Mine:
./manage.py load_project_datastore 1kg_project  
```

For other seqr options:
```bash
./manage.py help
```
