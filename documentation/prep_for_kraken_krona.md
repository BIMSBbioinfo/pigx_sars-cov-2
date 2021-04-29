## Prepare Database for Kraken2 

*After the Kraken2 is installed with guix the following steps need to be done in order to download the database (all genomes the reads should be classified against) and install it correctly*
*For the user I guess it is important to download the database and then provide the location to the settings-file (?), so that the rest can run via the pipeline <https://github.com/DerrickWood/kraken2/wiki/Manual>*

| Steps | Commands | Notes |
|-----|---------|----------|
| set $DBNAME | | the folder, where the database should be installed in |
| download db |. | So far we used the PlusPFP from <https://benlangmead.github.io/aws-indexes/k2> |
| extract db | tar -xvf [] | extract the db |
| download taxonomy | kraken2-build --download-taxonomy --db $DBNAME | builds the taxonomy in that folder|
| build db | kraken2-build --build --db $DBNAME | |
| sample sheet | | provide loacation to sample sheet ( KRAKEN_DIR) |


Now the kraken command can be run after the unaligned_reads.bam file is converted to a fastq file. This is done by snakefile.py:  
*samtools fastq [input sam] > [output fastq]*  
*kraken2 --report [output] --db $DBNAME --paired [input fastq R1 R2]*  

## Prepare Database for Krona 

*Before KronaTools can run there are also preparation steps. First load with Guix and then do the following*

| Steps | Commands | Notes |
|-----|---------|-----------|
| create $folder | | where to put the taxonomy files |
|update Taxonomy | ~/.guix-profile/share/krona-tools/updateTaxonomy.sh $folder | be sure to specify the folder |
|update accession | ~/.guix-profile/share/krona-tools/updateAccessions.sh $folder | |
|sample sheet | | provide loacation to sample sheet|

Now the Kronatool can run in the snakefile.py:  
*ktImportTaxonomy -m 3 -t 5 [kraken input] -tax $KRONA_DIR -o [output.html]*


## NOTES for renderSite

* (1) snakefile.py a script called generateSiteFiles.R is run, this will create the necessary files for the website  
* (2) snakefile runs the rendering of the .RMD scripts provided in a subfolder of scripts/report_scripts/
* (3) snakefile.py runs the renderSite() command to generate the website based on (1) and (2)

