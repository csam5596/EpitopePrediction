# Neo Epitope Prediction Workflow

The here presented workflow derives neo eptiopes from raw paired end sequence data from 
Tumor and Normal sample DNA as well as Tumor RNA.

This workflow is intended to be used with the [Cromwell Execution Engine](https://cromwell.readthedocs.io/en/stable/)

The provided configuration file is set up to use SGE as task scheduler and Singularity as container technology.

## Execution

In order for the pipeline to be executed the **Epi.simg** singularity file from 
https://github.com/csam5596/Epi needs to be compiled into
a container image with the name **epi.sif** and placed in the same folder from which the pipeline
is executed.

The pipeline can then be executed by calling:
```console
java -Dconfig.file=configuration.conf -jar cromwell.jar run epi.wdl --inputs input.json
```

## Input

The following parameters need to be specified in an json file which serves as input to the pipeline:

| Name               | Description |
|--------------------|-------------|
| seqNormal1         | Forward DNA reads from normal sample (**.fastq.gz**) |
| seqNormal2         | Reverse DNA reads from normal sample (**.fastq.gz**) |
| seqTumor1          | Forward DNA reads from tumor sample (**.fastq.gz**) |
| seqTumor2          | Reverse DNA reads from tumor sample (**.fastq.gz**) |
| seqRNA1            | Forward RNA reads from tumor sample (**.fastq.gz**) |
| seqRNA2            | Reverse RNA reads from tumor sample (**.fastq.gz**) |
| seqAdapters        | File containing adapter sequences introduced in the sequencing process (**.fa**) |
| cropLength         | Maximum length to trim all raw sequences to |
| refGenomeFolder    | Path to folder containing the reference genome |
| refGenomeName      | Filename of the reference genome used |
| refGenomeBWAFolder | Path to folder containing the reference genome file used by BWA |
| refGenomeBWAPrefix | Common name prefix for all files used by BWA |
| knownSitesFolder   | Path to folder containing the three subsequent files: |
| refSNP             | Reference SNV file used as **known-sites** for BQSR |
| refGoldIndels      | Reference Indel file used as **known-sites** for BQSR |
| refKnownIndels     | Reference Indel file used as **known-sites** for BQSR |
| cdnaIndex          | Index used by Kallisto for pseudo-alignment |
| gtf                | File containing genomic positions (**.gtf**) </br> (see http://www.ensembl.org/info/data/ftp/index.html) |
| proteins           | File containing reference proteins (**.fa**) </br> (see http://www.ensembl.org/info/data/ftp/index.html) |
| STARGenomeDir      | Path to reference genome for usage by STAR |
| vepCacheDir        | Path to cache used by VEP |

## Output

As output a single file: **epitopes.tsv** is generated containing the following columns:

| Parameter             | Type    | Description |
|-----------------------|---------|-------------|
| Peptide               | String  | Sequence of identified neo Epitope |
| Variants              | String  | Comma separated list of mutations that are contained in the sequence. Each mutation is given in the form: </br> **chr\<Id\>:\<Position\>:\<Reference\>:\<Mutated\>** |
| VariantTypes          | String  | Comma separated list of variant types. Each can be one of: </br> **SNV**, **Insertion**, **Deletion**, **Frameshift** |
| VariantCallers        | String  | Comma separated list of number of variant callers that confirmed the respective mutation |
| DNA VAF Normal        | String  | Comma separated list of VAFs in normal DNA sample (Each averaged over the reporting variant callers) |
| DNA VAF Tumor         | String  | Comma separated list of variant VAFs in tumor DNA sample (Each averaged over the reporting variant callers) |
| RNA VAF               | String  | Comma separated list of variant VAFs in tumor RNA |
| Read Count DNA Normal | String  | Comma separated list of variant read counts in normal DNA sample (Each averaged over the reporting variant callers) |
| Read Count DNA Tumor  | String  | Comma separated list of variant read counts in tumor DNA sample (Each averaged over the reporting variant callers) |
| Read Count RNA Tumor  | String  | Comma separated list of variant read counts in tumor RNA sample |
| gene                  | String  | Name of the gene the variants occurred in |
| Protein               | String  | Ensembl-ID of the mutated protein |
| HLA                   | String  | HLA type for which binding affinity is reported |
| BindingCore           | String  | 9-mer directly bound to MHC complex |
| ICore                 | String  | Interaction core |
| RawPredictionScore    | Float   | Raw binding affinity prediction score |
| Affinity (nM)         | Float   | Binding affinity reported in nano molar |
| %Rank                 | Float   | Rank of binding affinity |
| Exp                   | Integer | NetMHCpan associated value |
| Strength              | String  | **WB** if %Rank \< 2, **SB** if %Rank < 0.5, empty otherwise
| TPM                   | Float   | Transcripts per million of the respective gene |
| expressed             | String  | **yes** if TPM \> minTpmCount - **no** otherwise |
| ImmunogenicityScore   | Float   | Raw predicted immunogenicity score ranging from **0** to **1**
| Immunogenic           | String  | **yes** if ImmunogenicityScore \> 0.5 - **no** otherwise 

