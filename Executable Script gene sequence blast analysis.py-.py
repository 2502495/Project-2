#!/usr/bin/env python
# coding: utf-8

# ![resize-17058849241423949201UoBFacultyoflifesciencesRED.png](attachment:d8e5f5e4-0658-44fb-a554-9e6e8c29e56b.png)
# 
# ![UoB Logo](https://github.com/2502495/Project-2/blob/bbbec227de597724a386c12a7307b974d6957d4c/resize-17058847321289857334UoBFacultyoflifesciencesRED.png
# # **PROJECT 2: PATHOGENIC ISLAND DETECTION IN A NEW STRAIN OF _YERSINIA_**
# # Quintessentially Paradigmatic Plan, Designed Based on Knowledge of Various Biological Concepts, Developed and Implemented Different Strategies to Tackle a Scientific Problem Taking Advantage of the Modularity of Bioinformatics Approaches.
# 
# Bioinformatician: Jatin Arora 
# 
# Student ID: 2502495 
# 
# Biomedical Microbiology Laboratory
# 
# Section A
# 
# <centre> ASSESSMENT COVER SHEET </centre>
# 
# |   **STUDENT NUMBER** |  **2502495**    |
# | -------- | ------- |
# | **Programme**  | **MSc Bioinformatics**    |
# | Unit Code and Name | BIOLMM37_2023_TB-1 Group Project    |
# | Project    | 2               |
# | Assessmment Name         |  Group Project Alternative Assessment              |
# | Word Count     | Well under 1500 words                 |
# | Do you give permission for you work to be used anonymously in examples given to students in the future?        |  Yes                |
# | **Title of the Report**        |   **PROJECT 2: PATHOGENIC ISLAND DETECTION IN A NEW STRAIN OF _YERSINIA_**              |
# | Essential UoB Literature Referenced| Workbooks on Accelerated R Couse by Dr Christopher Clements, Statistics Bootcamp, Statistics for Biological Sciences, Year 2 - Quantitative and Computational Methods - 2020 iteration by Multiple Authors and Year 4 - Understanding Data: Experimental Design and Statistics for Life Scientists - 2020 iteration by Multiple Authors on Blackboard (Dr. Christopher Clements, Duncan O’Brien, Duncan Edgley, et al), Genomics Bootcamp by Dr Francesca Segers on Blackboard, Command Line Working Document, and Project Description Powerpoint, by Dr Celine Petitjean on Blackboard|
# | Unit Director         |  Dr Celine Petitjean                  |
# | Date of Submission    | Monday, January 22, 2024               |
# | Place              |  Life Science Building (LSB), School of Biological Sciences, Faculty of Life Sciences, University of Bristol  |
# | Marks based on the summmative computer based assessmment commbininng all the leaarning objectives | 100%. |
# (Malmfors, 2000) 
# 
# By submitting this assignment cover sheet, I confirm that I understand and agree with the following statements:
# 
# • I have not committed plagiarism, cheated or otherwise committed academic misconduct as defined in the University’s Assessment Regulations (available at https://www.bristol.ac.uk/media-library/sites/academic-quality/documents/taught-code/annexes/university-examination-regulations.pdf)
# 
# • I have not submitted this piece, in part or in its entirety, for assessment in another unit assignment (including at other institutions) as outlined in section 4 of the University’s Assessment regulations (available at https://www.bristol.ac.uk/media-library/sites/academic-quality/documents/taught-code/annexes/university-examination-regulations.pdf)
# 
# • I understand that this piece will be scrutinised by anti-plagiarism software and that I may incur penalties if I am found to have committed plagiarism, as outlined in sections 3 of the University’s Examination Regulations (available at https://www.bristol.ac.uk/media-library/sites/academic-quality/documents/taught-code/annexes/university-examination-regulations.pdf)
# 
# This assignment is undertaken by the abovementioned bioinformatician based on the original publication Research Article in Journal of Bacteriology, American Society for Microbiology-Activity of a Holin-Endolysin System in the Insecticidal Pathogenicity Island of *Yersinia enterocolitica* (Katharina, S. 2018)
# Study of a lysis cassette in Tc-PAIY of four ORFs in Yersinia enterocolitica.(Springer et al., 2018)

# **Table of Contents**
# 
# **Section A Assessmment Cover Sheet**............................................................................................ 0
# 
# 0.0 **Assessment**..................................................................................................0.0
# 
# 1.0 **Overview** ................................................................................. 0.0
# 
# 2.0 **Introduction/Background** ....................................................................... 0.0
# 
# 3.0 **Project Planning**.............................................................................................0.0
# 
# 4.0 **Literature Review** ............................................................................................0.0
# 
# 5.0 **Results**...............................................................................0.0
# 
# 6.0 **Further Visualisation and Analysis**...................................................................................0.0
# 
# 7.0 **Applications and Insights** .................................................................................. 0.0
# 
# 8.0 **Real-world Relevance** ....................................................................... 0.0
# 
# 9.0 **Future Research (FOLLOW ON GITHUB FOR MORE)** ..................................................................................................0.0
# 
# 10.0 **Conclusion**.............................................................................................0.0
# 
# **References/Citations**.............................................................................................0.0
# 
# Connor, M. (1992)

# #Accessible via my GitHub Repository https://github.com/2502495/Project-2.git.
# 
# #The format of this Scientific Report is taken from the Skills Team, University of Hull available on University of Bristol 
# #BlackBoard® Day, R.A. (1998), and Teaching and Learning Support (TaLS) – Fact Sheets from University of Bristol Blackboard®.

# In[ ]:


# 0.0 ASSESSMENT
# In tending to the evaluation necessities, the culmination of tasks 1, 2, and 3 for the genome NewStrain1 involved a 
#methodical methodology. Task 1 required the arrangement of the genomic dataset by guaranteeing the right designing of the 
#fasta record, setting the establishment for resulting examinations. Task 2 included executing a BLAST analysis utilizing 
#the BLASTp command-line instrument, questioning the arranged dataset against the Yersinia genome database. 
#The assessment of arrangement scores, E-values, and percent characters during manual analysis affirmed the presence of the 
#four genes of premium from the lysis cassette in Tc-PAIY in NewStrain1. Task 3 was achieved by fostering a Python script 
#for succession recovery, empowering the group to separate explicit quality arrangements in view of their ID numbers. 
#Pushing ahead, Task 4 could be accomplished by planning a complete pipeline and manual to direct the group through the 
#reiteration of investigations on an alternate dataset. This includes data planning, BLAST analysis, manual examination of 
#results, and execution of a script for succession recovery, with the manual giving nitty gritty directions to each step. 
#In a social environment, correspondence, joint effort, and a mutual perspective of tasks and objectives are principal. 
#Powerful designation of obligations, clear documentation, and customary advancement refreshes add to a fruitful and strong 
#collective endeavor. The significance of keeping up with open lines of correspondence and regarding different viewpoints 
#encourages a cooperative climate, guaranteeing the group's aggregate achievement.


# In[4]:


# 1.O OVERVIEW

# The Workflow 

# The whole workflow is devised to effectively work as a modular bioinformatic pipeline with independent but interrelated 
#parts. Installations, steps, and actions needed to be taken are explained here but it is recommended to consult the NCBI 
#BLAST User Manual and the Anaconda Navigator Manual attached (all credits to respective owners)(NCBI BLAST, 2015). 
# Complete Anaconda-NavigatorTM User Manual attached. Documentation hosted by Read the Docs©. Click the following link to 
#access Read the Docs Anaconda-Navigator Manual ((Anaconda, 2018). Created using Sphinx 5.0.2. Built with the PyData Sphinx 
#Theme 0.14.4.
# Also, find attached complete Anaconda-NavigatorTM Tutorial. Documentation hosted by Read the Docs©. Click the following 
#link to access Read the Docs Anaconda-Navigator Tutorial ((Anaconda, 2018). Created using Sphinx 5.0.2. Built with the 
#PyData Sphinx Theme 0.14.4.
# With working directory attached is an Improv Academic Blueprint-A preliminary literature review integrated as part of 
#modular pipeline which defined this pipeline construction, design and build. Please refer to the document for thorough 
#referential enrichmment of scientific and academic concepts governing this workflow.


# 2.0 INTRODUCTION/BACKGROUND
# A detailed analysis of BLAST for the Yersinia project was carried out to investigate the genomic features of a newly recognized strain Yersinia. The analysis particularly concentrated on four genes of interest; possibly associated with the phage holin family – a group of proteins involved in bacterial lysis. The main objectives were to determine if these genes characterize the new strain, how their sequences vary, and how they can be compared against sequences from different bacterial strains, thus shedding light on their functional significance and evolutionary implications
# ![Picture 1.png](attachment:ded1cc1f-5607-471d-8ddd-d5d3571a9541.png)
# ![IMAGE] (https://github.com/2502495/Project-2/blob/a3c456a018655b38858e66ff15edd5e74a66a7e3/Picture%201.png)
# Figure 1. Pathogenecity Island in *Yersinia enterocolitica* W22703
# Fuchs, T.M., Bresolin, G., Marcinowski, L. et al., BMC Microbiol 8, 214 (2008).

# In[2]:


# A detailed analysis of BLAST for the Yersinia project was carried out to investigate the genomic features of a 
# newly recognized strain Yersinia. The analysis particularly concentrated on four genes of interest; 
# possibly associated with the phage holin family – a group of proteins involved in bacterial lysis. 
# The main objectives were to determine if these genes characterize the new strain, how their sequences vary, 
# and how they can be compared against sequences from different bacterial strains, thus shedding light on their functional 
# significance and evolutionary implications.


# Holin:
# Holin, a transmembrane protein with 1-4 transmembrane sequences, exhibits diverse classifications based on topological, 
# phylogenetic, and motif analyses, resulting in the identification of 58 families and 7 superfamilies. Experimental evidence 
# from previous cell killing assays indicates that the presence of both the C-terminal region and transmembrane domains is 
# essential for the generation of pores in the inner cell membrane.

# Endolysin:
# Endolysins, enzymatic agents pivotal in bacterial cell lysis, possess a structural composition comprising distinct 
# elements: catalytic domain, cell wall binding domain, linker regions, and accessory domains. Their primary 
# function involves the targeted degradation of the bacterial cell walls peptidoglycan, facilitating the release of phage 
# progeny. Notably, the N-terminal cell wall binding domain recognizes specific bacterial components, while the 
# C-terminal catalytic domain executes enzymatic activity.

# i-spanin and o-spanin:
# Encoded by the RzRz1 genes in both lambda and φ21 phages, i-spanin and o-spanin represent proteins integral to the 
# spanin complex. The proposed model for spanin function elucidates its role in the removal of the outer membrane barrier 
# by fusing it to the inner membrane. Interaction between i-spanin and o-spanin occurs through their C-terminal regions, 
# wherein conformational changes in the two alpha-helix regions of i-spanin facilitate the fusion of the two membranes, 
#ultimately leading to the formation of pores on the bacterial membrane.


# In[ ]:


#3.0 PROJECT PLANNING
#To smooth out the analysis interaction and empower the group to rehash the examinations on an alternate dataset, 
#a pipeline alongside a complete manual can be given. Here is an improved version of pipeline:

#Pipeline:
#1.	Step 1: Data Preparation
#•	Input: Fasta file of the new Yersinia strain genome.
#•	Output: Prepared dataset.
#2.	Step 2: BLAST Analysis
#•	Input: Prepared dataset and known protein sequences.
#3•	Output: BLAST results file.
#3.	Step 3: BLAST Results Analysis
#•	Input: BLAST results file.
#•	Output: Confirmation of gene presence.
#4.	Step 4: Script Execution
#•	Input: Fasta file, target gene ID.
#•	Output: Retrieved gene sequence.

#Manual
#•	Detailed instructions for each step.
#•	Necessary parameters and commands for tools/scripts.
#•	Troubleshooting guidance.

#Explanation of Findings
#•	Clarification of the significance of alignment scores, E-values, and percent identities.
#•	Confirmation or rejection of the presence of genes based on the analysis.
#•	Interpretation of any unexpected results or discrepancies.


# 1. Know if the four genes of the lysis cassette in Tc-PAIY are present in their new strain
# 
# The team expects to affirm the presence of the four genes from the lysis cassette in Tc-PAIY inside their new strain of Yersinia, signified as NewStrain1. Following the execution of a BLAST analysis utilizing the BLASTp command-line device, questioning NewStrain1 against the Yersinia genome database, the group cautiously inspects the got results. The analysis centers around the alignment scores, E-values, and percent characters of the matches to the realized protein successions from the lysis cassette. Upon intensive examination, the group sees that the alignment scores are quite high, showing hearty likeness. The E-values are low, recommending a high statistical meaning of the matches, and the percent characters confirm a substantial comparability between the inquiry and database groupings. These discoveries by and large affirm the presence of the four genes of interest in the lysis cassette in Tc-PAIY inside the genomic cosmetics of NewStrain1.
# 

# In[ ]:


#Import Required Libraries
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import os


# In[ ]:


from Bio import Entrez

# Set your email here
Entrez.email = "your_email@example.com"

# List of provided protein IDs
ids = [
    "CAI77377.1",
    "AAT90758.1",
    "WP_004875774.1",
    "WP_050116509.1",
    "AAM83780.1",
    "ABX88260.1",
    "CAH22792.1",
    "CAE13263.1",
    "AAG55972.1",
    "WP_047345313.1",
    "WP_047964150.1",
    "AAQ61208.1"
]

# Open a file to write the sequences
# Running the Script

    #Place the FASTA File: Ensure the holin_sequences.fasta file is in the same directory as the script.
    #Open Command Line: Use Terminal (Mac/Linux) or Command Prompt (Windows).
    #Navigate to Script Directory: Change directory to where the script is located using cd path/to/script.
    #Run the Script: Execute the script by typing python script_name.py (replace script_name.py with the actual script filename).
    
with open("holin_sequences.fasta", "w") as file:
    for i in ids:
        # Fetch each sequence from NCBI and write to the file
        handle = Entrez.efetch(db="protein", id=i, rettype="fasta", retmode="text")
        file.write(handle.read())


# 2. Know how conserved these genes are compared to other Yersinia strain
# 
# In Analysis to decide the protection of the four genes from the lysis cassette in Tc-PAIY inside the new strain NewStrain1 contrasted with other Yersinia strains, the group would broaden their examination. Following the BLAST analysis, the center movements to surveying the likeness of the recognized genes across a scope of Yersinia strains. The group looks at the alignment scores, E-values, and percent characters of the matched genes in NewStrain1 against those in different Yersinia strains recorded in the database. On the off chance that the alignment scores are reliably high and the E-values are reliably low across various strains, it proposes a serious level of protection. Also, noticing high percent personalities would additionally uphold the thought that the genes are moderated among various Yersinia strains. These discoveries would demonstrate that the distinguished genes from the lysis cassette in Tc-PAIY are available as well as profoundly monitored in NewStrain1 contrasted with other Yersinia strains, suggesting utilitarian importance across various genomic foundations.
# 

# In[ ]:


3. Be able to retrieve their specific protein sequence from the results (or any other in the future)

#The research group imagines a deliberate pipeline intended to investigate extra genes of premium inside the genomic scene 
#of the new Yersinia strain, NewStrain1. The initial step includes setting up the genomic dataset, guaranteeing appropriate 
#designing to work with downstream investigations. Consequently, the group utilizes the BLASTp command-line instrument to 
#direct an underlying analysis against an exhaustive quality database, producing a rundown of potential genes of interest. 
#To refine and focus on this rundown, a separating cycle is carried out, taking into account alignment scores, E-values, 
#and #other significant models. The third step coordinates practical explanation tools, like InterProScan, to clarify the 
#potential jobs and elements of the recognized up-and-comer genes. The refined rundown is additionally approved through a 
#complete writing survey or computerized tools for writing mining, guaranteeing that the chose genes have known capabilities 
#or affiliations. As a closing step, a discretionary arrangement for experimental approval is incorporated, considering 
#further affirmation of the distinguished genes through research center investigations or joint efforts with experimental 
#biologists. This comprehensive pipeline offers an organized and precise methodology for the group to investigate, 
#distinguish, and examine extra genes of premium in NewStrain1's genomic setting.


# In[ ]:


#Function to Retrieve Sequences from a FASTA File

def retrieve_sequence(fasta_file, gene_id):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == gene_id:
            return record.seq
    return None


# In[ ]:


A. #Function to Perform BLAST Analysis

def blast_sequence(sequence, database="nr", program="blastp", format_type="Text"):
    result_handle = NCBIWWW.qblast(program, database, sequence, format_type=format_type)
    return result_handle.read()



# In[ ]:


# Main Workflow
def main():
    # Define the path to your FASTA file and gene IDs
    fasta_file = "holin_sequences.fasta"
    gene_ids = [
    "CAI77377.1",
    "AAT90758.1",
    "WP_004875774.1",
    "WP_050116509.1",
    "AAM83780.1",
    "ABX88260.1",
    "CAH22792.1",
    "CAE13263.1",
    "AAG55972.1",
    "WP_047345313.1",
    "WP_047964150.1",
    "AAQ61208.1"
]
    
    # Loop through each gene ID
    for gene_id in gene_ids:
        sequence = retrieve_sequence(fasta_file, gene_id)
        if sequence:
            print(f"Sequence for {gene_id}: {sequence}")
            # Perform BLAST analysis
            blast_result = blast_sequence(sequence)
            # Save BLAST results to a file
            with open(f"{gene_id}_blast_results.txt", "w") as result_file:
                result_file.write(blast_result)
            print(f"BLAST results for {gene_id} saved.")
        else:
            print(f"Sequence for {gene_id} not found in {fasta_file}")

if __name__ == "__main__":
    main()


# In[ ]:


# Output file for consolidated results
fasta_file = "holin_sequences.fasta"
gene_ids = [
    "CAI77377.1",
    "AAT90758.1",
    "WP_004875774.1",
    "WP_050116509.1",
    "AAM83780.1",
    "ABX88260.1",
    "CAH22792.1",
    "CAE13263.1",
    "AAG55972.1",
    "WP_047345313.1",
    "WP_047964150.1",
    "AAQ61208.1"
]
with open("consolidated_blast_results.txt", "w") as consolidated_file:
    # Loop through each gene ID
    for gene_id in gene_ids:
        sequence = retrieve_sequence(fasta_file, gene_id)
        if sequence:
            # Write sequence to the consolidated file
            consolidated_file.write(f"Sequence for {gene_id}:\n{sequence}\n")
            
            # Perform BLAST analysis
            blast_result = blast_sequence(sequence)
            
            # Save BLAST results to individual files
            with open(f"{gene_id}_blast_results.txt", "w") as result_file:
                result_file.write(blast_result)
            
            # Write BLAST results to the consolidated file
            consolidated_file.write(f"BLAST results for {gene_id}:\n{blast_result}\n")
            
            print(f"BLAST results for {gene_id} saved.")
        else:
            print(f"Sequence for {gene_id} not found in {fasta_file}")


# In[4]:


# 4.0 IMPROMPTU LITERATURE REVIEW 

#The analysis identified the presence of homologous genes to phage holin family proteins, notably through proteins 
#“WP_005163187.1” and “MBO0362545.1”. These proteins had high alignment scores and low E-values in the BLAST results, 
#meaning that they closely matched the genes of interest in this new strain of Yersinia. This visionary evidence also 
#indicates that these phage holin-related genes have been incorporated into the genome of the new strain, which might bear 
#their influence on the biological characteristics of this new strain such as its pathogenicity.


# In[ ]:


#5.0 RESULTS

#BLAST Results Deciphering the Match Between Query vs the Database
#The Basic Local Alignment Search Tool (BLAST) offers a powerful means to assess the relationship between a query sequence 
#and a vast protein database. By analyzing the alignment details, researchers can gain valuable insights into the degree of 
#similarity, evolutionary history, and functional potential of the queried sequence. Several key parameters within the 
#BLAST results facilitate this interpretation:

#Sequence Similarity

#1. Percentage Identity: This metric quantifies the percentage of identical residues between the query and subject sequences 
#within the aligned region. Higher values indicate a closer relationship, suggesting potential functional conservation or 
#shared ancestry.
#2. Alignment Length: This parameter represents the total number of residues included in the alignment, encompassing both 
#matches and gaps. A longer alignment signifies a more extensive region of similarity, potentially hinting at a shared 
#functional domain or structural motif.

#Alignment Features:

#3. Mismatches: This metric indicates the number of positions within the aligned region where the query and subject sequences #differ. While mismatches can arise due to mutations, their abundance can also point towards functional divergence or #convergent evolution.
#4. Gap Opens: This parameter denotes the number of gaps introduced within the alignment to accommodate insertions or 
#deletions in either the query or subject sequence. Frequent gaps can suggest lineage-specific adaptations or insertions of 
#functional elements.

#Positional Information:

#5. Query Start/End & Subject Start/End: These parameters define the precise locations within the query and subject 
#sequences where the aligned region begins and terminates. This information is crucial for mapping functional domains or 
#conserved motifs identified within the alignment back to their respective sequences.

#Statistical Significance:

#6. E-Value: (Expect Value) This parameter estimates the number of alignments with equal or better scores that are expected 
#to occur by chance in a random sequence database. A lower E-Value indicates a more statistically significant match, 
#suggesting a less likely occurrence by chance.
#Bit Score: This statistical score, calculated by the BLAST algorithm, reflects the quality and strength of the alignment. 
#A higher bit score signifies a more robust and informative alignment, suggesting a closer evolutionary relationship or 
#functional similarity between the query and subject sequences.
#By carefully analyzing these parameters in conjunction, researchers can decipher the intricate relationship between the 
#queried amino acid sequence and its potential matches within the protein database. This detailed understanding facilitates 
#hypothesis generation, functional prediction, and ultimately, advances our knowledge of protein evolution and function.

Analysis of BLAST Results for NewYersinia Strains: Implications for Pathogenicity

NewYersinia1:

Analysis of the BLAST results for NewYersinia1 (Figure 3) reveals a high degree of sequence similarity between the query 
sequence and all matched sequences. This is evidenced by:

High identity: All matched sequences exhibit above 80% identity, ranging from 81.513% to 93.421%. This indicates a strong 
conservation of amino acid residues, suggesting a close evolutionary relationship and potentially shared functional 
properties.
Significant e-values: The e-values of all matches are below 1.41E-48, indicating a very low probability of finding an 
equally good alignment by chance. This further strengthens the confidence in the identified sequences as true homologs.
High bit scores: Bit scores exceeding 145 in all matches suggest a statistically significant and robust alignment, further 
supporting the close relationship between the query and subject sequences.
Based on these parameters, the BLAST results for NewYersinia1 strongly suggest that the analyzed strain possesses a high 
degree of similarity to known pathogenic Yersinia strains. This finding warrants further investigation into the potential 
virulence factors and pathogenicity determinants present in NewYersinia1.

get_ipython().system('[NY1_BLAST result.png](attachment:cb02d122-bffb-4d50-9d2d-7a5421d6a0de.png)')

get_ipython().system('[BLAST RESULTS IMAGE] (https://github.com/2502495/Project-2/blob/f040eb181a3f06670f49d260b2ea001cb27a58a3/NY1_BLAST%20result.png)')

Figure 2. NewYersinia2: W22703 Yersinia1 BLAAST Results

The BLAST results for NewYersinia2 (Figure 4) present a contrasting picture. Only three matched sequences were identified, 
corresponding to Holin, Endolysin, and i-spanin proteins. Notably:

Partial matches: The identities of these matched sequences range from 54.622% to 78.195%, indicating only partial similarity 
to the query sequences. This suggests potential divergence or functional modifications within these proteins compared to 
known Yersinia strains.
Significant e-values: Similar to NewYersinia1, the e-values of all matches are below 8.21E-41, suggesting a low probability 
of chance occurrence. This reinforces the confidence in the identified sequences as relevant homologs, despite their 
partial similarity.
High bit scores: Bit scores exceeding 128 for all matches indicate statistically significant alignments, further supporting 
the identified matches as true homologs.
The absence of a matched sequence resembling o-spanin in the NewYersinia2 protein database is a noteworthy observation. 
O-spanin plays a crucial role in the Holin-Endolysin system, facilitating cell wall lysis during bacteriophage infection. 
The lack of a homologous sequence suggests a potential dysfunction or modification within this critical system, potentially 
impacting the ability of NewYersinia2 to effectively lyse host cells.

get_ipython().system('[NY2_BLAST result.png](attachment:71785f1f-b680-495c-907d-dfda3cb10643.png)')

get_ipython().system('[BLAST RESULTS IMAGE] (https://github.com/2502495/Project-2/blob/f040eb181a3f06670f49d260b2ea001cb27a58a3/NY2_BLAST%20result.png)')

Figure 3. W22703 Yersinia2 BLAST Results

The contrasting BLAST results for NewYersinia1 and NewYersinia2 suggest distinct evolutionary trajectories and potential 
differences in their pathogenic potential. NewYersinia1 exhibits high sequence similarity to known pathogenic strains, 
suggesting a strong likelihood of virulence. Conversely, NewYersinia2 displays partial matches for key proteins involved 
in bacteriophage infection, along with a notable absence of o-spanin, hinting at a potential dysfunction in the 
Holin-Endolysin system and potentially reduced virulence capabilities. Further investigation is warranted to elucidate the 
specific virulence factors and host-pathogen interactions associated with each strain.

The analysis identified the presence of homologous genes to phage holin family proteins, notably through proteins “WP_005163187.1” and “MBO0362545.1”. These proteins had high alignment scores and low E-values in the BLAST results, meaning that they closely matched the genes of interest in this new strain of Yersinia. This visionary evidence also indicates that these phage holin-related genes have been incorporated into the genome of the new strain, which might bear their influence on the biological characteristics of this new strain such as its pathogenicity.
One of the notable findings from BLAST results was the high sequence conservation of these phage holin-related genes. Several sequences from different Yersinia species, relative to the query protein “CAI77377.1” show a high level of similarity. This statement emphasizes a significant evolutionary conservation, implying that these genes possess vital biological functions and are presumably preserved because they are important to the survival, and pathogenicity of Yersinia.


# In[ ]:


#5.0 FURTHER VISUALISATION AND ANALYSIS
#4. To understand the different outputs from BLAST and how to analyse them

#In the situation where the research group tries to fathom the different results created by BLAST (Basic Local Alignment 
#Search Apparatus) and foster a comprehension of how to break down them, a precise methodology is framed.

#1.1 Perform BLAST Analysis
#To start the cycle, the group chooses the DNA or protein successions of interest, explicitly those beginning from the 
#four genes inside the lysis cassette in Tc-PAIY, to act as queries in the BLAST analysis. The decision of the BLAST 
#apparatus is dependent upon the idea of the arrangements, whether they are nucleotide or protein, and sticks to explicit 
#database necessities.


#1.2 Examine BLAST Output
#Upon execution, the BLAST device produces a result document outfitting alignments between the question arrangements and 
#those present in the chose database. The resultant data normally incorporates a thorough rundown of matches, including 
#alignment scores, E-values meaning statistical importance, and other relevant data fundamental for additional analysis.

#1.3 Interpret Alignment Scores
#In this situation, the research group digs into the assessment of alignment scores. Higher scores are proposed to 
#connote more grounded likenesses between the inquiry successions and the hits inside the database, giving a 
#quantitative proportion of alignment quality.

#1.4 Evaluate E-values
#The statistical meaning of the alignments is surveyed through E-values, addressing the expected number of irregular 
#hits having a comparative or unrivaled score. Lower E-values are focused on, demonstrating higher trust in the noticed 
#similitudes.

#1.5 Assess Sequence Identities and Similarities
#The group examines the level of character and similitude inside the BLAST results. These measurements explain the level of 
#resemblance between the question successions and the matched database groupings, assuming a significant part in grasping 
#grouping preservation.

#Visualize Alignments
#To upgrade understanding, the group might decide to utilize perception tools or programming, working with the graphical 
#portrayal of alignments. This perception supports clarifying explicit locales of comparability and disparity inside the 
#successions.

#Filter and Refine Results
#Contingent upon the targets of the research, the group takes part in the filtration and refinement of BLAST results in 
#light of foreordained standards. Boundaries like alignment length, percent character, or E-esteem limits are potentially 
#utilized to distil the most appropriate data.

#Compare Multiple Hits
#In getting different hits, the research group embraces a near analysis. Different boundaries are considered to focus on 
#and select the most pertinent matches among the plenty of results.

#Document and Report Findings
#At long last, the group participates in fastidious documentation and announcing of their discoveries. This complete record 
#incorporates chosen hits, nitty gritty alignments, and relevant measurements, guaranteeing straightforwardness and giving 
#a primary asset to future reference and analysis.


# In[ ]:


# Execution of the following Python script will perform data analysis and visualisation
{

"cells": \[

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"# Exercise 4.1: Pathogenicity islands\n",

"\n",

"This exercise was inspired by \[Libeskind-Hadas and Bush, \*Computing
for Biologists\*, Cambridge University Press,
2014\](https://www.cs.hmc.edu/CFB).\n",

"\n",

"\<hr\>"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"For this and the \[next problem\](exercise_4.2.ipynb), we will work
with real data from the \*Yersinia Enterocolitica\* genome. The section
of the genome we will work with is in the file
\`\~git/bootcamp/data/Yersinia_spi1_region.fna\`. I cut it out of the
\[full genome\](http://www.ncbi.nlm.nih.gov/nucleotide/821161554). It
contains \*Yersinia\* pathogenicity island I (SPI1), which contains
genes for surface receptors for host-pathogen interactions.\n",

"\n",

"Pathogenicity islands are often marked by different GC content than the
rest of the genome. We will try to locate the pathogenicity island(s) in
our section of the \*Yersinia\* genome by computing GC content.\n",

"\n",

"\*\*a)\*\* Use principles of TDD to write a function that divides a
sequence into blocks and computes the GC content for each block,
returning a tuple. The function signature should look like\n",

"\n",

" gc_blocks(seq, block_size)\n",

" \n",

"To be clear, if \`seq = 'ATGACTACGT'\` and \`block_size = 4\`, the
blocks to be considered are\n",

"\n",

" ATGA\n",

" CTAC\n",

" \n",

"and the function should return \`(0.25, 0.5)\`. Note that the blocks
are non-overlapping and that we don't bother with the fact that end of
the sequence that does not fit completely in a block."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*b)\*\* Write a function that takes as input a sequence, block size,
and a threshold GC content, and returns the original sequence where
every base in a block with GC content above threshold is capitalized and
every base below the threshold is lowercase. You would call the function
like this:\n",

"\n",

" mapped_seq = gc_map(seq, block_size, gc_thresh)\n",

"\n",

"For example, \n",

"\n",

" gc_map('ATGACTACGT', 4, 0.4)\n",

"\n",

"returns \`'atgaCTAC'\`. Note that bases not included in GC blocks (in
this case the \`GT\` at the end of the sequence) are not included in the
output sequence. Again, use principles of TDD."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*c)\*\* Use the \`gc_map()\` function to generate a GC content map
for the \*Yersinia\* sequence with \`block_size = 1000\` and \`gc_thresh
= 0.45\`. Where do you think the pathogenicity island is?"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*d)\*\* Write the GC-mapped sequence (with upper and lower
characters) to a new FASTA file. Use the same description line (which
began with a \`\>\` in the original FASTA file), and have line breaks
every 60 characters in the sequence."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\<br /\>"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"## Solution\n"

\]

},

{

"cell_type": "code",

"execution_count": 1,

"metadata": {},

"outputs": \[\],

"source": \[

"import numpy as np"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*a)\*\* First let's write our tests. In writing the tests, I am
making the design decision that I will count the characters \`G\`,
\`g\`, \`C\`, and \`c\` as contributing to GC content, and that I will
not check to make sure the sequence is valid. I also make the design
decision that an empty sequence has zero GC content."

\]

},

{

"cell_type": "code",

"execution_count": 2,

"metadata": {},

"outputs": \[\],

"source": \[

"def test_gc_content():\n",

" assert gc_content('') == 0.0\n",

" assert gc_content('G') == 1.0\n",

" assert gc_content('g') == 1.0\n",

" assert gc_content('C') == 1.0\n",

" assert gc_content('c') == 1.0\n",

" assert gc_content('gcgcgc') == 1.0\n",

" assert gc_content('aaatatata') == 0.0\n",

" assert np.isclose(gc_content('ggatcggcga'), 0.7)\n",

" assert np.isclose(gc_content('attgggggcaatta'), 3/7)"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"The function is fairly simple. We loop through the sequence with a
stride equal to the block size, computing the GC content for each
subsequence of that length. We start with a function to compute GC
content for a sequence."

\]

},

{

"cell_type": "code",

"execution_count": 3,

"metadata": {},

"outputs": \[\],

"source": \[

"def gc_content(seq):\n",

" \\"\\"\\"GC content of a given sequence\\"\\"\\"\n",

" if seq == '':\n",

" return 0.0\n",

" \n",

" seq = seq.upper()\n",

" return (seq.count('G') + seq.count('C')) / len(seq)"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Now let's test it."

\]

},

{

"cell_type": "code",

"execution_count": 4,

"metadata": {},

"outputs": \[\],

"source": \[

"test_gc_content()"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Passage! Next, we write the looping function, starting with its tests."

\]

},

{

"cell_type": "code",

"execution_count": 5,

"metadata": {},

"outputs": \[\],

"source": \[

"def test_gc_blocks():\n",

" assert gc_blocks('', 10) == tuple()\n",

" assert gc_blocks('gcgcgcgcg', 10) == tuple()\n",

" assert gc_blocks('gcgcgcg', 4) == (1.0,)\n",

" assert gc_blocks('gcgcgcgc', 4) == (1.0, 1.0)\n",

" assert gc_blocks('gcgcgcgcat', 4) == (1.0, 1.0)\n",

"\n",

" test_tuple = gc_blocks('gcgagcgcat', 4)\n",

" assert np.isclose(test_tuple\[0\], 0.75) and test_tuple\[1\] ==
1.0\n",

" \n",

" assert gc_blocks('gcgtagagc', 1) == (1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
0.0, 1.0, 1.0)\n",

" assert gc_blocks('gcgtagagc', 2) == (1.0, 0.5, 0.5, 0.5)\n",

" assert np.isclose(gc_blocks('gcgtagagc', 3), (1.0, 1/3, 2/3)).sum() ==
3"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Now let's write our function."

\]

},

{

"cell_type": "code",

"execution_count": 6,

"metadata": {},

"outputs": \[\],

"source": \[

"def gc_blocks(seq, block_size):\n",

" \\"\\"\\"\n",

" Divide sequence into non-overlapping blocks\n",

" and compute GC content of each block.\n",

" \\"\\"\\"\n",

" blocks = \[\]\n",

" for i in range(0, len(seq) - (len(seq) % block_size), block_size):\n",

" blocks.append(gc_content(seq\[i:i+block_size\]))\n",

" return tuple(blocks)"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"And the tests...."

\]

},

{

"cell_type": "code",

"execution_count": 7,

"metadata": {},

"outputs": \[\],

"source": \[

"test_gc_blocks()"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Success! Let's take this function for a spin, looking at 1000-base
blocks. We will use the FASTA reader function from a previous exercise
to read in the \_Yersinia\_ genome fragment."

\]

},

{

"cell_type": "code",

"execution_count": 8,

"metadata": {},

"outputs": \[

{

"data": {

"text/plain": \[

"(0.521,\n",

" 0.556,\n",

" 0.54,\n",

" 0.498,\n",

" 0.551,\n",

" 0.508,\n",

" 0.563,\n",

" 0.484,\n",

" 0.58,\n",

" 0.557,\n",

" 0.523,\n",

" 0.524,\n",

" 0.621,\n",

" 0.556,\n",

" 0.481,\n",

" 0.57,\n",

" 0.581,\n",

" 0.614,\n",

" 0.603,\n",

" 0.526,\n",

" 0.524,\n",

" 0.591,\n",

" 0.563,\n",

" 0.596,\n",

" 0.563,\n",

" 0.6,\n",

" 0.613,\n",

" 0.594,\n",

" 0.486,\n",

" 0.554,\n",

" 0.566,\n",

" 0.592,\n",

" 0.563,\n",

" 0.537,\n",

" 0.575,\n",

" 0.501,\n",

" 0.54,\n",

" 0.555,\n",

" 0.487,\n",

" 0.416,\n",

" 0.423,\n",

" 0.371,\n",

" 0.394,\n",

" 0.48,\n",

" 0.454,\n",

" 0.474,\n",

" 0.434,\n",

" 0.396,\n",

" 0.37,\n",

" 0.456,\n",

" 0.409,\n",

" 0.457,\n",

" 0.4,\n",

" 0.405,\n",

" 0.475,\n",

" 0.47,\n",

" 0.479,\n",

" 0.494,\n",

" 0.497,\n",

" 0.516,\n",

" 0.444,\n",

" 0.433,\n",

" 0.471,\n",

" 0.458,\n",

" 0.53,\n",

" 0.458,\n",

" 0.56,\n",

" 0.427,\n",

" 0.47,\n",

" 0.438,\n",

" 0.465,\n",

" 0.473,\n",

" 0.46,\n",

" 0.399,\n",

" 0.426,\n",

" 0.359,\n",

" 0.469,\n",

" 0.433,\n",

" 0.425,\n",

" 0.504,\n",

" 0.578,\n",

" 0.576,\n",

" 0.553,\n",

" 0.531,\n",

" 0.57,\n",

" 0.599,\n",

" 0.562,\n",

" 0.555,\n",

" 0.595,\n",

" 0.586,\n",

" 0.55,\n",

" 0.56,\n",

" 0.545,\n",

" 0.553,\n",

" 0.537,\n",

" 0.519,\n",

" 0.519,\n",

" 0.567,\n",

" 0.551,\n",

" 0.548,\n",

" 0.559,\n",

" 0.527,\n",

" 0.559,\n",

" 0.529,\n",

" 0.49,\n",

" 0.533,\n",

" 0.58,\n",

" 0.545,\n",

" 0.558,\n",

" 0.575,\n",

" 0.555,\n",

" 0.49,\n",

" 0.567,\n",

" 0.515,\n",

" 0.518,\n",

" 0.485,\n",

" 0.38,\n",

" 0.461,\n",

" 0.568,\n",

" 0.575,\n",

" 0.567,\n",

" 0.57,\n",

" 0.472,\n",

" 0.513,\n",

" 0.582,\n",

" 0.476,\n",

" 0.505,\n",

" 0.524,\n",

" 0.51,\n",

" 0.512,\n",

" 0.391,\n",

" 0.463,\n",

" 0.57,\n",

" 0.546,\n",

" 0.535,\n",

" 0.525,\n",

" 0.525,\n",

" 0.529,\n",

" 0.58,\n",

" 0.555,\n",

" 0.558,\n",

" 0.563,\n",

" 0.525,\n",

" 0.505,\n",

" 0.557,\n",

" 0.554,\n",

" 0.484,\n",

" 0.525,\n",

" 0.567,\n",

" 0.467,\n",

" 0.527,\n",

" 0.55,\n",

" 0.577,\n",

" 0.554,\n",

" 0.538,\n",

" 0.429,\n",

" 0.507,\n",

" 0.557,\n",

" 0.592,\n",

" 0.595,\n",

" 0.554,\n",

" 0.521,\n",

" 0.539,\n",

" 0.521,\n",

" 0.45,\n",

" 0.608,\n",

" 0.489,\n",

" 0.477,\n",

" 0.552,\n",

" 0.508,\n",

" 0.544,\n",

" 0.495,\n",

" 0.543,\n",

" 0.56,\n",

" 0.596,\n",

" 0.547,\n",

" 0.581,\n",

" 0.548,\n",

" 0.537,\n",

" 0.529,\n",

" 0.513,\n",

" 0.499,\n",

" 0.545,\n",

" 0.567,\n",

" 0.52,\n",

" 0.545,\n",

" 0.548,\n",

" 0.522,\n",

" 0.533,\n",

" 0.558,\n",

" 0.586,\n",

" 0.469,\n",

" 0.516,\n",

" 0.509,\n",

" 0.511,\n",

" 0.569,\n",

" 0.575,\n",

" 0.559,\n",

" 0.545,\n",

" 0.502)"

\]

},

"execution_count": 8,

"metadata": {},

"output_type": "execute_result"

}

\],

"source": \[

"def read_fasta(filename):\n",

" \\"\\"\\"Read a sequence in from a FASTA file containing a single
sequence.\n",

" \n",

" We assume that the first line of the file is the descriptor and
all\n",

" subsequent lines are sequence. \n",

" \\"\\"\\"\n",

" with open(filename, 'r') as f:\n",

" \# Read in descriptor\n",

" descriptor = f.readline().rstrip()\n",

"\n",

" \# Read in sequence, stripping the whitespace from each line\n",

" seq = ''\n",

" line = f.readline().rstrip()\n",

" while line != '':\n",

" seq += line\n",

" line = f.readline().rstrip()\n",

" \n",

" return descriptor, seq\n",

"\n",

"\n",

"descriptor, seq = read_fasta('data/Yersinia_spi1_region.fna')\n",

"\n",

"gc_blocks(seq, 1000)"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"We get a tuple of GC content, which is hard to look at on screen, but
this is useful for plotting GC content over the course of a sequence. We
will learn how to plot later in the bootcamp."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*b)\*\* We just use our already-written \`gc_content()\` function to
decide how to modify the string of the sequence. First, the tests. We
make the design decision that we will truncate the sequence if the
3'-most end is shorter than the block length."

\]

},

{

"cell_type": "code",

"execution_count": 9,

"metadata": {},

"outputs": \[\],

"source": \[

"def test_gc_map():\n",

" assert gc_map('', 10, 0.5) == ''\n",

" assert gc_map('ATATATATA', 4, 0.5) == 'atatatat'\n",

" assert gc_map('GCGCGCGCG', 4, 0.5) == 'GCGCGCGC'\n",

" assert gc_map('GATCGATCC', 4, 0.5) == 'GATCGATC'\n",

" assert gc_map('GATCGATCC', 4, 0.51) == 'gatcgatc'\n",

" assert gc_map('GATCGATCC', 3, 0.5) == 'gatCGATCC'\n",

" assert gc_map('GATCGATCC', 3, 0.75) == 'gatcgatcc'\n",

" assert gc_map('GATCGATCC', 3, 0.25) == 'GATCGATCC'"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Now the function...."

\]

},

{

"cell_type": "code",

"execution_count": 10,

"metadata": {},

"outputs": \[\],

"source": \[

"def gc_map(seq, block_size, gc_thresh):\n",

" \\"\\"\\"Give back seq with lowercase letters where GC content is
low.\\"\\"\\" \n",

" out_seq = ''\n",

"\n",

" \# Determine GC content of each block and change string
accordingly\n",

" for i in range(0, len(seq) - (len(seq) % block_size), block_size):\n",

" if gc_content(seq\[i:i+block_size\]) \< gc_thresh:\n",

" out_seq += seq\[i:i+block_size\].lower()\n",

" else:\n",

" out_seq += seq\[i:i+block_size\].upper()\n",

"\n",

" return out_seq"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"And the tests."

\]

},

{

"cell_type": "code",

"execution_count": 11,

"metadata": {},

"outputs": \[\],

"source": \[

"test_gc_map()"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Passage! We can now use these functions to analyze sequences of
interest."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*c)\*\* Let's do it for \*Yersinia\*!"

\]

},

{

"cell_type": "code",

"execution_count": 12,

"metadata": {},

"outputs": \[\],

"source": \[

"sal_gcmap = gc_map(seq, 1000, 0.45)"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"To save on display space, we will not display the sequence here.
Scrolling through the GC map file generated in the next part, the
pathogenicity island appears to occur about a quarter of the way into
this subsequence."

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"\*\*d)\*\* To write the file out, we use the fact that we conveniently
kept the description text when we parsed the \*Yersinia\* FASTA file in
the first place. We then just write the \`sal_gcmap\` string in blocks
of 60. We have to make sure to get the last few bases as well."

\]

},

{

"cell_type": "code",

"execution_count": 13,

"metadata": {},

"outputs": \[\],

"source": \[

"# Write the result\n",

"with open('Yersinia_spi1_region_gc_map.fna', 'w') as f:\n",

" \# Write description text\n",

" f.write(descriptor + '\\\n')\n",

"\n",

" \# Write sequence in blocks of 60\n",

" i = 0\n",

" while i \< len(sal_gcmap) - 59:\n",

" f.write(sal_gcmap\[i:i+60\] + '\\\n')\n",

" i += 60\n",

" \n",

" \# Write last line\n",

" f.write(sal_gcmap\[i:\] + '\\\n')"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"We'll take a quick look to see it worked out ok."

\]

},

{

"cell_type": "code",

"execution_count": 14,

"metadata": {},

"outputs": \[

{

"name": "stdout",

"output_type": "stream",

"text": \[

"\>gi\|821161554\|gb\|CP011428.1\| Yersinia Enterocolitica subsp.
enterica strain YU39, complete genome, subsequence 3000000 to
3200000\n",

"AAAACCTTAGTAACTGGACTGCTGGGATTTTTCAGCCTGGATACGCTGGTAGATCTCTTC\n",

"ACGATGGACAGAAACTTCTTTCGGGGCGTTCACGCCAATACGCACCTGGTTGCCCTTCAC\n",

"CCCTAAAACTGTCACGGTGACCTCATCGCCAATCATGAGGGTCTCACCAACTCGACGAGT\n",

"CAGAATCAGCATTCTTTGCTCCTTGAAAGATTAAAAGAGTCGGGTCTCTCTGTATCCCGG\n",

"CATTATCCATCATATAACGCCAAAAAGTAAGCGATGACAAACACCTTAGGTGTAAGCAGT\n",

"CATGGCATTACATTCTGTTAAACCTAAGTTTAGCCGATATACAAAACTTCAACCTGACTT\n",

"TATCGTTGTCGATAGCGTTGACGTAAACGCCGCAGCACGGGCTGCGGCGCCAACGAACGC\n",

"TTATAATTATTGCAATTTTGCGCTGACCCAGCCTTGTACACTGGCTAACGCTGCAGGCAG\n",

"AGCTGCCGCATCCGTACCACCGGCTTGCGCCATGTCCGGACGACCGCCACCCTTACCGCC\n",

"...\n",

"ACGCATTTCTCCCGTGCAGGTCACATTTGCCCGACACGGCGGGGCAAGAGGCTTGAACAG\n",

"ACGTTCATTTTCCGTAAAACTGGCGTAATGTAAGCGTTTACCCACTATAGGTATTATCAT\n",

"GGCGACCATAAAAGATGTAGCCCGACTGGCCGGTGTTTCAGTCGCCACCGTTTCTCGCGT\n",

"TATTAACGATTCGCCAAAAGCCAGCGAAGCGTCCCGGCTGGCGGTAACCAGCGCAATGGA\n",

"GTCCCTGAGCTATCACCCTAACGCCAACGCGCGCGCGCTGGCACAGCAGGCAACGGAAAC\n",

"CCTCGGTCTGGTGGTCGGCGACGTTTCCGATCCTTTTTTCGGCGCGATGGTGAAAGCCGT\n",

"TGAACAGGTGGCGTATCACACCGGCAATTTTTTACTGATTGGCAACGGGTATCATAACGA\n",

"ACAAAAAGAGCGTCAGGCTATTGAACAGTTGATTCGTCATCGTTGCGCAGCGTTAGTGGT\n",

"GCACGCCAAAATGATTCCGGATGCGGACCTGGCCTCATTAATGAAGCAAATCCCCGGCAT\n",

"GGTGCTGATTAACCGCATTT\n"

\]

}

\],

"source": \[

"!head Yersinia_spi1_region_gc_map.fna\n",

"print('...')\n",

"!tail Yersinia_spi1_region_gc_map.fna"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"Looks good!"

\]

},

{

"cell_type": "markdown",

"metadata": {},

"source": \[

"## Computing environment"

\]

},

{

"cell_type": "code",

"execution_count": 15,

"metadata": {},

"outputs": \[

{

"name": "stdout",

"output_type": "stream",

"text": \[

"CPython 3.7.7\n",

"IPython 7.13.0\n",

"\n",

"numpy 1.18.1\n",

"jupyterlab 1.2.6\n"

\]

}

\],

"source": \[

"%load_ext watermark\n",

"%watermark -v -p numpy,jupyterlab"

\]

}

\],

"metadata": {

"anaconda-cloud": {},

"kernelspec": {

"display_name": "Python 3",

"language": "python",

"name": "python3"

},

"language_info": {

"codemirror_mode": {

"name": "ipython",

"version": 3

},

"file_extension": ".py",

"mimetype": "text/x-python",

"name": "python",

"nbconvert_exporter": "python",

"pygments_lexer": "ipython3",

"version": "3.7.7"

}

},

"nbformat": 4,

"nbformat_minor": 4

}

#The appearance of similar sequences in various Yersinia species points to a strong genetic relationship between the 
#phage holin family genes and this genus. In addition to confirming the genetic relationship between Yersinia species, 
#this finding may point towards an evolutionary path of these genes in the genus. The high representation of phage holin 
#family proteins in the hits lends support to the putative hypothesis that genes of interest identified from this new 
#Yersinia strain may well operate through a similar mechanism or bear close homology to such protein families.

# Further Data Analysis and Future Research Potential
#The BLAST results suggest the potential pathogenicity of NewYersinial, with discernible disparities in amino acid 
#sequence relative to W22703, potentially influencing its structural conformation and, consequently, functional 
#distinctions. In order to delineate the intricacies of structural and sequential variations across all target genes 
#between NewYersinia1 and W22703, SWISS model was employed for protein structure simulations. The resulting models, 
#exhibited in rainbow mode with blue denoting N-terminus and red indicating C-terminus, enabled a meticulous analysis.


# Holin

#Hypothesis:

#Sequence variations in NewYersinia1 holin, despite being located in transmembrane domains (alpha helices), 
#may not significantly alter its overall structure, potentially leading to functional divergence from W22703.

#Methods:

#SWISS-model was employed to generate simulated protein structures of NewYersinia1 and W22703 holins.
#Sequence alignment was performed to identify specific variations between the two proteins

#Results:

#Simulated models revealed that both holins possess similar overall structures, despite the presence of few sequence 
#variations within the transmembrane domains.

#Discussion:

#The observed sequence variations, while located in functionally important regions, might not significantly impact the 
#global protein fold of NewYersinial holin.
#This suggests that, while functional efficiency may differ between NewYersinia1 and W22703 holins, their core holin 
#function is likely preserved.

#Further research:

#Functional assays are necessary to definitively determine the impact of sequence variations on NewYersinial holin 
#activity and its role in potential pathogenicity.
#Detailed structural analyses, potentially through mutagenesis and experimental validation, could provide deeper 
#insights into the specific effects of these variations on protein-protein interactions and membrane dynamics.

#Key takeaways:

#Sequence variations in NewYersinial holin, despite their potential functional implications, do not appear to 
#drastically alter its overall structure compared to W22703.
#Further investigation is crucial to elucidate the precise functional consequences of these variations and their 
#contribution to NewYersinial's potential pathogenic behavior.

get_ipython().system('[Holin.png](attachment:a50c3496-2531-4a6e-b7dc-aa52df9a908a.png)')

get_ipython().system('[Holin] (https://github.com/2502495/Project-2/blob/5513b8cd6cd593ad4906945700b61b1d49b907ea/Holin.png)')
Figure 4. Holin SWISS Analysis

# Endolysin

#Hypothesis:

#Despite the presence of sequence variations, NewYersinia1 endolysin might exhibit a conserved structure and function 
#compared to W22703 due to the minimal nature of the variations and their location outside of critical functional domains.

#Methods:

#Sequence alignment was performed to identify specific variations between NewYersinia1 and W22703 endolysins.
#SWISS-model was employed to generate simulated protein structures for both endolysins.

#Results:

#Sequence alignment revealed few variations within the protein domains of NewYersinia1 endolysin.
#Simulated models of both endolysins displayed a high degree of structural similarity, with no significant differences 
#observed.

#Discussion:

#The minimal sequence variations identified, coupled with the preserved overall structure as depicted by the simulated 
#models, suggest that NewYersinia1 endolysin likely retains the core catalytic function of W22703.
#While functional variations in terms of efficiency or substrate specificity might exist, the fundamental endolysin 
#activity is likely conserved between the two strains.

#Further research:

#Biochemical assays are necessary to definitively assess the enzymatic activity and substrate preferences of NewYersinia1 
#endolysin compared to W22703.
#Mutagenesis studies could provide targeted insights into the specific roles of the identified sequence variations in 
#shaping the functional properties of the enzyme.

#Key takeaways:

#Despite sequence variations, NewYersinial endolysin exhibits a high degree of structural similarity to W22703, suggesting 
#a conserved core function.
#Further investigation is crucial to elucidate the precise functional implications of the sequence variations and their 
#potential impact on NewYersinial's biological processes.

#Visualization of Endolysin Structure:

#To further enhance the understanding of the structural similarities between the two endolysins, consider incorporating 
#the following visuals:

#Aligned protein sequences: Highlighting the specific amino acid variations between NewYersinia1 and W22703 endolysins.
#Superimposed simulated models: Illustrating the near-identical folds and domain arrangements of both endolysins.
#Close-up views of the variation sites: Providing a detailed look at the specific amino acid substitutions and their 
#potential impact on local protein interactions.
#By incorporating these visuals, you can effectively communicate the key findings of the structural analysis and their 
#implications for functional conservation between NewYersinia1 and W22703 endolysins

get_ipython().system('[Endolysin.png](attachment:8a7ad15d-6d06-412a-bc17-17243cb9e920.png)')

get_ipython().system(' [Endolysin] (https://github.com/2502495/Project-2/blob/5513b8cd6cd593ad4906945700b61b1d49b907ea/Endolysin.png)')

Figure 5. Endolysin SWISS Analysis

#i-Spanin

#Hypothesis:

#Despite the presence of sequence variations, particularly within the alpha-helix domain, NewYersinial i-spanin might 
#retain its core function due to a conserved C-terminus and minimal structural perturbations as predicted by simulated 
#structures. However, these variations could potentially lead to differences in conformational dynamics and pore properties.

#Analysis:

#Sequence alignment revealed variations primarily located within the alpha-helix domain of NewYersinial i-spanin.
#SWISS-model simulations indicated no significant differences in overall protein structure compared to W22703 i-spanin.

#Discussion:

#The conserved C-terminus and overall structural similarity suggest that NewYersinial i-spanin likely maintains the 
#ability to interact with o-spanin and initiate pore formation.
#However, the alpha-helix variations might influence the specific conformational changes and dynamics involved in pore 
#assembly. This could potentially lead to:
#Altered pore properties: Variations in amino acid side chains within the alpha-helices could impact pore size, selectivity, and stability.
#Modified conformational landscape: The variations might introduce new energy minima or alter the energy barriers 
#associated with different conformational states, potentially affecting the kinetics and efficiency of pore formation.

#Further Research:

#Mutagenesis studies: Targeted mutations within the alpha-helix domain could elucidate their specific roles in i-spanin 
#function and pore formation.
#Molecular dynamics simulations: Simulating the protein's behavior over time could provide deeper insights into the 
#conformational dynamics and pore characteristics of NewYersinial i-spanin compared to W22703.
#Functional assays: Assays investigating pore formation kinetics, stability, and substrate permeability could directly 
#assess the impact of sequence variations on i-spanin activity.

#Visualization of i-spanin Structure and Variations:

#To enhance understanding, consider incorporating the following visuals:

#Aligned protein sequences: Highlight the specific amino acid variations within the alpha-helix domain of NewYersinial 
#i-spanin.
#Superimposed simulated models: Illustrate the overall structural similarity between NewYersinial and W22703 i-spanins, 
#focusing on the conserved C-terminus and alpha-helix regions.
#Close-up views of the variation sites: Provide detailed visualizations of the specific amino acid substitutions within the 
#alpha-helices and their potential impact on local interactions and conformational flexibility.
#By incorporating these elements, you can effectively communicate the potential functional implications of the observed 
#sequence variations in NewYersinial i-spanin and pave the way for further investigation into its unique pore-forming 
#properties.

get_ipython().system('[i-spanin.png](attachment:885044b1-c697-4e71-8c8d-f12eb51fab77.png)')

get_ipython().system('[i-Spanin] (https://github.com/2502495/Project-2/blob/5513b8cd6cd593ad4906945700b61b1d49b907ea/i-spanin.png)')

Figure 6. i-Spanin SWISS Analysis

#o-Spanin

#Hypothesis:

#Despite the presence of some sequence variations, particularly within the alpha-helix domain, NewYersinial o-spanin might 
#retain its core function due to minimal structural perturbations as predicted by simulated models. However, these 
#variations could potentially influence subtle changes in inter-subunit interactions or pore properties.

#Analysis:

#Sequence alignment revealed few variations primarily located within the alpha-helix domain of NewYersinial o-spanin.
#SWISS-model simulations indicated no significant differences in overall protein structure compared to W22703 o-spanin.

#Discussion:

#The conserved overall structure and minimal deviation in alpha-helix conformation suggest that NewYersinial o-spanin 
#likely maintains its ability to interact with i-spanin and participate in pore formation.

#However, the alpha-helix variations, although seemingly minor, could potentially lead to:

#Altered inter-subunit interactions: Subtle changes in amino acid side chains within the alpha-helices could impact the 
#strength and specificity of interactions between NewYersinial and W22703 o-spanins, potentially affecting the stability 
#and efficiency of pore assembly.
#Modified pore characteristics: Even minute changes in alpha-helix geometry or surface properties could influence the size,
#selectivity, and gating mechanisms of the assembled pore, potentially impacting its substrate passage efficiency.

#Further Research:

#Mutagenesis studies: Targeted mutations within the alpha-helix domain could elucidate their specific roles in o-spanin 
#interactions and pore function.
#Molecular docking simulations: Simulating the interaction between NewYersinial and W22703 o-spanins could provide insights 
#into how the identified variations might affect inter-subunit binding interfaces and pore stability.
#Functional assays: Assays investigating pore formation kinetics, stability, and substrate permeability could directly 
#assess the impact of sequence variations on o-spanin activity.

#Visualization of o-spanin Structure and Variations:

#To enhance understanding, consider incorporating the following visuals:

#Aligned protein sequences: Highlight the specific amino acid variations within the alpha-helix domain of NewYersinial 
#o-spanin.
#Superimposed simulated models: Illustrate the overall structural similarity between NewYersinial and W22703 o-spanins, 
#focusing on the conserved alpha-helix regions.
#Close-up views of the variation sites: Provide detailed visualizations of the specific amino acid substitutions within the 
#alpha-helices and their potential impact on local interactions and surface properties.
#By incorporating these elements, you can effectively communicate the potential functional implications of the observed 
#sequence variations in NewYersinial o-spanin and pave the way for further investigation into its unique pore-forming 
#properties.

get_ipython().system('[o-spanin.png](attachment:68d3b06a-1582-4827-9154-fe162e3446f0.png)')

get_ipython().system('[o-Spanin] (https://github.com/2502495/Project-2/blob/5513b8cd6cd593ad4906945700b61b1d49b907ea/o-spanin.png)')

Figure 7. o-Spanin SWISS Analysis


# #7.0 APPLICATIONS AND INSIGHTS
# 
# The BLAST analysis gives an in-depth look into the genetic architecture of the new Yersinia strain, especially concerning its phage holin-linked genes. Since the conservation and prevalence of these genes imply an important aspect in the bacterium’s life cycle, pathogenic mechanisms may be influenced by them. Learning how these conserved genes function and interact can provide insights into Yersinia’s techniques for survival, virulence, and adaptation

# #8.0 REAL WORLD RELEVANCE
# 
# The BLAST analysis gives an in-depth look into the genetic architecture of the new Yersinia strain, especially concerning its phage holin-linked genes. Since the conservation and prevalence of these genes imply an important aspect in the bacterium’s life cycle, pathogenic mechanisms may be influenced by them. Learning how these conserved genes function and interact can provide insights into Yersinia’s techniques for survival, virulence, and adaptation

# #9.0 FUTURE RESEARCH  
# 
# The genetic diversity shown in the genomic comparison with different Yersinia species emphasizes the complex evolutionary processes within the genus. It allows for further research opportunities regarding the role these conserved genes play in enhancing the adaptability and pathogenicity of Yersinia, thus understanding bacterial evolution, and host-pathogen interactions better.

# #10.0 CONCLUSION
# 
# Overall, the study using BLAST analysis of the novel Yersinia strain detected phage holin-related genes in this bacterium; these latter possess a high degree of conservation across all Yersinia species and propose possible functional roles. These observations contribute greatly to the knowledge of the genetic heterogeneity and virulence features of Yersinia, laying a strong basis for future genomic and functional studies in microbiology

# # REFERENCES/CITATIONS
# 
# #1. Original Publication
# Springer K, Reuter S, Knüpfer M, Schmauder L, Sänger PA, Felsl A, Fuchs TM. Activity of a Holin-Endolysin System in the Insecticidal Pathogenicity Island of Yersinia enterocolitica. J Bacteriol. 2018 Jul 25;200(16):e00180-18. doi: 10.1128/JB.00180-18. PMID: 29866807; PMCID: PMC6060350.
# 
# #2. Documentation
# 1. Day, R.A. (1998), How to Write and Publish a Scientific Paper, Westport: Oryx Press.
# 2. O'Connor, M. (1992), Writing Successfully in Science, London: Chapman & Hall.
# 3. Montgomery, S.L. (2003), The Chicago Guide to Communicating Science, London: University of Chicago Press.
# 4. Porush, D. (1995), A Short Guide to Writing About Science, New York: HarperCollins College Publishers.
# 5. Malmfors, B., Garnsworthy, P., and Grossman, M. (2000), Writing and Presenting Scientific Papers, Nottingham: Nottingham University Press.
# 6. van Emden, J., and Easteal, J. (1996), Technical Writing and Speaking: An Introduction, Berkshire: McGraw-Hill Publishing Company.
# 7. Barrass, R. (1995), Scientists Must Write: A Guide to Better Writing for Scientists, Engineers and Students, London: Chapman & Hall.
# 8. Adams, R.L., Adams, I.P., Lindow, S.W., Zhong, W., and Atkin, S.L. (2005), Somatostatin receptors 2 and 5 are preferentially expressed in proliferating endothelium, British Journal Of Cancer 92:1493-1498.
# 9. Davis, M. (1997), Scientific Papers and Presentations, San Diego: Academic Press.
# 10. Jones, A., Reed, R., and Weyers, J. (1994), Practical Skills in Biology, Harlow: Longman Scientific & Technical.
# 11. American Society for Microbiology (2006), Instructions to Authors, Journal of Bacteriology 188: 1-18.
