#!/usr/bin/env python
'''
This code works in conjunction with other python scripts to determine if lab-developed PCR assays can 
theoretically detect all circulating strains of a target organism (i.e. findAmpliconInSequence.py, 
extractAmpliconFromEmbossOutput.py, ampliconSequevarsToDict.py) and is meant to be used as part of a Galaxy
workflow.
Input: (1) fasta file of amplicon sequevars (and an optional 'safe list' of previously vetted sequences)
	   (2) fasta file of primers and probes
	   (3) tab-delim text results of blastn-short searches of amplicon sequevars against a primer/probe fasta
	   (4) output file name
Output:(1) fasta file of sequevars of extracted and concatenated hybridization regions, containing ONLY 
		primer/probe hybridization sites in the order F primer - probe - R primer.
		(2) search results summary text file detailing how many sequevars were found, what they were, and 
		which amplicon sequences fell into the various sequevar groupings.

usage without safelist: python extract_hyb_regions_from_amplicon_sequevars.py amplicon_sequevars.fasta
	oligos.fasta blast_hit_results.tabular hyb_site_sequevars
usage with safelist: python extract_hyb_regions_from_amplicon_sequevars.py -s hybsite_safelist.fasta
	amplicon_sequevars.fasta oligos.fasta blast_hit_results.tabular hyb_site_sequevars
Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory, June 2018
'''

import sys,string,os, time, Bio, re, argparse
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

#parse cmd line args - user can specify a list of vetted primer/probe hybridization sites sequences
parser = argparse.ArgumentParser()
parser.add_argument("-s", dest = "hyb_sites_safelist", 
	help = "specify filename of previously vetted hyb-site sequences",action = "store")
parser.add_argument("fastaToParse")
parser.add_argument("oligo_fasta")
parser.add_argument("blast_hits")
parser.add_argument("outFileHandle")
arguments = parser.parse_args()
#print("ARGUMENTS:")
#print(arguments)
#print('safelist = ', arguments.hyb_sites_safelist)

#construct output filenames from user-specified output file name
outputFastaHandle = arguments.outFileHandle + "_hyb_site_sequevars.fasta"
results_handle = arguments.outFileHandle + "_hyb_sites_summary.txt" #output file of GB access #'s of sequences requiring further investigation
results_summary = open(results_handle,'w')
all_extracted_hyb_regions_handle = arguments.outFileHandle + "_all_extracted_hyb_regions.fasta"

def output_information_to_files(hyb_site_sequevars):
	'''Output search summary text file detailing files used in analysis, how many sequevars were found, and
	sequence ID's associated with these sequevars.'''
	localtime = time.asctime(time.localtime(time.time())) #date and time of analysis
	results_summary.write("---------------------------------------------------------------------------\n")
	results_summary.write("RESULTS OF SEARCH FOR EXACT MATCHES TO OLIGONUCLEOTIDES IN QUERY SEQUENCES:\n")
	results_summary.write("---------------------------------------------------------------------------\n\n")
	#information about files used in analysis
	results_summary.write("Analysis as of: " + localtime + "\n\n")
	results_summary.write("Input Fasta: %s\n" % (arguments.fastaToParse))
	results_summary.write("Oligonucleotides searched for: %s\n" % (arguments.oligo_fasta))
	results_summary.write("Blast results file: %s\n" % (arguments.blast_hits))
	if arguments.hyb_sites_safelist: #if safelist used, record what file was used
		safelist_handle = str(arguments.hyb_sites_safelist)
		results_summary.write("Safelist file: %s\n" % (arguments.hyb_sites_safelist))
	else:
		results_summary.write("Safelist file: N/A\n")
	results_summary.write("---------------------------------------------------------------------------\n")
	#print each unique concatenated hyb site sequence, followed by GB # of records with that sequence
	if len(hyb_site_sequevars) > 0: #if there are sequevars to investigate, print their id's
		output_sequevars_to_fasta(hyb_site_sequevars) #output fasta file of sequevars i.e.representative record
		results_summary.write("%i Unique Sequevar(s) of Concatenated Oligo Hybridization Sites Found: \n" % (len(hyb_site_sequevars)))
		for sequevar in hyb_site_sequevars:
			recordList = hyb_site_sequevars[sequevar] #access the list of SeqRecords with the sequevar
			#write sequevar to file, along with id's of sequences with the sequevar
			results_summary.write("\n\n%i Record(s) with sequevar: %s" % (len(recordList),sequevar))
			#print number of records have the sequevar to console
			print("--> %i Record(s) with sequevar: %s" % ((len(recordList)), sequevar))
			for record in recordList:
				results_summary.write("\n\t" + record.id) #write records with this sequence to file
	else: #if no sequences require further analysis
		print("No sequences requiring further investigation")
		results_summary.write("No sequences requiring further investigation")
	print("\nEND OF RESULTS\n")
	results_summary.write("\n---------------------------------------------------------------------------")
	results_summary.write("\n\nEND OF RESULTS\n")
	#print locations of all output files for user's reference
	results_summary.write("\nResults Summary (this file): " + results_handle) #write output filenames
	print("Search Results file: " + results_handle) #print names of info and fasta output files
	#print output filepaths to console and to search results file
	if len(hyb_site_sequevars) > 0:
		results_summary.write("\nSequevars of concatenated oligo hybridization regions: " + outputFastaHandle + "\n")
		print("Sequevars of concatenated oligo hybridization regions: " + outputFastaHandle + "\n")
	else:
		results_summary.write("\nSequevars of concatenated oligo hybridization regions: N/A\n")
		print("Sequevars of concatenated oligo hybridization regions: N/A\n")
	results_summary.close() #close the search results file
	return

def extract_hyb_sites_using_blast_results(ampliconSequevar, hitlines, oligoList):
	'''Check hitlines for query id matching that of ampliconSequevar. Use the start and end points of the
	subject alignment to extract corresponding parts of the ampliconSequevar sequence and concatenate them'''
	fprimerSite = "" #regions of ampliconSequevar sequence corresponding to the respective oligos
	rprimerSite = ""
	probeSite = ""
	delimiter = "" #placeholder to visually separate oligo hydridization sites
	reverse = False #by default, the oligo orientation on the template is forward
	for hit in hitlines: #assign blast hit line parameters
		hit_parameters = hit.rstrip().split('\t')
		templateSeqID = hit_parameters[0] #query sequence id (i.e which ampliconSequevar was blasted)
		oligoSeqID = hit_parameters[1] #subject sequence id (i.e. which oligo)
		align_length = int(hit_parameters[2]) #sequence match between ampliconSequevar and oligo
		templateStart = int(hit_parameters[3]) #alignment start on ampliconSequevar
		templateEnd = int(hit_parameters[4]) #alignment end on ampliconSequevar
		oligoStart = int(hit_parameters[5]) #alignment start on oligo i.e. subj start in blast
		oligoEnd = int(hit_parameters[6]) #alignment end on oligo i.e. subj end in blast hit
		oligoForward = oligoStart < oligoEnd #determine directionality of oligo to template

		if templateSeqID == ampliconSequevar.id: #if blast hit line is for this ampliconSequevar
			oligo_length = check_oligo_length(oligoSeqID, oligoList) #find oligo length in this hit line
			tStart = templateStart
			tEnd = templateEnd
			if align_length == oligo_length: # if entire oligo hybridizes, slice region from template
				hybSlice = extract_hyb_region(ampliconSequevar, templateStart, templateEnd)
			else:
				missing = oligo_length - align_length #determine how much of oligo isn't aligning
				#check for oligo directionality first!
				if oligoForward: #forward oligo
					if (oligoEnd < oligo_length) and (oligoStart > 1): #oligo not hybridizing at both ends
						#move start and end points of template
						tEnd = templateEnd + (oligo_length - oligoEnd)
						tStart = templateStart - (oligoStart -1)
					elif (oligoEnd < oligo_length): #oligo not hybridizing at 3' end
						#move template end point right to where the oligo end point SHOULD be
						tEnd = templateEnd + (oligo_length - oligoEnd)
						tStart = templateStart
					elif(oligoStart > 1): #oligo not hybridizing at 5' end
						tStart = templateStart - (oligoStart - 1)
						tEnd = templateEnd
				else: #reverse oligo
					if (oligoStart < oligo_length) and (oligoEnd > 1): #not hybridizing at both ends
						#move start and end points of template
						tStart = templateStart - (oligo_length - oligoStart)
						tEnd = templateEnd + (oligoEnd - 1)
					elif (oligoStart < oligo_length): #oligo not hybridizing at 3' end
						#move template start point right to where the oligo end point SHOULD be
						tStart = templateStart - (oligo_length - oligoStart)
						tEnd = templateEnd
					elif (oligoEnd > 1): # 3' end of reverse primer not hybridizing
						#move template start point left to where oligo should hybridize
						tStart = templateStart
						tEnd = templateEnd + (oligoEnd - 1)
					else: #reverse oligo
						print("ERROR: I truly hadn't account for this happening...")
						#move template endpoint forward to where oligo start point SHOULD be
						tEnd = templateEnd + (oligoStart - 1)
						tStart = templateStart
				#extract region from ampliconSequevar
				hybSlice = extract_hyb_region(ampliconSequevar,tStart,tEnd)

			#determine which oligo the hyb region sliceit belongs to
			if oligoSeqID.endswith("_F"):
				if fprimerSite == "": #to avoid the first (i.e. best) blast hit being over-written by a subsequent blast hit,
					fprimerSite = hybSlice # insist that the variable doesn't already contain a sequence!
				else:
					print("\tMultiple blast hits of %s to %s oligo --> only first hit accepted!" % (ampliconSequevar.id,oligoSeqID))
			elif oligoSeqID.endswith("_P"):
				if probeSite == "":
					probeSite = hybSlice
				else:
					print("\tMultiple blast hits of %s to %s oligo --> only first hit accepted!" % (ampliconSequevar.id,oligoSeqID))
			elif oligoSeqID.endswith("_R"):
				if rprimerSite == "":
					rprimerSite = hybSlice
				else:
					print("\tMultiple blast hits of %s to %s oligo --> only first hit accepted!" % (ampliconSequevar.id,oligoSeqID))
			else:
				print("ERROR: SEQUEVAR MATCH TO INVALID PROBE!")
	#concatenate hyb slices in correct order
	return fprimerSite + delimiter + probeSite + delimiter + rprimerSite

def check_oligo_length(name,oligoList):
	'''Search oligoList for an oligo corresponding to 'name' and return sequence length.'''
	for oligo in oligoList:
		if oligo.id == name:
			return len(oligo.seq)

def extract_hyb_region(ampliconSequevar, tStart, tEnd):
	'''Return a subsequence from ampliconSequevar usingthe start and end indices.'''
	#since we use blastn query start/end points to slice from ampliconSequevar, query start is ALWAYS less than query end.
	return ampliconSequevar.seq[tStart - 1:tEnd]
	
def output_sequevars_to_fasta(hyb_site_sequevars):
    '''Produce an output fasta file of one representative sequence for each different sequevar found.'''
    output_fasta = open(outputFastaHandle, 'w')
    hyb_site_sequevars_list = []
    for sequevar in hyb_site_sequevars:
        first_record = hyb_site_sequevars[sequevar][0] #grab the first SeqRecord in the list sharing the sequevar
        hyb_site_sequevars_list.append(first_record)  #add it to the list of SeqRecords
    SeqIO.write(hyb_site_sequevars_list, output_fasta, "fasta")
    output_fasta.close() #close the fasta file when done writing to it
    return
  
print("\nEXTRACTING OLIGONUCLEOTIDE HYBRIDIZATION SITES...\n")
#Step 1 - read blast hit results
with open(arguments.blast_hits, 'r') as blastFile:
	hitlines = blastFile.readlines() #the list of tab-separated lines
	
#Step 2 - read oligos from fasta file into SeqRecords, uppercase them, and add to a list
with open(arguments.oligo_fasta, 'r') as oligoFile:
	oligoList = [oligo.upper() for oligo in list(SeqIO.parse(oligoFile,"fasta",alphabet=IUPAC.ambiguous_dna))]
	for rec in oligoList:
		print("%s:  Length: %i  Sequence(5' -> 3'): %s" % (rec.id, len(rec.seq),rec.seq))
	print()

#Step 3 - read unique amplicon sequevars from fasta file into SeqRecords, uppercase, and add to a list
with open(arguments.fastaToParse,'r') as inFile:
	ampliconSequevarList = [rec.upper() for rec in list(SeqIO.parse(inFile,"fasta",alphabet=IUPAC.ambiguous_dna))]
	hyb_site_sequevars = {} # empty dictionary of ampliconSequevars comprised only of primer/probe hyb sites
	safe_hyb_sites_list = []  #empty list of prevetted 'safe' hyb site ampliconSequevars
	allHybSitesList = [] #list to store concatenated extracted hybridization sites

	if arguments.hyb_sites_safelist: #if safelist specified by user, read these sequences into a list
		safe_hyb_sites_list = [record.seq for record in SeqIO.parse(arguments.hyb_sites_safelist,"fasta",
			alphabet = IUPAC.ambiguous_dna)]
		print("Comparing extracted hybridization regions to %i previously-tested sequevar(s)..." % (len(safe_hyb_sites_list)))
		'''for seq in safe_hyb_sites_list:
			print(seq) #print all the vetted unique hybridization site sequences'''

	for ampliconSequevar in ampliconSequevarList: #extract primer/probe binding sites from each unique ampliconSequevar
		concatenated_hyb_site = extract_hyb_sites_using_blast_results(ampliconSequevar,hitlines, oligoList) #concatenated, extracted primer/probe hyb sites
		#print("%s HybSites (%i): %s" % (ampliconSequevar.id,len(concatenated_hyb_site),concatenated_hyb_site))
		#create SeqRecords of concatenated_hyb_site and add to a list of SeqRecords
		hybSiteSeqRec = SeqRecord(concatenated_hyb_site, ampliconSequevar.id, description ="concat_oligo_hyb_sites")
		allHybSitesList.append(hybSiteSeqRec)

#Step 4 - create dict of concatenated hyb region sequevars, excluding sequevars present in safelist
for record in allHybSitesList:
	sequence = str(record.seq)
	#check to see if this sequence is present in the safelist
	if (sequence in safe_hyb_sites_list):
		print("--> Sequence from %s previously tested in assay." % (record.id))
	else: #add sequence as new dict key, mapping to a list of SeqRecords with record as the first list item
		if sequence in hyb_site_sequevars: #if the sequence is already a key in the dict
			hyb_site_sequevars[sequence].append(record) #add corresponding SeqRecord to the list of SeqRecords this key maps to
		else: #this sequence is unique to the dictionary so add new dict element consisting of sequence: <list> SeqRecords
			hyb_site_sequevars[sequence] = [record] # this record will element zero in the theSeqRecord list

#get sorted list of sequevar keys
sorted_unique_sequence_keys = sorted(hyb_site_sequevars.keys())
#process each list of SeqRecords represented by a sequevar and write to an information file for follow-up
output_information_to_files(hyb_site_sequevars) #output information .txt file



		


