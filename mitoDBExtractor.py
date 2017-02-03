# mitoDBExtractor.py
# Reads in each line of a file containing id,family,genus,species,subspecies with no spaces.
# For a each organism, finds the closest relative in the database
# Check if genus exists in mitoDB.sqlite, Do same for Subfamily, etc.
# Does not look in NCBI database for closest relative, will make another program to do this.

################ Note ############################
# consider possibly saving the full taxonomy for an organism
# so we don't have to look it up twice
# and maybe save the closest reference in a new table as well?
# or will this be a lot of extra info we don't want in the end


from collections import OrderedDict
from bs4 import BeautifulSoup
from ast import literal_eval
import sqlite3
import httplib
from mitoDBmaker import get_response
import sys

dbConnection = sqlite3.connect('./mitoDB.sqlite')
dbConnection.row_factory = sqlite3.Row


def checkSynonyms(genus, species):
    trans_genus = None
    trans_species = None
    cur = dbConnection.cursor()
    organism = str(genus + " " + species)

    select = 'SELECT * FROM Synonyms WHERE Synonym=?'
    cur.execute(select, (organism,))
    fetched = cur.fetchone()

    if fetched is not None:  # not a synonym
        genus, species = str(fetched["Synonym"]).split()
        trans_genus, trans_species = str(fetched["ReferenceName"]).split()
        print str(genus) + " translated to " + str(trans_genus) + " and " + str(species) + " translated to " + str(trans_species)

        # trans_genus=trans_genus.replace("u'", "")

    return trans_genus, trans_species


def get_taxonomy(genus, species):  # different return than DBMaker, uses dict; returns 3 fewer results

    # AllOtherRank=[]
    Superkingdom = None
    Kingdom = None
    Superphylum = None
    Phylum = None
    Subphylum = None
    Class = None
    Superorder = None
    Order = None
    Suborder = None
    Infraorder = None
    Parvorder = None
    Superfamily = None
    Family = None
    Subfamily = None

    #####################################################################
    # get id from taxonomy database for specified genus and species
    #####################################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=Eumetopias+jubatus
    url = '/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + str(genus) + '+' + str(species)
    response = get_response(url)
    rData = response.read()
    id_soup = BeautifulSoup(rData, 'xml')
	 
    if id_soup is None:
        print "No taxonomy xml returned"
        return


    if id_soup.eSearchResult is None or id_soup.eSearchResult.IdList.Id is None:
        print "No id in taxonomy xml for: " + str(genus) + " " + str(species)
        return

    id = id_soup.eSearchResult.IdList.Id.string

    #####################################################
    # Get taxonomy xml for the id
    #####################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=34886&retmode=xml
    url = '/entrez/eutils/efetch.fcgi?db=taxonomy&id=' + id + '&retmode=xml'
    response = get_response(url)
    rData = response.read()
    fastaSoup = BeautifulSoup(rData, 'xml')

    ###################################################
    # check that the genus and species names are correct
    ###################################################

    orgList = fastaSoup.find("ScientificName")  # start value in list format
    organism = orgList.contents[0]  # Genus species
    org_name = str(organism).split()

    Taxon = fastaSoup.find_all("Taxon")

    if org_name[0] != genus:
        correct_organism = False
        Synonyms = Taxon[0].find_all("Synonym")
        Synonyms += Taxon[0].find_all("GenbankSynonym")

        for synonym in Synonyms:
            if synonym is None or synonym.string is None:
                continue
            if len(org_name) < 2:  # missing species name, genus only
                # check for synonym
                if genus in synonym.string:
                    correct_organism = True
                    break
            else:
                scientific_name = genus + " " + species
                if scientific_name in synonym.string:
                    correct_organism = True
                    break
        if correct_organism == False:
            print "synonym not found for " + str(genus) + " " + str(species)
            return

    ###############################################
    # Get the lineage
    ###############################################

    for taxa in Taxon:

        rank = taxa.find("Rank")
        # if rank.string == "no rank":
        # AllOtherRank += [str(taxa.ScientificName.string)]
        if rank.string == "superkingdom":
            Superkingdom = taxa.ScientificName.string
        elif rank.string == "kingdom":
            Kingdom = taxa.ScientificName.string
        elif rank.string == "superphylum":
            Superphylum = taxa.ScientificName.string
        elif rank.string == "phylum":
            Phylum = taxa.ScientificName.string
        elif rank.string == "subphylum":
            Subphylum = taxa.ScientificName.string
        elif rank.string == "class":
            Class = taxa.ScientificName.string
        elif rank.string == "superorder":
            Superorder = taxa.ScientificName.string
        elif rank.string == "order":
            Order = taxa.ScientificName.string
        elif rank.string == "suborder":
            Suborder = taxa.ScientificName.string
        elif rank.string == "infraorder":
            Infraorder = taxa.ScientificName.string
        elif rank.string == "parvorder":
            Parvorder = taxa.ScientificName.string
        elif rank.string == "superfamily":
            Superfamily = taxa.ScientificName.string
        elif rank.string == "family":
            Family = taxa.ScientificName.string
        elif rank.string == "subfamily":
            Subfamily = taxa.ScientificName.string
        elif rank.string == "genus":
            continue
        elif rank.string == "species":
            continue
        # Subspecies is taken from the organism name in the xml, usually None is given
        # else:
        elif rank.string != "no rank":
            print "unknown rank " + str(rank.string) + " " + str(taxa.ScientificName.string)
            # AllOtherRank += [str(taxa.ScientificName.string)]

            # if rank.string=="no rank" and taxa.ScientificName.string=="Cetartiodactyla":
            #	order = "Cetartiodactyla"

    return OrderedDict([
        ("Subfamily", Subfamily),
        ("Family", Family),
        ("Superfamily", Superfamily),
        ("Parvorder", Parvorder),
        ("Infraorder", Infraorder),
        ("Order", Order),
        ("Suborder", Suborder),
        # ("AllOtherRank", str(AllOtherRank)),
        ("Superorder", Superorder),
        ("Class", Class),
        ("Subphylum", Subphylum),
        ("Phylum", Phylum),
        ("SuperPhylum", Superphylum),
        ("Kingdom", Kingdom),
        ("Superkingdom", Superkingdom)])


def checkDb(rank, rankValue, geneList):
    missing_gene_list = []
    cur = dbConnection.cursor()
    fasta = ""

    for gene in geneList:
        # check if the gene is present at the given taxonomic level
        if rank == "Species":  # look for genus and species, since two Genus can potentially have same species
            select = 'SELECT * FROM RefGenes WHERE Genus=? AND Species=? AND GeneName=?'
            cur.execute(select, (genus, species, gene))
        else:
            # select = 'SELECT * FROM RefGenes WHERE Genus=? AND GeneName=?'
            select = 'SELECT * FROM RefGenes WHERE "' + rank + '"=? AND GeneName=?'
            cur.execute(select, (rankValue, gene))
        results = cur.fetchone()

        if results is not None:  # returns one row, assumes cellCount=24 but this will change if attributes change
            #print "fetching results for " + str(results["Genus"]) + " " + str(results["Species"]) + "for gene " + gene

            sequence = str(results["GeneSequence"])
            acc_num = str(results["AccessionNumber"])
            genus_used = str(results["Genus"])
            species_used = str(results["Species"])
            partial_flag = str(results["PartialFlag"])
            if partial_flag == "0":
                partial_or_complete = "complete"
            if partial_flag == "1":
                partial_or_complete = "partial"
            #print "fasta-ing"
            header = ">" + gene + '_' + genus_used + '_' + species_used + '_' + acc_num + '_' + rank + '_' + \
                     rankValue + '_' + partial_or_complete
            fasta += header + "\n" + sequence + "\n"

        else:  # e.g., no gene present for organism of the given rank
            missing_gene_list += [gene]  # still missing, keep in list
            #print "making the missing"
    #print "geneList is " + str(geneList) + "and missing_gene_list is " + str(missing_gene_list)
    return missing_gene_list, fasta


if __name__ == '__main__':

    try:
        sample_list = sys.argv[1]  # raw data list from Alan's lab [id, family, genus, species]

    except IndexError:
        # print "please specify a filename"
        print "using default filename"
    	#sample_list = "Zaher18Samples.csv"
    	sample_list = "SampleList.csv"

    try:
        geneList = literal_eval(sys.argv[2])  # i.e. geneList=[\'COX1\', \'ND2\',\'12S\',\'16S\',\'ND5\']
        #populateRelatives = sys.argv[3]
        # Check that genes are qualified geneList names, to avoid redundant DB entries
        for gene in geneList:
            if gene.upper() == "COX1":
                continue
            if gene.upper() == "ND2":
                continue
            if gene.upper() == "12S":
                continue
            if gene.upper() == "16S":
                continue
            if gene.upper() == "ND5":
                continue
            if gene.upper() == "CYTB":
                continue
            # if gene.upper()=="";
            #	continue
            else:
                print "Please use one of the following gene names: ['COX1', 'ND2','12S','16S','ND5']"

    except IndexError:
        print "Using default gene list and no relative population in database"
        geneList = ['COX1', 'CYTB', 'ND2', '12S', '16S', 'ND5']
        #populateRelatives = "No"




    with open(sample_list, 'r') as sample:
        lines = sample.readlines()

        for line in lines:
            # read in list of Genus species names
            # Get the genus and species name from comma delimited file (e.x. ZaherRaw.csv)
            # line = line.strip() #should not have any spaces
            line = line.replace('\n', '')
            cells = line.split(',')
            id = cells[0]
            family = cells[1]
            genus = cells[2]
            species = cells[3]
            trans_genus = None
            trans_species = None
            fasta = ""
            missing_gene_list = list(geneList)  # LISTS ARE MUTABLE


            ##### Look for same subspecies in sqlite first????###
            if species != "sp.":
                rank0 = "Species"  # the rank as it appears in the Database
                rankValue0 = species
                missing_gene_list, fasta_lines = checkDb(rank0, rankValue0, missing_gene_list)
                if fasta_lines is not None:
                    fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)
                    # if populateRelatives=="yes" AND fasta_lines is None:
                    # giIDs, trans_genus, trans_species = mitoDBmaker.get_ids(genus, species, subspecies, geneOrGenome)
                else:
                    # Check synonyms in DB, if exists then use reference name
                    trans_genus, trans_species = checkSynonyms(genus, species)
                    if trans_genus is not None:
                        print " translated " + str(genus) + " to " + str(trans_genus)
                        genus = trans_genus
                        print " translated " + str(species) + " to" + str(trans_species)
                        species = trans_species
                        missing_gene_list, fasta_lines = checkDb(rank0, rankValue0, missing_gene_list)
                        if fasta_lines is not None:
                            fasta += fasta_lines

            if len(missing_gene_list) != 0:
                rank1 = "Genus"  # the rank as it appears in the Database
                print "Checking genus for " + str(genus) + " " + str(species)
                rankValue1 = genus
                missing_gene_list, fasta_lines = checkDb(rank1, rankValue1, missing_gene_list)
                if fasta_lines is not None:
                    fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)
                    # else:
                    #	giIDs, transRankValue =  mitoDBRelativeMaker.get_ids(rankValue, gene)

            if len(missing_gene_list) != 0:  # No genus present, search for full taxonomy
                # If genes remain to be found, find taxonomy and use relative for those genes
                # AllOtherRank,Superkingdom,Kingdom,Superphylum,Phylum,Subphylum,Class,Superorder,Order,Suborder,Infraorder,Parvorder,Superfamily,Family,Subfamily = get_full_taxonomy(genus,species)
                # ranks = [Subfamily, Family, Superfamily, Parvorder, Infraorder, Order, AllOtherRank, Superorder, Class, Subphylum, Phylum, SuperPhylum, Kingdom, SuperKingdom]

                temp_taxonomy = get_taxonomy(genus, species)
                if temp_taxonomy is None:
                    print "####### No taxonomy for " + str(genus) + " " + str(species)
                    temp_taxonomy = get_taxonomy(genus, "")
                    if temp_taxonomy is None:
                        print "## No taxonomy for " + str(genus)
                        if temp_taxonomy is None:
                            temp_taxonomy = get_taxonomy(family, "")
                            if temp_taxonomy is None:
                                print "No taxonomy for " + str(family) + " " + str(genus) + " " + str(species)
                                temp_taxonomy = taxonomy = OrderedDict([("Family", family)])
                                print "update seed file manually above family level"
                            else:
                                print "family found"
                temp_taxonomy = taxonomy = OrderedDict([("Family", family)])

                #print temp_taxonomy
                for rank, rankValue in temp_taxonomy.items():

                    if rankValue is None:
                        continue
                    else:
                        print "iteration: " + str(rank) + " " + str(rankValue) + " " + str(missing_gene_list[0]) + " for " + str(genus) + " " + str(species)
                        missing_gene_list, fasta_lines = checkDb(rank, rankValue, missing_gene_list)
                        if fasta_lines is not None:
                            fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)

                        if len(missing_gene_list) == 0:
                            # save (append) gene in fasta file for that species
                            # line of info followed by line of sequence
                            # filePrefix = genus + '_' + species #Want to use id instead
                            # f_out= open(filePrefix + ".seeds", 'w')
                            f_out = open("./seeds/" + id + ".seeds", 'w')
                            f_out.write(fasta)
                            f_out.close()
                            break

            else:  # missing_gene_list is None, species or Genus was present
                #print "writing to file"
                # filePrefix = genus + '_' + species
                f_out = open("./seeds/" + id + ".seeds", 'w')
                f_out.write(fasta)
                f_out.close()
                #### NOTE: Still need to add an optional check to NCBI when not found in database ####
                ##### Will write over fasta files of same species with a different subspecies
                #### Need to handle AllOtherRank for optional NCBI search

        dbConnection.close()


        # Create an option to search for relatives in NCBI and populate database if none found
        # modify superextender to use only the line of the gene you want (sys args)
