# mitoDBRelativeMaker.py
# reads [SampleList]NeedReference.csv
# for each item in list looks for closest relative in the NCBI database
import sys
import sqlite3
import time
import string
from collections import OrderedDict
from bs4 import BeautifulSoup
from mitoDBmaker import get_response
from mitoDBmaker import get_full_taxonomy
from mitoDBExtractor import get_taxonomy
from mitoDBExtractor import checkSynonyms
from difflib import SequenceMatcher

lastRequestTime = 0


def get_ids(rank_value, gene):
    global lastRequestTime

    #############################################################
    # get the GenBank accession number via esearch, ex 256557273
    #############################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Atractaspis+bibronii+COX1
    # Note: check that if subspecies isn't found, still finds species

    url = '/entrez/eutils/esearch.fcgi?db=nuccore&term=' + rank_value + '+' + gene

    time_since_last_request = time.time() - lastRequestTime
    if time_since_last_request < 0.3:
        time.sleep(0.3 - time_since_last_request)
    response = get_response(url)
    lastRequestTime = time.time()

    rData = response.read()
    giSoup = BeautifulSoup(rData, 'xml')

    giIdList = giSoup.IdList
    translationSet = giSoup.TranslationSet.Translation

    if translationSet is None:  # Genus not found, same geneOrGenome returned
        print "            No query translation for " + str(rank_value) + " " + str(gene)
        return None, None, None

    translationString = translationSet.To.string
    translation = translationString.split('"')
    translation = translation[1].split(' ')
    # How does this translation work for genus, family, etc?
    if len(translation) == 1:
        transRankValue = translation[0]

    else:
        print "        Improper query translation " + str(translation) + "for " + str(rank_value) + " " + str(gene)
        return None, None, None

    if giIdList.Id is None:  # No id returned for gene but name translation available
        return None, transRankValue

    # giIDs = giSoup.IdList.Id.string #only returns first id
    giIdList = giSoup.IdList
    giIDs = giIdList.find_all("Id")  # return a list of IDs
    giIDs = [giId.string for giId in giIDs]
    return giIDs, transRankValue  # returns "translated" rank Value, for subfamily level and higher


def checkDb(rank, rankValue, gene):
    cur = dbConnection.cursor()
    missing_gene = False

    # select = 'SELECT * FROM RefGenes WHERE Genus=? AND GeneName=?'
    select = 'SELECT * FROM RefGenes WHERE "' + rank + '"=? AND GeneName=?'

    cur.execute(select, (rankValue, gene))

    results = cur.fetchone()

    if results is not None:  # returns one row, assumes cellCount=24 but this will change if attributes change
        # may want to check partial flag and only accept complete references at some point?
        # partial_flag = str(results["PartialFlag"])
        # acc_num = str(results["AccessionNumber"])
        # if partial_flag == "0"
        # boolean asks if partial. # zero is false: thus complete, 1 is true: thus partial
        missing_gene = False
        # print "Relative found in database."

    else:  # e.g., no gene present for organism of the given rank
        missing_gene = True  # still missing, keep in list

    return missing_gene


def get_xml(rank_value, gene_or_genome, giID):
    global lastRequestTime
    is_tax_correct = False
    ####################################################
    # get sequence for specified gene or genome
    ####################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=256557273&retmode=xml
    url = '/entrez/eutils/efetch.fcgi?db=nuccore&id=' + str(giID) + '&retmode=xml'
    print "searching entrez for " + str(rank_value) + " " + str(gene_or_genome)

    time_since_last_request = time.time() - lastRequestTime
    if time_since_last_request < 0.3:
        time.sleep(0.3 - time_since_last_request)
    response = get_response(url)  # get_full_taxonomy
    lastRequestTime = time.time()

    r_data = response.read()

    fasta_soup = BeautifulSoup(r_data, 'xml')
    # NOTE: BeautifulSoup does not allow you to parse dashes in xml names

    # check that rank_value is correct
    org_list = fasta_soup.find("GBSeq_taxonomy")
    # print "full organism name is " + str(org_list)
    taxonomy_list = org_list.contents[0]  # Genus species
    if rank_value in taxonomy_list:
        is_tax_correct = True

    if is_tax_correct == False:
        ##### organism is incorrect ####
        print str(rank_value) + " not found in " + str(taxonomy_list) + ' for ' + str(gene_or_genome)
        return None  # means it failed

    return fasta_soup


def get_gene_from_xml(fasta_soup, gene, db_connection):  # fasta_Soup is an object that contains the XML with the gene

    # check for synonyms
    # matching ratio is too high when comparing "NADH 1" and "NADH 2",
    # so must use exact match for longer words
    # Lists originally taken from http://www.genecards.org/ which has multiple sources including HUGO
    gene_synonym = None

    if gene == "COX1":
        gene_synonym = ["CO1", "COI", "COXI", "Cytochrome C Oxidase I",
                        "Cytochrome C Oxidase Subunit I", "MTCO1",
                        "Mitochondrially Encoded Cytochrome C Oxidase I",
                        "Cytochrome C Oxidase Polypeptide I", "EC 1.9.3.1"]

    if gene == "CYTB":
        gene_synonym = ["COB", "Cytochrome B", "cyt b",
                        "Mitochondrially Encoded Cytochrome B", "Complex III Subunit 3",
                        "Complex III Subunit III", "Cytochrome B-C1 Complex Subunit 3",
                        "Ubiquinol-Cytochrome-C Reductase Complex Cytochrome B Subunit", "MTCYB"]
    if gene == "ND2":
        # must search for exact match here, seqMatcher will match "NADH" with "NADH2"
        gene_synonym = ["NADH2", "NADH Dehydrogenase 2", "NADH Dehydrogenase Subunit 2",
                        "Complex I ND2 Subunit", "NADH-Ubiquinone Oxidoreductase Chain 2", "MTND2",
                        "Mitochondrially Encoded NADH Dehydrogenase 2", "EC 1.6.5.3"]
    if gene == "12S":  # "12S rRNA gene" AND "12S ribosomal RNA" should be captured by 12S alone
        gene_synonym = ["12S RNA", "s-rRNA", "Mitochondrially Encoded 12S RNA", "MTRNR1", "RNR1"]  # short
        # note: gene for 12S is MT-RNR1 in animals
    if gene == "16S":
        gene_synonym = ["l-rRNA", "Mitochondrially Encoded 16S RNA", "MTRNR2", "RNR2",
                        "Humanin Mitochondrial", "Formyl-Humanin", "Humanin", "HNM", "HN"]
    if gene == "trnV":
        gene_synonym = ["tRNA-Val", "TRNA Valine", "Mitochondrially Encoded TRNA Valine", "MTTV"]
        # Note: TRNA is a known synonym but there are other types of TRNA. Check for this!!!

    gene_loc = fasta_soup.find("GBSeq_locus")
    gene_locus = gene_loc.contents[0]
    # insert gene_locus into database, this is the accession number

    genus_species = fasta_soup.find("GBSeq_organism")
    genus_species = genus_species.contents[0]
    genus_species = genus_species.split()
    genus = genus_species[0]
    species = genus_species[1]
    try:
        subspecies = genus_species[2]
    except IndexError:
        subspecies = None

    seqList = fasta_soup.find("GBSeq_sequence")

    if seqList is None:  # File contains only GBSeq_contig
        print "        The xml file returned has no sequence. May be a file listing contigs."
        return None
    sequence = seqList.contents[0]  # full sequence

    # Find closest match to gene we are looking for
    best_match = 0  # Start with no match, search for best
    best_feature = None
    best_word = None
    features = fasta_soup.find_all("GBFeature")

    for feat in features:
        if feat.GBFeature_quals is None:
            continue
        # qualifiers = feat.GBFeature_quals.GBQualifier only finds first
        qualifiers = feat.find_all("GBQualifier")
        for qualifier in qualifiers:
            if qualifier.GBQualifier_value is None:
                continue
            value = qualifier.GBQualifier_value.string

            # Check entire value for a match before individual words
            if value.upper() == gene.upper():
                best_match = 1.0
                best_feature = feat
                best_word = value
                break
            if gene_synonym is not None:  # look for exact match of synonym
                for syn in gene_synonym:
                    if value.upper() == syn.upper():
                        best_match = 1.0
                        best_feature = feat
                        best_word = value
                        # print "##### gene synonym match ####" + str(syn)
                        break  # want to break out of outer loop, make function?

            # for each word in value
            # check the matching ratio
            words = value.split()
            # print "words are " + str(words)

            # if no exact match, check individual words
            for word in words:
                s = SequenceMatcher(None, word.upper(), gene.upper())
                matchRatio = s.ratio()  # 0 is no match, 1 is perfect match
                # if ratio is better that best match, save as best_feature and make new best_match score
                if matchRatio > best_match:
                    best_match = matchRatio
                    best_feature = feat
                    best_word = word

    if best_feature is None:
        print str(gene) + " not found in XML file for " + genus + " " + species
        return None

    if best_match is None:
        print str(gene) + " has no match found in XML file for " + genus + " " + species
        return None

    gene_match = best_feature.GBFeature_quals.GBQualifier.GBQualifier_value.string

    if 1 > best_match > 0.8:  # else below (prints closest) or perfect match (no print)
        print "closest match for " + gene + " is " + best_word

    if best_match <= 0.8:
        # print "############### No close match: " + str(best_word)
        # print " for " + str(gene) + " of " + str(genus) + ' ' + str(species)
        print "############### No close match: " + str(gene_match) + " for " + str(gene) + " of " + str(genus) + ' ' \
              + str(species) + ' ' + str(gene_locus)

        return None

    # get start and stop
    location = best_feature.GBFeature_location.string
    print "LOCATION IS " + str(location)
    if "join" in location or "order" in location:
        try:

            print "JOIN OR ORDER FOUND"
            print "LOCATION IS " + str(location)
            first, next = location.rsplit(",", 1)
            print "NEXT IS " + str(next)
            start, stop = next.split("..")
            stop = stop.strip(")")
            print "START IS " + str(start)
            print "STOP IS " + str(stop)
        except:
            print "EXCEPTION IN JOIN"
            return None
    else:
        start, stop = location.split("..")
    try:
        print "START IS " + str(start)
        print "STOP IS " + str(stop)
    except:
        print "FAILURE IN START OR STOP"
        return None

    partial = False
    is_complement = False

    if "complement" in start:  # check for complement() and remove
        start = start.strip("complement(")
        # stop = stop.strip(")")
        is_complement = True

    stop = stop.strip(")")

    if "<" in start or ">" in stop:  # check for &lt, &gt which indicates partial, and remove
        partial = True
        # insert partial into database
        start = start.strip("<")
        stop = stop.strip(">")

    start = int(start) - 1  # NCBI uses 1 for first item, python requires zero
    stop = int(stop)
    sequence = sequence[start:stop]

    # Complements require a==>t, g==> c, etc.
    translation_table = string.maketrans("atugcyrswkmbdhvn", "taacgryswmkvhdbn")

    if is_complement:
        translated_sequence = string.translate(str(sequence), translation_table)  # take complement
        sequence = translated_sequence[::-1]  # reverse order

    # Save whether it is a gene or cds, rna, etc
    feature_type = best_feature.GBFeature_key.string

    # save organism for header, and insert as a field or column in database
    # organism = fasta_soup.GBSeq_organism.string

    # submitToDatabase:
    # genus=>Genus, species=>Species, gene=>GeneName , gene_locus=>AccessionNumber, partial=>PartialFlag,
    # and sequence=>GeneSequence
    con = db_connection
    cur = con.cursor()

    # Check if the Genus, Species, Gene, and accession number and sequence already exist in the database
    # If it returns something, print the thing and do not update
    # otherwise update

    select = 'SELECT * FROM RefGenes WHERE Genus=? AND Species=? AND GeneName=? AND AccessionNumber=?'

    cur.execute(select, (genus, species, gene, gene_locus))
    # figure out how many rows it returned.

    rowcount = len(cur.fetchall())

    if rowcount < 0:
        print "invalid rowcount"

    if rowcount > 0:  # update
        # print "updating table for " + str(organism) + ' ' + str(gene)
        AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Superorder, Order, Suborder, \
        Infraorder, Parvorder, Superfamily, Family, Subfamily = get_full_taxonomy(genus, species)
        if AllOtherRank == None:
            return "FAILURE"
        AllOtherRank = AllOtherRank.replace("\'", "")
        update = ''' UPDATE RefGenes SET
            AllOtherRank=?, Superkingdom=?, Kingdom=?, Superphylum=?, Phylum=?, Subphylum=?, Class=?,
            Superorder=?, "Order"=?, Suborder=?, Infraorder=?, Parvorder=?,
            Superfamily=?, Family=?, Subfamily=?, Genus=? , Species=?, subspecies=?,
            GeneName=? , GeneType=?, AccessionNumber=? , PartialFlag=?, GeneSequence=?
            WHERE Genus=? AND Species=? AND GeneName=? AND AccessionNumber=? AND GeneSequence=?
            '''

        cur.execute(update, (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Superorder,
                             Order, Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily, genus, species,
                             subspecies, gene, feature_type, gene_locus, int(partial), sequence, genus, species, gene,
                             gene_locus, sequence))

    if rowcount == 0:  # insert
        AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Superorder, Order, Suborder, \
        Infraorder, Parvorder, Superfamily, Family, Subfamily = get_full_taxonomy(genus, species)
        if AllOtherRank == None:
            return "FAILURE"
        AllOtherRank = AllOtherRank.replace("\'", "")
        insert = 'INSERT INTO RefGenes (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, ' \
                 'Superorder, "Order", Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily, Genus, ' \
                 'Species, subspecies, GeneName, GeneType, AccessionNumber, PartialFLag, GeneSequence) ' \
                 'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur.execute(insert, (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Superorder,
                             Order, Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily, genus, species,
                             subspecies, gene, feature_type, gene_locus, int(partial), sequence))

    con.commit()

    return "successful"


if __name__ == '__main__':
    try:
        sample_list = sys.argv[1]  # raw data list from Alan's lab [gene, id, family, genus, species]

    except IndexError:
        # print "please specify a filename"
        print "Using default sample list."
        sample_list = "Zaher1SampleNeedReference.csv"
        # sample_list = "Zaher18SamplesNeedReference.csv"

    with open(sample_list, 'r') as sample:
        lines = sample.readlines()
        dbConnection = sqlite3.connect('./mitoDB.sqlite')

        for line in lines:
            line = line.replace('\n', '')
            cells = line.split(',')
            gene = cells[0]
            id = cells[1]
            family = cells[2]
            genus = cells[3]
            species = cells[4]
            # genus = cells[2]
            # species = cells[3]

            ########################################################
            ############## check for genus first ###################
            ########################################################
            # We already know the gene for this species is not in the database
            # since this uses an output file from MitoDBmaker
            # best to start with full genome of relative of same genus, but will disregard for now

            # missing_gene = checkDb("Genus", genus, gene)

            trans_genus, trans_species = checkSynonyms(genus, species)
            if trans_genus is not None:
                print " translated " + str(genus) + " to " + str(trans_genus)
                genus = trans_genus
                print " translated " + str(species) + " to" + str(trans_species)
                species = trans_species

            taxonomy = OrderedDict([("Genus", genus)])  # add genus
            temp_taxonomy = get_taxonomy(genus, species)
            if temp_taxonomy is None:
                print "No taxonomy for " + str(genus)
                temp_taxonomy = get_taxonomy(family, " ")
                if temp_taxonomy is None:
                    print "No taxonomy for " + str(family)
                else:
                    print "family found"

            # temp_taxonomy = get_taxonomy(family, " ")
            taxonomy.update(temp_taxonomy)  # add everything else
            for rank, rank_value in taxonomy.items():
                print "rank is " + str(rank) + " " + str(rank_value)
                if rank_value is None:
                    continue
                missing_gene = checkDb(rank, rank_value, gene)
                print missing_gene

                if missing_gene:  # still missing, search NCBI for gene
                    # gi_ids, translated_genus = get_ids(genus, gene)
                    gi_ids, translated_rank_value = get_ids(rank_value, gene)

                    if translated_rank_value is None or translated_rank_value != rank_value:
                        print " Improper value returned: " + translated_rank_value + " for " + str(rank_value)
                        # Will the translated_rank_value ever be None?
                        continue  # up to next rank?

                    else:  # we have the correct Genus (or rank)
                        if gi_ids is not None:
                            entry = None
                            for gi_id in gi_ids:
                                xml = get_xml(rank_value, gene, gi_id)  # get the xml for that genus (rank)

                                if xml is not None:  # Genome present
                                    entry = get_gene_from_xml(xml, gene, dbConnection)

                                    if entry is not None:  # successful entry
                                        break  # break from ids

                                    else:
                                        # gene name or gene location not found in xml
                                        print "gene name doesn't match or gene location not given in xml"
                                        # goes to the next gi_id

                                if entry is not None:
                                    break  # break out of rank loop
                else:
                    break