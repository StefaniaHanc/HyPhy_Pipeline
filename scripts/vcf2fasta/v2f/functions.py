import re
import sys
import time
import collections

# test for pysam
try:
    import pysam
except ImportError:
    sys.exit("pysam is not installed. Try: pip install pysam")

def getSequences(intervals, gene, ref, vcf, ploidy, phased, samples, args):

    seqs = collections.defaultdict(dict)
    
    # prep sequence dictionary
    if args.blend:
        seqs[gene] = collections.defaultdict(str)
    else:
        feat_ind = 0
        if args.feat:
            for _ in intervals[gene]:
                featname = f"{gene}_{args.feat}_{feat_ind}"
                seqs[featname] = collections.defaultdict(str)
                feat_ind += 1
        else:
            seqs[gene] = collections.defaultdict(str)
                
    if phased:
        for key in seqs.keys():
            for sample in samples:
                for i in range(ploidy):
                    seqs[key][f"{sample}_{i}"] = ''
            if args.addref:
                seqs[key]['REF_0'] = ''
    else:
        for key in seqs.keys():
            for sample in samples:
                seqs[key][sample] = ''
            if args.addref:
                seqs[key]['REF'] = ''
                
   
    # initial values before looping through VCF slice
    varsites = 0
    feat_ind = 0
    featname = gene
    codon_start = collections.defaultdict(list)
    

    # for rec in intervals[gene]:
    #     if not args.blend and args.feat:
    #         featname = f"{gene}_{args.feat}_{feat_ind}"
    #         feat_ind += 1

    # Sort exons/CDS

    records = intervals[gene]

    if args.bed:
        records = sorted(records, key=lambda x: int(x[1]))
        strand_test = records[0][5] if len(records[0]) >= 6 else "+"
    else:
        records = sorted(records, key=lambda x: int(x[3]))
        strand_test = records[0][6]

    # reverse order for minus strand
    if strand_test == "-":
        records = records[::-1]

    # process each exon/CDS block
    for rec in records:

        if not args.blend and args.feat:
            featname = f"{gene}_{args.feat}_{feat_ind}"
            feat_ind += 1

        tmpseqs = collections.defaultdict(str)

        # coordinates
        if args.bed:
            try:
                chrom = rec[0]
                start = int(rec[1])  # BED already 0-based
                end   = int(rec[2])
                strand = rec[5]
                cs = rec[4]
            except IndexError:
                chrom = rec[0]
                start = int(rec[1])
                end   = int(rec[2])
                strand = "+"
                cs = "."
        else:
            chrom = rec[0]
            start = int(rec[3]) - 1  # GFF 1-based -> 0-based
            end   = int(rec[4])
            strand = rec[6]
            cs = rec[7]

        codon_start[featname].append(cs)

        # fetch reference exon
        refseq = ref.fetch(chrom, start, end).upper()

        # initialize all samples with reference
        for sample in seqs[featname]:
            tmpseqs[sample] = refseq

        # insert SNPs
        posadd = 0

        for vrec in vcf.fetch(chrom, start, end):
        #  added
            if len(vrec.ref) != 1 or any(len(a) != 1 for a in vrec.alts):
                continue
            varsites += 1

            pos = vrec.pos - start - 1 + posadd

            alleles, max_len = getAlleles(
                vrec, ploidy, phased, args.addref
            )

            ref_len = len(vrec.ref)

            for sample in alleles:
                tmpseqs[sample] = UpdateSeq(
                    alleles, sample, pos, ref_len, tmpseqs[sample]
                )

            posadd += max_len - ref_len

        # strand handling
        if strand == "-":
            for sample in seqs[featname]:
                seqs[featname][sample] += revcomp(tmpseqs[sample])
        else:
            for sample in seqs[featname]:
                seqs[featname][sample] += tmpseqs[sample]


    #  phase trimming
    if args.inframe:
        if args.blend and codon_start[featname][0] != ".":
            if strand_test == "+" and codon_start[featname][0] != "0":
                trim = int(codon_start[featname][0])
                for k in seqs[featname]:
                    seqs[featname][k] = seqs[featname][k][trim:]


            # elif strand_test == "-" and codon_start[featname][0] != "0":
            #     trim = int(codon_start[featname][0])     
            elif strand_test == "-" and codon_start[featname][-1] != "0":
                trim = int(codon_start[featname][-1])
                for k in seqs[featname]:
                    seqs[featname][k] = seqs[featname][k][trim:]
                    
    return seqs, varsites



def getAlleles(rec, ploidy, phased, addref):
    alleles = { i[0]:i[1].alleles for i in rec.samples.items() }
    # pysam can handle ?
    segregating = list(set(sum([ [ x for x in alleles[i] ] for i in alleles.keys() ], [])))
    
    if addref:
        segregating = list(set(segregating+[rec.ref]))
        
    max_len = max([ len(i) for i in segregating if i is not None ])
    dict_expanded = { i:(i + '-' * (max_len - len(i))) for i in segregating if i is not None }
    alleles_expanded = { i:[dict_expanded[j] for j in alleles[i]] for i in alleles.keys() if alleles[i][0] is not None }
    
    if addref:
        alleles_expanded['REF'] = [dict_expanded[rec.ref]]
        
    if None in segregating:
        dict_expanded[''] = '-' * max_len
    # if None in segregating:
    #     dict_expanded[''] = rec.ref + '-' * (max_len - len(rec.ref))        
        alleles_missing = { i:[dict_expanded[''] for j in range(ploidy)] for i in alleles.keys() if alleles[i][0] is None }
        for i in alleles_missing.keys(): 
            alleles_expanded[i] = alleles_missing[i]
            
    if phased:
        alleles_expanded = makePhased(alleles_expanded)
    else:
        alleles_expanded = {key:getIUPAC(alleles_expanded[key]) for key in alleles_expanded.keys()}
        ## if want to get first allele only : 
        #alleles_expanded = {key:x[0] for key,x in alleles_expanded.items()}

    return alleles_expanded, max_len

def makePhased(alleles):
    alleles_phased = {}
    for samp in alleles.keys():
        for i in range(len(alleles[samp])):
            alleles_phased[f"{samp}_{i}"] = alleles[samp][i]
    return alleles_phased

def UpdateSeq(alleles, samp, pos, ref_len, seq):
    return seq[:pos] + alleles[samp] + seq[pos+ref_len:]

def getIUPAC(x):
    '''
    Collapses two or more alleles into a single IUPAC string
    '''
    if len(x) == 1 or x[0][0] == '?':
        return x[0]
    elif len(list(set(x))) == 1:
        return x[0]
    else:
        iupacd = ''
        for i in range(len(x[0])):
            nuc = list(set([ y[i] for y in x ]))
            
            if len(nuc) == 1:
                iupacd += nuc[0]
            else:
                nuc = [ j for j in nuc if j != '-' ]
                
                if len(nuc) == 0:
                    iupacd += '-'
                elif len(nuc) == 1:
                    iupacd += nuc[0]
                elif len(nuc) == 2:
                    if   'A' in nuc and 'G' in nuc: iupacd += 'R'
                    elif 'A' in nuc and 'T' in nuc: iupacd += 'W'
                    elif 'A' in nuc and 'C' in nuc: iupacd += 'M'
                    elif 'C' in nuc and 'T' in nuc: iupacd += 'Y'
                    elif 'C' in nuc and 'G' in nuc: iupacd += 'S'
                    elif 'G' in nuc and 'T' in nuc: iupacd += 'K'
                    else:                           iupacd += '?'
                elif len(nuc) == 3:
                    if   'A' in nuc and 'T' in nuc and 'C' in nuc: iupacd += 'H'
                    elif 'A' in nuc and 'T' in nuc and 'G' in nuc: iupacd += 'D'
                    elif 'G' in nuc and 'T' in nuc and 'C' in nuc: iupacd += 'B'
                    elif 'A' in nuc and 'G' in nuc and 'C' in nuc: iupacd += 'V'
                    else:                                          iupacd += '?'
                elif len(nuc) == 4:
                    if 'A' in nuc and 'G' in nuc and 'C' in nuc and 'T' in nuc:
                        iupacd += 'N'
                    else:
                        iupacd += '?'
                else:
                    iupacd += '?'
        return iupacd

def printFasta(seqs, out, remove_gap):
    '''
    writes FASTA to file
    '''
    for head in seqs.keys():
        out.write(">" + head + "\n")
        if not remove_gap:
            out.write(seqs[head] + "\n")
        else:
            out.write(seqs[head].replace('-', '') + "\n")

def getFeature(file):
    features = collections.defaultdict()
    with open(file, "r") as f:
        for line in f:
            if line[0] != "#":
                fields = line.rstrip().split("\t")
                features[fields[2]] = None
    return list(features.keys())

def getGeneNames(file, format):
    geneNames = collections.defaultdict()
    with open(file, "r") as f:
        if format == "gff":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGFF(fields[8])
                    if last.get('Name'):
                        geneNames[last['Name']] = None
                    elif last.get('Parent'):
                        geneNames[last['Parent']] = None
                    elif last.get('ID'):
                        geneNames[last['ID']] = None
        elif format == "gtf":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGTF(fields[8])
                    if last.get('transcript_id'):
                        geneNames[last['transcript_id']] = None
                    elif last.get('gene_id'):
                        geneNames[last['gene_id']] = None
    return list(geneNames.keys())

def processGeneNameGFF(lastfield):
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if "=" in i:
            x = i.split("=")
            last[re.sub("\"| ","",x[0])] = re.sub("\"| ","",x[1])
    return last

def processGeneNameGTF(lastfield):
    last = collections.defaultdict()
    for i in lastfield.split(";"):
        if " " in i:
            x = i.split(" ")
            last[x[0]] = re.sub("\"| ","",x[1])
    return last

def ReadBED(file):
    bed = {}
    with open(file) as f:
        i = 1
        for line in f:
            if not line.strip().startswith('#'):
                tmp = line.strip().split('\t')
                if len(tmp) == 3:
                    bed[f'g{i}'] = tmp
                    i += 1
                elif len(tmp) == 6:
                    if tmp[4] in ['.', '0', '1', '2']:
                        bed[tmp[3]] = tmp
                    else:
                        sys.exit('The fifth column in Bed6 file is considered as phase value here\nand therefore should be the one from (".", "0", "1", "2")')
                else:
                    sys.exit('please input Bed3 or Bed6')
    return bed

def ReadGFF(file, parser):
    if file.split('.')[-1].lower() in ["gff", "gff3"]:
        format = "gff"
    elif file.split('.')[-1].lower() == "gtf":
        format = "gtf"
    else:
        print('Cannot figure out GFF/GTF format. File should end with .gff or .gtf')
        sys.exit(parser.print_help())
        
    geneNames = getGeneNames(file, format)
    features  = getFeature(file)
    gff = collections.defaultdict(lambda: collections.defaultdict(list))
    
    with open(file, "r") as f:
        if format == "gff":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGFF(fields[8])
                    if last.get('Name'):
                        gff[last['Name']][fields[2]].append(fields)
                    elif last.get('Parent'):
                        gff[last['Parent']][fields[2]].append(fields)
                    else:
                        gff[last['ID']][fields[2]].append(fields)
        elif format == "gtf":
            for line in f:
                if line[0] != "#":
                    fields = line.rstrip().split("\t")
                    last = processGeneNameGTF(fields[8])
                    if last.get('transcript_id'):
                        gff[last['transcript_id']][fields[2]].append(fields)
                    elif last.get('gene_id'):
                        gff[last['gene_id']][fields[2]].append(fields)
    return gff

def filterFeatureInGFF(gff, feat):
    filtered_gff = collections.defaultdict()
    for gene in gff.keys():
        if len(gff[gene][feat]) != 0:
            filtered_gff[gene] = gff[gene][feat]
    return filtered_gff

def getPloidy(vcf):
    var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    p = [ len(v.get('GT')) for v in var if v.get('GT')[0] is not None ]
    return p[0]

def getPhased(vcf):
    var = [ y for x,y in next(vcf.fetch()).samples.items() ]
    p = any([ not v.phased for v in var ])
    return not p

def revcomp(seq):
    table = str.maketrans(
        "ACGTRYMKSWBDHVN?-",
        "TGCAYRKMWSVHDBN?-"
    )
    return seq[::-1].translate(table)


def make_progress_bar(rec, total, t1, width):
    i = (rec/total * 100) % 100
    if i != 0:
        plus = "+" * int(i * (width/100))
        dots = "." * (width-int(i * width/100))
    else:
        plus = "+" * width
        dots = ""
        i = 100
    t2 = time.time()
    elapsed = t2-t1

    return f"[{plus}{dots}] {i:5.2f}% {elapsed:7.2f} s"





