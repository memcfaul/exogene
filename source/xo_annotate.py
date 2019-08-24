from urllib import request
from xo_ape import ape
from collections import defaultdict as dd
from regex import compile as re_compile,IGNORECASE as re_ig
from pathlib import Path as plp

#path to primer list excel book

class annotate(object):
    def __init__(self,symbol,species="danio_rerio",splice="cds", peptide=False,\
            cDNA=True,gDNA=True,flank=500,do_open=True,save_path='',auto_start=True):
        """splice: "cds", "cdna", or "all" 
        """
        
        self.save_path = save_path
        self.splice = splice.lower()
        self.do_open = do_open
        #gene name/symbol
        self.species = species
        self.symbol = symbol 
        #mRNA transcript IDs and info
        self.transcripts = []
        self._transcript_info = []
        self.best_transcript = None       
        self.sequences = dd(lambda x=dd:x(list))
        self.genomic = None
        self.cdna = {}
        self.peptide = {}
        # controls
        self.flank = flank
        self.get_gDNA = gDNA
        self.get_cDNA = cDNA
        self.get_peptide = peptide
        if auto_start:
            self._initialize()
    
    def _initialize(self):
        a = self._fasta_clean(self.fetch_gene(self.symbol,flank=self.flank,get_cDNA=True))
        if not a:
            raise IndexError("Gene '{}' could not be found for organism '{}'".format(self.symbol,self.species))
        b = self._parse_fasta(a)
        self._parse_sorter(b)
        self._ape_annotate()  
             
    def prime_and_write(self,path=None,do_open=False):
        if not path:
            path = self.save_path
        self.import_primers()
        self.write_ape(path=path,do_open=do_open)
        
    def fetch_gene(self,symbol=None,flank = 500,output = False,get_peptide = None,get_gDNA=None,get_cDNA = None):
        gdna = ''
        cdna = ''
        peptide = ''
        if symbol is None:
            symbol = self.symbol
        if get_cDNA or (get_cDNA is None and self.get_cDNA):
            cdna = 'param=cdna;param=coding;param=utr5;param=utr3;param=exon;'
        if get_gDNA or (get_gDNA is None and self.get_gDNA):
            gdna = 'genomic=unmasked;'
        if get_peptide or (get_peptide is None and self.get_peptide):
            peptide = 'param=peptide;'.format(self.best_transcript)
        url = '''http://uswest.ensembl.org/{0}/Export/Output/Gene?db=core;flank3_display={1};flank5_display={1};g={2};output=fasta;strand=feature;{3}{4}{5}_format=Text'''.format(self.species,flank,symbol,cdna,gdna,peptide)
        fasta = request.urlopen(url).read()
        fasta = fasta.decode('utf-8')
        if output:   
            with open('{}.txt'.format(self.symbol),'w') as out:
                out.write(fasta)
        return [''.join(['>',h]) for h in fasta.split('>') if h] 
    
    def _fasta_clean(self,fasta):
        cleaned = []
        for i in fasta:
            splits = i.split('\r\n')
            if len(splits)>1 and splits[1]:
                head = splits[0]
                seq = ''.join(splits[1:])
                cleaned.append('{}\n{}'.format(head,seq.strip()))
        return "".join(cleaned)
                           
        
    def _parse_fasta(self,fasta):
        pat = re_compile(r'>(?P<hGene>[\w\d\.]+ ){0,1}(?P<hSplice>[\w\.]+)[ :]{1,2}(?P<id>[\w\.]+) (?P<type>[\w\.]+)[ :]{1,2}(?P<info>[\w\.]+)[ :]{0,1}(?P<chrom>\w+:\w+:*\d*:*[- 0-9]*){0,1}.*?\n(?P<seq>[^>]+)',flags=re_ig)
        fasta_parse = pat.finditer(fasta)
        return fasta_parse
    
    def _parse_sorter(self,fasta_parse):
        params = ['cds','utr5', 'exon', 'utr3' ,'cdna','dna','peptide']   
        for i in fasta_parse:
            splice = i.group('id')
            kind = i.group('type').replace(' ','') 
            if kind == params[6]:
                self.peptide[i.group('info')] = i.group('seq')      
            elif splice in ('chromosome','primary_assembly','dna'):
                if kind in ('chromosome','primary_assembly'):
                    info = i.group('chrom').split(':')
                    self.sequences["dna"]["sequence"] = i.group('seq')
                    self.sequences["dna"]["version"] = 'Genome Version = {}'.format(i.group('info'))
                    self.sequences["dna"]["chromosome"] = 'Chromosome = {}'.format(info[0])
                    self.sequences["dna"]["coordinates"] = 'Coordinates = {},{}; Orientation = {}'.format(info[1],info[2],info[3])
            else:          
                if 'exon' in i.group('info'):
                    kind = 'exon'     
                self.sequences[splice][kind].append(i.group('seq'))
    
        if self.splice == "cdna":
            self.best_transcript = self._best_transcript("cdna")
        else:
            self.best_transcript = self._best_transcript("cds")
        self.sequences["dna"]["annotation"] = "Annotated by transcript {}"\
            .format(self.best_transcript)
    
    def _best_transcript(self,par):
        if par not in ("cds","cdna"):
            seqs = [(k,len(v.get(par,[""])[0])) for k,v in self.sequences.items()]
            return max(seqs,key=lambda x: len(x[1]))[0]
        else:
            primary = [1,2][par=="cdna"]
            secondary = [2,1][primary-1]
            seqs = [(k,len(v.get("cds",[""])[0]),len(v.get("cdna",[""])[0])) for k,v in self.sequences.items()]
            p_max_val = max(seqs,key=lambda x:x[primary])[primary]
            primary_max = [i for i in seqs if i[primary]>=p_max_val]
            return max(primary_max,key=lambda x: x[secondary])[0]           
              
    def _ape_annotate(self):
        if self.get_gDNA:
            self.genomic = ape(name="{}_gDNA".format(self.symbol),\
                sequence=self.sequences["dna"].pop("sequence"))
            [self.genomic.add_comment(self.sequences["dna"].get(i)) \
                for i in self.sequences["dna"].keys()]
            self.sequences.pop("dna")
        if self.splice == "all":
            self.cdna.update(**{k:ape(name="{}_{}_cDNA".format(self.symbol,k),\
                sequence=v.get("cdna")[0])
            for k,v in self.sequences.items()})
        else:
            k,v =  self.best_transcript,\
                self.sequences[self.best_transcript]
            self.cdna.update(**{k:ape(name="{}_{}_cDNA".format(self.symbol,k),\
                sequence=v.get("cdna")[0])})
                
        for trans in self.cdna.keys():
            exon_count = 1
            params = ['cds','exon','utr5','utr3']
            color_dict = {
            'exon':['#FFCC99', '#5EF1F2', '#426600'],
            'utr5':'#FFA8BB',
            'utr3':'#FFA8BB',
            'cds':'#0075DC'
            }
            color_pos = 0
            for par in params:  
                dict_shot = self.sequences[trans]
                features = dict_shot.get(par,'')
                if features:
                    for seq in features:
                        name = par
                        try:
                            if par == 'exon':
                                name += str(exon_count)
                                exon_count += 1                            
                                c_color = color_dict[par][color_pos]
                                g_color = color_dict[par][0]
                                color_pos+=1
                                if color_pos == len(color_dict[par]):
                                    color_pos =0
                            else:
                                c_color = color_dict[par]
                                g_color = color_dict[par]
                            self.cdna[trans].add_feature(seq,name = name,fwd_color = c_color,feature_type=par)
                        except self.cdna[trans].FeatureError:
                            try:
                                self.cdna[trans].add_feature(seq,name = name,fwd_color = c_color,feature_type=par,mismatch=int(len(seq)*0.02))
                            except self.cdna[trans].FeatureError as e:
                                print(e,'- cdna[{}].add_feature({})'.format(trans,par))  
                        
                        if trans == self.best_transcript:                                             
                            try:
                                if self.genomic:  
                                    self.genomic.add_feature(seq,name = name,fwd_color = g_color,feature_type=par)
                            except self.genomic.FeatureError as e:
                                if e.arg!=3:
                                    print(e, '- genomic.add_feature({})'.format(par))
                                if 'utr' not in par:
                                    break
                                try:    
                                    seq_ape = ape(seq)
                                    for ex in dict_shot['exon']:
                                        if len(seq_ape.sequence) >= len(ex):
                                            seq_search = seq_ape.degenerate_search(ex,int(len(ex)*0.02),search_reverse=False)
                                            if seq_search:
                                                g_search = self.genomic.degenerate_search(ex,int(len(ex)*0.02),search_reverse=False)
                                                self.genomic.add_feature(*min(g_search),name = name,fwd_color = g_color,feature_type=par)
                                                seq_ape.sequence = seq_ape.sequence[max(min(seq_search)):]
                                            
                                        else:
                                            ex_ape = ape(ex)
                                            seq_copy = seq_ape.sequence
                                            search_seq = ex_ape.degenerate_search(seq_copy,int(len(seq_copy)*0.02),search_reverse=False)
                                            g_search = self.genomic.degenerate_search(ex,int(len(ex)*0.02),search_reverse=False)
                                            if search_seq and g_search:
                                                start = min(min(g_search))
                                                end = start+max(min(search_seq))-1
                                                self.genomic.add_feature(start,end,name = name,fwd_color = g_color,feature_type=par)
                                                break                                    
                                except self.genomic.FeatureError as e:
                                    print(e, '- genomic.add_feature({})'.format(par))
   
    def _add_feat(self,ape_obj,primer,name,feature_type,fwd_color,rev_color,mismatch,partial_primer,sequence_splits):
            try:
                ape_obj.add_feature(seq_start=primer,name=name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch)
            except ape_obj.FeatureError:
                if partial_primer:
                    try:
                        ape_obj.add_partial_feature(sequence=primer,sequence_splits=sequence_splits,name=name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch)
                    except ape_obj.FeatureError:
                        pass
                        
    def determine_file_type(self,path):
        ext = str(path).rsplit(".")[-1]
        if ext in ["xlsx","xls","xltx","xlt"]:
            return "excel"
        else:
            return "txt"
            
    def get_primer_seqs(self,path,primer_name_index=0,primer_sequence_index=1,sheet_index=0):
        ape = ape()
        file_type = self.determine_file_type(path)
        if file_type=="excel":
            primers = ape.get_primers_excel(path,sheet_index,primer_name_index,primer_sequence_index)
        else:
            primers = ape.get_primers_txt(path,primer_name_index,primer_sequence_index)  
        return primers   
        

    def import_primers(self,path='primers.xlsx',feature_type='primer',fwd_color='#2BCE48',rev_color ='#FF5005',mismatch=0,sheet_index=0,primer_name_index=0,primer_sequence_index=1,partial_primer=False,sequence_splits='default'):             
        primers = self.get_primer_seqs(path,sheet_index=sheet_index,primer_name_index=primer_name_index,primer_sequence_index=primer_sequence_index)
        self.add_primers(primers,feature_type,fwd_color,rev_color,mismatch,partial_primer,sequence_splits)
        
    def add_primers(self,primers,feature_type='primer',fwd_color='#2BCE48',rev_color ='#FF5005',mismatch=0,partial_primer=False,sequence_splits='default'):             
        for item in primers.items():
            primer_name = item[0]            
            primer = item[1]
            if primer and primer_name:
                if self.get_cDNA:
                    for trans in self.cdna.values():
                        self._add_feat(trans,primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch,partial_primer=partial_primer,sequence_splits=sequence_splits)
                if self.get_gDNA:
                    self._add_feat(self.genomic,primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch,partial_primer=partial_primer,sequence_splits=sequence_splits)
    
    def import_crisprs(self,path='crisprs.xlsx',feature_type='Cas9_sgRNA',fwd_color='#8000FF',rev_color ='#FF00FF',pam_color='#FFFF00',sheet_index=2,primer_name_index=0,primer_sequence_index=1,sequence_splits='default',pam='NGG'):
        primers = self.get_primer_seqs(path,primer_name_index,primer_sequence_index)
        ape = ape()
        for item in primers.items():
            primer_name = item[0] 
            primer = item[1]
            primer_strip = ape.strip_features(sequence=primer,sequence_splits=sequence_splits)
            primer = primer_strip
            if primer and primer_name:
                if self.get_cDNA:
                    for trans in self.sequences.keys():
                        self.cdna[trans]._add_crispr(primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,pam=pam,pam_color=pam_color)
                if self.get_gDNA:
                    self.genomic._add_crispr(primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,pam=pam,pam_color=pam_color)
    
    def write_ape(self):            
        path=self.save_path            
        if self.get_cDNA:
            for i in self.cdna.keys():
                self.cdna[i].create_ape_file("{}_{}_cDNA_{}".format(self.symbol,i,self.species),path,do_open=self.do_open)
        if self.get_gDNA:
            self.genomic.create_ape_file('{}_gDNA_{}'.format(self.symbol,self.species),path,do_open=self.do_open)
        if self.get_peptide:
            file_name = "{}_protein_{}.txt".format(self.symbol,self.species)
            pf = '{}/{}'.format(path,file_name).replace('//','/')
            with open(pf,'w') as out:
                for name,seq in self.peptide.items():
                    out.write(">{}_{}_protein\n".format(self.symbol,name))
                    for i in range(0,len(seq),100): 
                        out.write("{}".format(seq[i:i+100]))
                    out.write("\n\n")
            if self.do_open:
                ape = ape()
                ape.open_ape(pf,"open")
                
                                                                                   
if __name__ == '__main__':
    gene = annotate("etv5a",species = "danio_rerio",save_path=plp.home() /"Desktop" )
    gene.write_ape()
        