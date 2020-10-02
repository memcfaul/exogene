from os import remove,path as os_path
from pathlib import Path as plp
import random
import xlrd
from subprocess import run
from collections import OrderedDict as od
from fuzzywuzzy import process as fwpro
from regex import compile as r_compile, IGNORECASE as ig, sub as r_sub,findall
from sys import platform
class ape(object):
    def __init__(self,sequence='',name='DNA',validate_sequence=False):
        self.parial_dict = {
            't7':(r_compile('\w*?'+self.deg_replace('TAATACGACTCACTATAGN'),flags=ig),'NN'),
            'sp6':(r_compile('\w*?'+self.deg_replace('ATTTAGGTGACACTATAGN'),flags=ig),'N'),
            't3':(r_compile('\w*?'+self.deg_replace('AATTAACCCTCACTAAAGN'),flags=ig),'N'),
            'sgRNA_overlap':(r_compile(self.deg_replace('GTTTTAGAGCTAGAAN')+'*',flags=ig),''),
            'attB1R_primer': (r_compile('\w*?'+self.deg_replace('ACTGCTTTTTTGTACAAACTTG'),flags=ig),''),
            'attB1_primer': (r_compile('\w*?'+self.deg_replace('ACAAGTTTGTACAAAAAAGCAGGCT'),flags=ig),''),
            'attB2R_primer': (r_compile('\w*?'+self.deg_replace('ACCACTTTGTACAAGAAAGCTGGGT'),flags=ig),''),
            'attB2_primer': (r_compile('\w*?'+self.deg_replace('ACAGCTTTCTTGTACAAAGTGG'),flags=ig),''),
            'attB3R_primer': (r_compile('\w*?'+self.deg_replace('ACAACTTTGTATAATAAAGTTG'),flags=ig),''),
            'attB4_primer': (r_compile('\w*?'+self.deg_replace('ACAACTTTGTATAGAAAAGTTG'),flags=ig),'')}
        self.sequence = r_sub('[\W]','',sequence)
        self.valid_characters = {'A','T','C','G','N','H','D','V','B','K','M','Y','R','W','S'}
        if validate_sequence and not self.validate_sequence(self.sequence):
            raise self.SequenceError("Sequence contains invalid characters")
        self.name = name
        self.file_name = None
        self._save_path = plp.cwd()
        self.features = od()
        self.comments = []
        self.length = len(self.sequence)
        self._ape_header1 = 'LOCUS\nACCESSION\nVERSION\n'
        self._ape_header2 = 'COMMENT     ApEinfo:methylated:1\nFEATURES{}Location/Qualifiers\n'.format(' '*13)                 
        self._ape_feature = '     {7}{5}{0[0]}..{0[1]}{6}\n{4}/label={1}\n{4}/ApEinfo_fwdcolor="{2}"\n{4}/ApEinfo_revcolor="{3}"\n{4}/ApEinfo_graphicformat="arrow_data {{{{0 1 2 0 0 -1}} {{}} 0}}\n{4}width 5 offset 0"\n'
        self.color_dict = od(
        [('Cayenne', '#800000'), ('Asparagus', '#808000'), ('Clover', '#008000'), ('Teal', '#008080'), 
        ('Midnight', '#000080'), ('Plum', '#800080'), ('Tin', '#7F7F7F'), ('Nickel', '#808080'), 
        ('Mocha', '#804000'), ('Fern', '#408000'), ('Moss', '#008040'), ('Ocean', '#004080'), 
        ('Eggplant', '#400080'), ('Maroon', '#800040'), ('Steel', '#666666'), ('Aluminum', '#999999'), 
        ('Marascino', '#FF0000'), ('Lemon', '#FFFF00'), ('Spring', '#00FF00'), ('Turquoise', '#00FFFF'), 
        ('Blueberry', '#0000FF'), ('Magenta', '#FF00FF'), ('Iron', '#4C4C4C'), ('Magnesium', '#B3B3B3'), 
        ('Tangerine', '#FF8000'), ('Lime', '#80FF00'), ('SeaFoam', '#00FF80'), ('Aqua', '#0080FF'), 
        ('Grape', '#8000FF'), ('Strawberry', '#FF0080'), ('Tungsten', '#333333'), ('Silver', '#CCCCCC'), 
        ('Salmon', '#FF6666'), ('Banana', '#FFFF66'), ('Flora', '#66FF66'), ('Ice', '#66FFFF'), 
        ('Orchid', '#6666FF'), ('Bubblegum', '#FF66FF'), ('Lead', '#191919'), ('Mercury', '#E6E6E6'), 
        ('Cantaloupe', '#FFCC66'), ('Honeydew', '#CCFF66'), ('Spindrift', '#66FFCC'), ('Sky', '#66CCFF'), 
        ('Lavender', '#CC66FF'), ('Carnation', '#FF6FCF'), ('Licorice', '#000000'), ('Snow', '#FFFFFF'),
        ('black', 'black'), ('blue', 'blue'), ('brown', 'brown'), ('cyan', 'cyan'), ('green', 'green'), 
        ('magenta', 'magenta'), ('orange', 'orange'), ('purple', 'purple'), ('red', 'red'), 
        ('yellow', 'yellow'), ('white', 'white')]
        )    
        self._color_series = ['#FF0000', '#FFFF00', '#00FFFF', '#80FF00', '#0000FF', '#FF00FF', '#0080FF', '#FF8000', '#00FF80', '#8000FF', '#FF0080']
        self._color_pos = 0
    class SequenceError(Exception):
        pass           
    class FeatureError(Exception):
        def __init__(self,arg=''):
            self.arg = arg
            self.fail_dict = {1 : ', input wrong format',
            2 : ', out of range',
            3: ', not found in sequence'}
        def __str__(self):
            return 'Could not add feature {}'.format(self.fail_dict.get(self.arg,self.arg)) 
    def validate_sequence(self,seq):
        return all((i.upper() in self.valid_characters for i in seq))
    def __repr__(self):
        return '"ape({})"'.format(self.name)                  
    def __str__(self):
        return "{}: {}bp".format(self.name,self.length)       
    def __len__(self):
        return len(self.sequence)     
    def __nonzero__(self):
        return bool(self.sequence)
                  
    class results(od):
        def __str__(self):
            prnt = []
            for i in self.iteritems():
                prnt.append("{}:{}".format(i[0],i[1]))
            return "\n".join(prnt)
        
    def deg_replace(self,subs):
        deg_substring = ''
        replace_dict = {'B':'[CGTBSKY]','D':'[AGTDRWK]','H':'[ACTHMYW]','K':'[GTK]','M':'[ACM]','N':'[ACGTBDHKMNRSVWY]','R':'[AGR]','S':'[CGS]','V':'[ACGVMSR]','W':'[ATW]','Y':'[CTY]',}
        change = True
        for i in subs:
            if i == '[':
                change = False
            elif i == ']':
                change = True
            if change: 
                i = replace_dict.get(i.upper(),i)
            deg_substring += i
        return deg_substring   
             
    def degenerate_search(self,substring,mismatch,overlap=False,search_reverse=True,correction=(1,0)):
        def sub_search(substring,string,mismatch,overlap,correction):
            match = []
            p = r"(?=({})){{s<={}}}".format(substring,mismatch)
            pat = r_compile(p,ig)
            for m in pat.finditer(string):
                match.append((m.start()+correction[0],m.start()+correction[1]+len(m.group(1))))
            return match   

        deg_substring = self.deg_replace(substring)
        match = sub_search(deg_substring,self.sequence,mismatch,overlap,correction=correction)
        
        if search_reverse:
            r_substring = self.reverse_complement(substring)
            deg_r_substring = self.deg_replace(r_substring)
            reverse = sub_search(deg_r_substring,self.sequence,mismatch,overlap,correction=correction)
            reverse = [i[::-1] for i in reverse]
            match = match+reverse
            match.sort()
        if match:
            return match
            
    def strip_features(self,sequence,sequence_splits='default',min_len_pre=30,min_len_post=16):
        if sequence_splits == 'default':
            sequence_splits = self.parial_dict.values()
        elif not isinstance(sequence_splits,list):
            raise self.FeatureError('Partial sequence collection must be list or dictionary')
        if len(sequence)<=min_len_pre:
            return sequence
        for i in sequence_splits:
            reverse_sequence = self.reverse_complement(sequence)
            sequence = r_sub(i[0],i[1],sequence)
            r_sequence = r_sub(i[0],i[1],reverse_sequence)
            sequence = min((sequence,r_sequence),key=lambda x: len(x))
        if len(sequence)>=min_len_post: 
            return sequence
        return None 
                
    def add_partial_feature(self,sequence,sequence_splits='default',name='partial_feature',mismatch=0,feature_type='partial_feature',fwd_color='',rev_color='',fade=True):
        orig_sequence = sequence
        sequence = self.strip_features(sequence,sequence_splits)
        loc = self.add_feature(sequence,name=name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch) 
        if len(sequence)<=len(orig_sequence):
            self.add_comment('{} full sequence is {}'.format(name,orig_sequence.upper()))
        return loc  
                   
    def reverse_complement(self,seq,remove=True,verbose=True):
        def non_letter(sequence,letter,talk=True):
            if talk:
                print('Non-ATG letter "{}" detected'.format(letter),repr(str(sequence)))
                    
        seq = seq.replace(' ','')
        seq_dict_small = {'A':'T','C':'G','N':'N','H':'D','V':'B','K':'M','Y':'R','W':'S'}
        seq_dict = {}
        for i in seq_dict_small.items():
            seq_dict[i[0]] = i[1]
            seq_dict[i[0].lower()] = i[1].lower()
            seq_dict[i[1]] = i[0]
            seq_dict[i[1].lower()] = i[0].lower()
        rev_comp = ''
        for i in seq[::-1]:
            letter =  seq_dict.get(i,'x')
            if letter =='x':
                non_letter(seq,i,verbose)
                if not remove:
                    rev_comp += i
            else:
                rev_comp += letter                
        return rev_comp
        
    def get_best_match(self,query,choice_list=None,cutoff=90,verbose=True,len_error=1,get_scores=False):
        if choice_list is None:
            choice_list = self.features.keys()
        matches = fwpro.extract(query,choice_list)
        if matches:
            best = matches[0]
            if best[1] >= cutoff and len(best[0])-len_error <= len(query) <= len(best[0])+len_error:
                return best[0]
        if verbose:
            print("Could not find '{}'".format(query))
            poss = [i[0] for i in matches if i[1] >= 45]
            if poss:
                print(("Did you mean "+"{}, "*len(poss))[:-2].format(*poss))
        if get_scores:
            return matches
 
    def get_seq(self,start,stop=None):
        reverse = False
        if stop is None:
            feat_name = self.get_best_match(start)
            if feat_name is None:
                raise self.FeatureError(3) 
            start,stop,rev = self.features.get(feat_name)[1:4]
            reverse = rev == -1
            
        if isinstance(start,int) and isinstance(stop,int):
            if max(start,stop) <= len(self.sequence) and min(start,stop) > 0 :
                if start>stop:
                    reverse = True
                    start,stop = stop,start 
                seq = self.sequence[start-1:stop]
                if reverse:
                    seq = self.reverse_complement(seq,verbose=False)
                return seq        
            else:
                raise self.FeatureError(2)      
        else:
            raise self.FeatureError(1)
              
  
    def get_primers_txt(self,path,primer_name_index=0,primer_sequence_index=1,sep_regex = "[\t;,]"):
        primers = od()
        sep = r_compile("{}".format(sep_regex))
        with open(path,'r') as f:
            for line in f:
                splits = sep.split(line)
                if len(splits)>1:
                    name = splits[primer_name_index]
                    seq = splits[primer_sequence_index]
                    primers[name.strip()] = seq.strip()
        return primers
                 
    def determine_file_type(self,path):
        ext = str(path).rsplit(".")[-1]
        if ext in ["xlsx","xls","xltx","xlt"]:
            return "excel"
        else:
            return "txt"    
            
    def get_primers_excel(self,path,sheet_index=0,primer_name_index=0,primer_sequence_index=1):
        book = xlrd.open_workbook(path)
        sheet = book.sheet_by_index(sheet_index)
        rows = sheet.nrows
        primers = od()
        for i in range(rows):
            primer_name = sheet.cell_value(i,primer_name_index).strip()
            primer = r_sub('[\W]','',sheet.cell_value(i,primer_sequence_index))
            if primer and self.validate_sequence(primer):
                primers[primer_name] = primer
        return primers
        
    def _tint_hex(self,hex_color,offset=50,exception=True):
        hex_split = findall('\w\w',hex_color)   
        if len(hex_split)==3:
            hex_split = map(lambda x: int(x,16)+offset,hex_split)
            hex_split = map(lambda x: hex(min((255,max(0,x))))[2:].zfill(2),hex_split)
            return "#"+"".join(hex_split)
        if exception:
            raise self.FeatureError('{} is not valid hex code'.format(hex_color))
        else:
            return hex_color
                
    def import_primers(self,path,feature_type='primer',fwd_color='#2BCE48',rev_color='#FF5005',sheet_index=0,primer_name_index=0,primer_sequence_index=1,mismatch=0,partial_primer=False,sequence_splits='default',fade=True,verbose=False,circular=False):
        primers = self.get_primers_excel(path,sheet_index,primer_name_index,primer_sequence_index)
        return self.add_primers(primers,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,sheet_index=sheet_index,primer_name_index=primer_name_index,primer_sequence_index=primer_sequence_index,mismatch=mismatch,partial_primer=partial_primer,sequence_splits=sequence_splits,fade=fade,verbose=verbose,circular=circular)
    
    def add_primers(self,primers,feature_type='primer',fwd_color='#2BCE48',rev_color='#FF5005',sheet_index=0,primer_name_index=0,primer_sequence_index=1,mismatch=0,partial_primer=False,sequence_splits='default',fade=True,verbose=False,circular=False):                  
        aligned = []
        for item in primers.items():
            primer_name = item[0]            
            primer = item[1]
            if primer and primer_name:
                try:
                    alignment = self.add_feature(primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,mismatch=mismatch )
                    if alignment:
                        for match in alignment: 
                            aligned.append([primer_name,match[1],primer,match[-1]])  
                except self.FeatureError:
                    if partial_primer:
                        try:
                            alignment = self.add_partial_feature(sequence=primer,sequence_splits=sequence_splits,name=primer_name,mismatch=mismatch,feature_type=feature_type,fwd_color=fwd_color,rev_color=rev_color,fade=fade)
                            if alignment:
                                for match in alignment: 
                                    aligned.append([primer_name,match[1],primer,match[-1]])                             
                        except self.FeatureError:
                           pass
        return aligned  
        
    def _add_crispr(self,primer,name,feature_type,fwd_color,rev_color,pam_color,pam,mismatch=0):
        loc = self.degenerate_search(substring=primer+pam,mismatch=mismatch)
        pam_len = len(pam)
        if loc:
            for i in loc:
                i = self.orientation(name,*i)
                ori = i[-1]
                start = i[0:3][ori]
                pam_stop = i[0:3][ori*-1]
                stop = pam_stop-(pam_len*ori)
                pam_start = pam_stop-((pam_len-1)*ori)
                self.add_feature(start,stop,name,feature_type,fwd_color,rev_color,mismatch=mismatch)
                self.add_feature(pam_start,pam_stop,name+"_pam",'PAM',pam_color,pam_color)
                
    def import_crisprs_excel(self,path,feature_type='Cas9_sgRNA',fwd_color='#8000FF',rev_color ='#FF00FF',pam='NGG',pam_color='#FFFF00',mismatch=0,sheet_index=3,primer_name_index=0,primer_sequence_index=1,sequence_splits='default'): 
        primers = self.get_primers_excel(path,sheet_index,primer_name_index,primer_sequence_index)
        for item in primers.items():
            primer_name = item[0] 
            primer = item[1]           
            primer_strip = self.strip_features(sequence=primer,sequence_splits=sequence_splits)
            primer = primer_strip
            if primer:
                self._add_crispr(primer,name=primer_name,feature_type=feature_type,fwd_color=fwd_color,mismatch=mismatch,rev_color=rev_color,pam=pam,pam_color=pam_color)
           
                 
    def orientation(self,name,start,stop):
        if start <= stop:
            feat = [name,start,stop,1]
        else:
            feat = [name,stop,start,-1]
        return feat  
                                               
    def add_feature(self,seq_start,stop=None,name='new_feature',feature_type='misc_feature',fwd_color='',rev_color='',mismatch=0,overlap=True):
        '''add feature to sequence with tuple of coordinates, or with tuple of name and string to search sequence and add if present
        input: (1,100,"name"), (100,1,"name")-reverse complement, OR ('ATCG...',name = "name'") 
        first nucleotide coordinate is 1'''
        
        feat = ''
        if isinstance(seq_start,str) and seq_start.isalpha():
            loc = self.degenerate_search(seq_start,mismatch,overlap=overlap)
            if not loc:
                raise self.FeatureError(3)
                
        elif isinstance(seq_start,int) and isinstance(stop,int):
            if max(seq_start,stop) <= len(self.sequence) and min(seq_start,stop) > 0:
                loc = [(seq_start,stop)]          
            else:
                print(seq_start,stop,name)
                raise self.FeatureError(2)    
           
        else:
            raise self.FeatureError(1)
        
        fwd_color = self._color_check(fwd_color)
        
        if rev_color:
            rev_color = self._color_check(rev_color)
        else:
            rev_color = fwd_color
        feat_list =[]    
        for add_feat in loc: 
            feat = self.orientation(name,*add_feat)
            new_name=tuple(feat)
            feat_list.append(new_name)
            feat.extend([fwd_color,rev_color,feature_type])
            self.features[new_name]=feat
        return feat_list
             
    def _get_color(self):
        color = self._color_series[self._color_pos]
        self._color_pos +=1
        if self._color_pos>len(self._color_series)-1:
            self._color_pos =0
        return color  
         
    def _color_check(self,color):
        hex_color = r_compile("^#[0-9a-f]{6}$",flags=ig)
        if type(color) == str: 
            if color == '':
                pass
            elif hex_color.match(color):
                    return color
            else:
                dict_color = self.get_best_match(color,self.color_dict.keys(),verbose=False)
                if dict_color:
                    return self.color_dict[dict_color]
                else:
                    print('Color format is incorrect:',color)
                    print()
        else:
            print("Color format is incorrect!:",color)
        
        return self._get_color()
                    
    def add_comment(self,comment):
        if comment not in self.comments:
            self.comments.append(comment)
        else:
            print("Duplicate Comment: {}".format(comment))
            
    def open_ape(self,file_path,program_path=None):
        if program_path is not None:
            program_path=plp(program_path).resolve()
            try:
                if platform == "darwin":
                    run(["open","-a",program_path,file_path])
                elif platform == "win32":
                    run([str(program_path),str(file_path)])
                return
            except Exception as e:
                print(e)
                print("Could not open file {} with {}\nOpening with system default".format(file_path,str(program_path)))
        try:
            if platform == "darwin":
                run(["open",file_path])
            elif platform == "win32":
                run(["start","/wait",str(file_path)])
        except Exception as e:
            print(e)
            raise OSError("Could not open file {}, ensure ApE is installed".format(file_path))
                   
    def create_ape_file(self,ape_name='',path=None,do_open=False,program_path=None):
        if path is None:
            path = self._save_path
        path = plp(path)
        if not ape_name:
            ape_name = self.name
        self.file_name = ape_name+'.ape'.replace(".ape.ape",".ape")
        with open(path / self.file_name,'w') as out:
            comp1 = ['','','complement(']
            comp2 = ['','',')']
            out.write(self._ape_header1)
            for i in self.comments:
                out.write('COMMENT{}{}\n'.format(' '*5,i))
            out.write(self._ape_header2)
            for i in self.features.values():
                out.write(self._ape_feature.format([i[1],i[2]],i[0],i[4],i[5],' '*21,comp1[i[3]],comp2[i[3]],i[6][:15]+' '*(16-len(i[6][:15]))))
            out.write('ORIGIN\n{}\n//'.format(self.sequence))
        if do_open:
            file_path = os_path.join(path,self.file_name)
            self.open_ape(file_path,program_path)

def _parse_features(feat_list):
    features = od()
    for i in feat_list:
        key = (i[2],min(i[:2]),max(i[:2]),[-1,1][i[0]<i[1]])
        value = key+(i[4],i[5],i[3])
        features[key] = value
    return features
            
def parse_ape(ape_file):
    params = {
    "features":[],
    "comments":[],
    "ape_header1":[],
    "sequence":[]
    }
    
    with open(ape_file) as file_object:
        feature = False
        origin = False
        for line in file_object:
            line = line.strip()
            if origin:  
                params["sequence"].append(filter(lambda x: x.isalpha(),line))
            elif 'ORIGIN' in line:
                    feature = False
                    origin = True
            elif feature:
                feat_split = line.split('..')
                if len(feat_split) > 1:
                    r_complement = 'complement' in line
                    feat_start = feat_split[-2].split(' ')
                    feature_type = feat_start[0]
                    start = feat_start[-1]
                    end = feat_split[-1]
                    feat_start=[]
                    feat_end=[]
                    for i in [start,end][r_complement]:
                        if i.isdigit():
                            feat_start.append(i)
                    for i in [end,start][r_complement]:
                        if i.isdigit():
                            feat_end.append(i)
                    feat_start = ''.join(feat_start)
                    feat_end = ''.join(feat_end)
                elif '/label' in line:
                    line.split('=')
                    feat_label = line.split('=')[-1]
                elif 'fwdcolor=' in line:
                    fwd_color = line.split('=')[-1].replace('"','').replace("'",'')
                elif 'revcolor=' in line:
                    rev_color = line.split('=')[-1].replace('"','').replace("'",'')
                    params["features"].append((int(feat_start),int(feat_end),feat_label,feature_type,fwd_color,rev_color))                 
        
            elif 'COMMENT' in line:
                if 'ApEinfo:methylated' not in line:
                    params["comments"].append(line[len('COMMENT'):].strip())
                
            elif not feature and 'FEATURES' in line:
                feature = True
                
            elif not feature and not origin:
                params["ape_header1"].append(line)
    params["features"] = _parse_features(params["features"])            
    params["ape_header1"] = "\n".join(params["ape_header1"])
    params["sequence"] = "".join(params["sequence"])
    return params


def load_ape(**kwargs):
    ap = ape()
    [setattr(ap,*i) for i in kwargs.items()]
    return ap
              
def load_fasta(fa_file):
    fa = r_compile('>([^(\r\n|\r|\n)]*)[\r\n|\r|\n]([^>]*)')
    ap = []
    with open(fa_file,'rU') as sfile:
        seq = fa.findall(sfile.read())
    for i in seq:
        ap.append(ape(*i[::-1]))
    return ap
                
def load_file(ape_file, name=''):
    ape_file = plp(ape_file).resolve()
    err_message = 'Could not process {}, it\'s not an supported file'.format(ape_file)
    
    def check(ext):
        return ext.lower() in ['.fa','.fasta','.txt']
             
    ext = ape_file.suffix
    save_path = ape_file.parents[1]
    if not name:
        name = ape_file.parts[-1].replace(ext,"")
    if ext:
        if ext.lower() in ('.ape','.gb'):
            params = parse_ape(ape_file)
            if params.get("sequence",None):
                params.update({"name":name,"_save_path":save_path,"length":len(params["sequence"])})
                ap =  load_ape(**params)
            else:
                print(err_message)
                return None
        else:
            ap = load_fasta(ape_file)
            if not ap:
                params = parse_ape(ape_file)
                if params.get("sequence",None):
                    params.update({"name":name,"_save_path":save_path,"length":len(params["sequence"])})
                    ap =  load_ape(**params)                                    
                    if not ap.sequence:
                        ap = ape(ape_file,name,validate_sequence=True)
                        if not ap.sequence:
                            print(err_message)
                            return None
    try:             
        message = 'ApE file {} was loaded'.format(ap.name)
    except AttributeError:
        if len(ap)==1:
            ap = ap[0]
            message = 'ApE file {} was loaded'.format(ap.name)
        else:
            message = '{} ApE files were loaded'.format(len(ap))
    print(message)        
    return ap  
    
def test_color_file(colors,names=[]):
    from time import sleep
    if names:
        if len(colors)!=len(names):
            raise IndexError('names and colors lists must be same length')
    c = ape('','')
    seq = []
    for color in colors:
        segment = []
        for j in range(25):
            segment.append(random.choice('agct'))
        seq.append(''.join(segment))
        
    c.sequence = ''.join(seq) 
    for num,seg in enumerate(seq):
        color = colors[num]
        if names:
            name = names[num]
        else:
            name = color
        c.add_feature(seq_start=seg,name=name,fwd_color=color)
        c.add_comment('{} was added'.format(name))
        
    c.create_ape_file('whereami',do_open=True)
    sleep(2.5)
    remove(c.file_name)
    return seq                                                    
                                                
if __name__ =='__main__': 
    colors = list(ape().color_dict.values())
    hex_dict = dict([(i[1],i[0]) for i in ape().color_dict.items()])
    names = ["{}({})".format(hex_dict[i],i) for i in colors]
    sequence = test_color_file(colors,names)
   